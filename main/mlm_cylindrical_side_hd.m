%% Direction-position MLM for cylindrical-side HD cells
%
% Description:
%   Decouples direction and position tuning components from neuronal activity
%   data using iterative MLE.
%
% Author: Xin Yuan, Xuan Zhang / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]
clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
addpath(genpath(functions_dir));
data_dir = fullfile(project_root, 'example_data');
screening_file = fullfile(project_root, 'results', 'cylindrical_side_hd_results.mat');
output_dir = fullfile(project_root, 'results', 'MLM_side_HD');
export_figures = true;

nx = 36;
ny = 20;
max_iter = 1000;
tol = 1e-6;

loaded_data = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), ...
    'Spilt_behave_calcium_data');
decode_para = loaded_data.Spilt_behave_calcium_data;
required_fields = {'behav_pos_cylinder', 'hd_cylinder', ...
    'cylinder_calcium_time', 'cylinder_calcium_event'};
assert(isstruct(decode_para) && isscalar(decode_para) && ...
    all(isfield(decode_para, required_fields)), ...
    'Spilt_behave_calcium_data is missing required side-recording fields.');
assert(isfile(screening_file), ...
    'Run hd_cylindrical_side_screening.m before running this analysis.');
selected = load(screening_file, 'cylindrical_side_hd_cells', 'cell_ids', ...
    'cylinder_MVL', 'cylinder_stability', 'MVL_threshold');
head_direction_cell = selected.cylindrical_side_hd_cells(:);
n_original_cell = size(decode_para.cylinder_calcium_event, 2);
assert(all(isfinite(head_direction_cell)) && ...
    all(head_direction_cell == fix(head_direction_cell)) && ...
    all(head_direction_cell >= 1 & head_direction_cell <= n_original_cell) && ...
    numel(unique(head_direction_cell)) == numel(head_direction_cell), ...
    'Selected HD cell IDs must be unique original calcium-event column indices.');

[found, screening_index] = ismember(head_direction_cell, selected.cell_ids(:));
assert(all(found), 'Selected cell IDs are missing from the original screening result.');
original_MVL = selected.cylinder_MVL(screening_index);
original_stability = selected.cylinder_stability(screening_index);
original_MVL_threshold = selected.MVL_threshold(screening_index);
assert(all(original_MVL > original_MVL_threshold & original_stability > 0.4), ...
    'Initial side HD cells must pass MVL > shuffle 95th and stability > 0.4.');

time_vec = decode_para.cylinder_calcium_time(:);
n_frame = numel(time_vec);
assert(n_frame >= 2 && all(isfinite(time_vec)) && all(diff(time_vec) > 0));
assert(isequal(size(decode_para.behav_pos_cylinder), [n_frame, 3]) && ...
    isequal(size(decode_para.hd_cylinder), [n_frame, 3]) && ...
    size(decode_para.cylinder_calcium_event, 1) == n_frame);
assert(all(abs(decode_para.behav_pos_cylinder(:, 1) - time_vec) < 1e-8) && ...
    all(abs(decode_para.hd_cylinder(:, 1) - time_vec) < 1e-8), ...
    'Side behavior, heading and calcium timestamps must align frame by frame.');
delta_t = nanmean(diff(time_vec));
behav_pos = [decode_para.behav_pos_cylinder(:, 1:2), decode_para.hd_cylinder(:, 3)];
idx_valid = all(isfinite(behav_pos), 2);
behav_pos = behav_pos(idx_valid, :);
assert(size(behav_pos, 1) >= 2, 'Too few valid side position and heading samples.');

position_percentiles = prctile(behav_pos(:, 2), [2, 98]);
x_range = position_percentiles(2) - position_percentiles(1);
assert(isfinite(x_range) && x_range > 0, 'The side circumference-position range must be positive.');
x_norm = mod(behav_pos(:, 3), 360) / 360;
y_norm = mod(behav_pos(:, 2), x_range) / x_range;
x_bins = linspace(0, 1, nx + 1);
y_bins = linspace(0, 1, ny + 1);
x_idx = discretize(x_norm, x_bins);
y_idx = discretize(y_norm, y_bins);
assert(all(isfinite(x_idx)) && all(isfinite(y_idx)));

dir_angles_deg = ((0:nx-1)' + 0.5) * 360 / nx;
position_centers = ((0:ny-1)' + 0.5) / ny;
t_mat = zeros(nx, ny);
for ix = 1:nx
    for iy = 1:ny
        t_mat(ix, iy) = sum(x_idx == ix & y_idx == iy);
    end
end

num_cells = numel(head_direction_cell);
cell_ids = head_direction_cell;
tc = cell(n_original_cell, 1);
direction_factor = nan(num_cells, nx);
position_factor = nan(num_cells, ny);
direction_tc = cell(num_cells, 1);
position_tc = cell(num_cells, 1);
raw_tc_statistics = cell(num_cells, 1);
corrected_tc_statistics = cell(num_cells, 1);
n_mat_all = zeros(num_cells, nx, ny);
joint_rate_per_frame = nan(num_cells, nx, ny);
event_frame_count = sum(decode_para.cylinder_calcium_event(:, cell_ids) > 0, 1)';
retained_event_frame_count = zeros(num_cells, 1);
raw_MVL = nan(num_cells, 1);
corrected_MVL = nan(num_cells, 1);
raw_pfd = nan(num_cells, 1);
corrected_pfd = nan(num_cells, 1);

if ~isfolder(output_dir)
    mkdir(output_dir);
end
if export_figures
    figure_dir = fullfile(output_dir, 'figures');
    if ~isfolder(figure_dir)
        mkdir(figure_dir);
    end
end

for iselect = 1:num_cells
    icell = head_direction_cell(iselect);
    spike_events = decode_para.cylinder_calcium_event(:, icell);
    spike_times = time_vec(spike_events > 0);
    retained_event_frame_count(iselect) = sum(spike_events(idx_valid) > 0);
    if retained_event_frame_count(iselect) == 0
        warning('No events with valid side position and heading for cell %d. Skipping.', icell);
        continue;
    end

    Mxy = nan(nx, ny);
    n_mat = zeros(nx, ny);
    for ix = 1:nx
        for iy = 1:ny
            idx_time = (x_idx == ix & y_idx == iy);
            time_in_bin = behav_pos(idx_time, 1);
            n_time = numel(time_in_bin);
            if n_time > 0
                n_spike = sum(ismember(spike_times, time_in_bin));
                Mxy(ix, iy) = n_spike / n_time;
                n_mat(ix, iy) = n_spike;
            end
        end
    end

    [dir_factor, pos_factor] = mlm_estimate(n_mat, t_mat, delta_t, max_iter, tol);
    raw_rate = zeros(nx, 1);
    dir_occupancy = sum(t_mat, 2);
    occupied = dir_occupancy > 0;
    raw_rate(occupied) = sum(n_mat(occupied, :), 2) ./ (dir_occupancy(occupied) * delta_t);
    tc{icell} = [dir_angles_deg, raw_rate];
    tc_stat_0 = analyses.tcStatistics(tc{icell}, 360 / nx, 50);
    tc_stat_sep = analyses.tcStatistics([dir_angles_deg, dir_factor], 360 / nx, 50);

    direction_factor(iselect, :) = dir_factor';
    position_factor(iselect, :) = pos_factor;
    direction_tc{iselect} = [dir_angles_deg, dir_factor];
    position_tc{iselect} = [position_centers, pos_factor(:)];
    raw_tc_statistics{iselect} = tc_stat_0;
    corrected_tc_statistics{iselect} = tc_stat_sep;
    n_mat_all(iselect, :, :) = n_mat;
    joint_rate_per_frame(iselect, :, :) = Mxy;
    raw_MVL(iselect) = tc_stat_0.r;
    corrected_MVL(iselect) = tc_stat_sep.r;
    raw_pfd(iselect) = tc_stat_0.peakDirection;
    corrected_pfd(iselect) = tc_stat_sep.peakDirection;

    if export_figures
        fig = plot_tuning_curves(tc{icell}, dir_factor, pos_factor, tc_stat_0, tc_stat_sep, icell);
        saveas(fig, fullfile(figure_dir, sprintf('cell_%d_tuning.png', icell)));
        close(fig);
    end
    if mod(iselect, 25) == 0 || iselect == num_cells
        fprintf('Processed %d / %d selected side HD cells.\n', iselect, num_cells);
    end
end

raw_tc = tc(cell_ids);
is_mlm_side_hd = corrected_MVL > original_MVL_threshold;
mlm_side_hd_cells = cell_ids(is_mlm_side_hd);
screening_table = table(cell_ids, event_frame_count, retained_event_frame_count, ...
    original_MVL, original_stability, original_MVL_threshold, ...
    raw_MVL, corrected_MVL, raw_pfd, corrected_pfd, is_mlm_side_hd, ...
    'VariableNames', {'CellID', 'EventFrames', 'RetainedEventFrames', ...
    'OriginalMVL', 'OriginalStability', 'OriginalShuffle95th', ...
    'RawMVL', 'CorrectedMVL', 'RawPFD_deg', 'CorrectedPFD_deg', 'IsMLMSideHD'});
parameters = struct('direction_bins', nx, 'position_bins', ny, 'max_iter', max_iter, ...
    'tol', tol, 'sample_time_s', delta_t, 'position_period_cm', x_range, ...
    'position_period_percentiles', [2, 98], 'position_phase_origin_cm', 0, ...
    'position_reference', 'circumference_arc', 'heading_reference', 'tangent_plane', ...
    'speed_filter_applied', false, 'raw_curve_smoothing_bins', 0, ...
    'valid_sample_count', size(behav_pos, 1), 'original_cell_count', n_original_cell, ...
    'initial_stability_threshold', 0.4, 'shuffle_percentile', 95, ...
    'post_mlm_threshold', 'original_cellwise_MVL_shuffle_95th', ...
    'screening_file', screening_file);
save(fullfile(output_dir, 'mlm_side_hd_results.mat'), 'cell_ids', 'mlm_side_hd_cells', ...
    'is_mlm_side_hd', 'original_MVL', 'original_stability', 'original_MVL_threshold', ...
    'direction_factor', 'position_factor', 'raw_tc', 'direction_tc', 'position_tc', ...
    'raw_tc_statistics', 'corrected_tc_statistics', 'n_mat_all', 't_mat', ...
    'joint_rate_per_frame', 'screening_table', 'parameters', '-v7.3');
writetable(screening_table, fullfile(output_dir, 'mlm_side_hd_summary.csv'));
fprintf('Retained %d / %d initially selected side HD cells after MLM.\n', ...
    numel(mlm_side_hd_cells), num_cells);
disp(mlm_side_hd_cells');
