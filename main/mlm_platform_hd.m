%% Direction-position MLM for platform HD cells
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
screening_file = fullfile(project_root, 'results', 'platform_hd_results.mat');
output_dir = fullfile(project_root, 'results', 'MLM_platform_HD');
export_figures = true;

n_dir = 36;
n_pos_x = 4;
n_pos_y = 4;
max_iter = 1000;
tol = 1e-6;

loaded_data = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), ...
    'Spilt_behave_calcium_data');
decode_para = loaded_data.Spilt_behave_calcium_data;
required_fields = {'behav_pos_plane', 'hd_dir_plane', ...
    'platform_calcium_time', 'platform_calcium_event'};
assert(isstruct(decode_para) && isscalar(decode_para) && ...
    all(isfield(decode_para, required_fields)), ...
    'Spilt_behave_calcium_data is missing required platform-recording fields.');
assert(isfile(screening_file), ...
    'Run hd_platform_screening.m before running this analysis.');
selected = load(screening_file, 'platform_hd_cells', 'cell_ids', ...
    'plane_MVL', 'plane_stability', 'MVL_threshold');
cell_ids = selected.platform_hd_cells(:);
n_original_cell = size(decode_para.platform_calcium_event, 2);
assert(all(isfinite(cell_ids)) && all(cell_ids == fix(cell_ids)) && ...
    all(cell_ids >= 1 & cell_ids <= n_original_cell) && ...
    numel(unique(cell_ids)) == numel(cell_ids), ...
    'Selected HD cell IDs must be unique original calcium-event column indices.');
[found, screening_index] = ismember(cell_ids, selected.cell_ids(:));
assert(all(found), 'Selected cell IDs are missing from the original screening result.');
original_MVL = selected.plane_MVL(screening_index);
original_stability = selected.plane_stability(screening_index);
original_MVL_threshold = selected.MVL_threshold(screening_index);
assert(all(original_MVL > original_MVL_threshold & original_stability > 0.4), ...
    'Initial platform HD cells must pass MVL > shuffle 95th and stability > 0.4.');

time_vec = decode_para.platform_calcium_time(:);
n_frame = numel(time_vec);
assert(n_frame >= 2 && all(isfinite(time_vec)) && all(diff(time_vec) > 0));
assert(isequal(size(decode_para.behav_pos_plane), [n_frame, 3]) && ...
    isequal(size(decode_para.hd_dir_plane), [n_frame, 2]) && ...
    size(decode_para.platform_calcium_event, 1) == n_frame);
assert(all(abs(decode_para.behav_pos_plane(:, 1) - time_vec) < 1e-8) && ...
    all(abs(decode_para.hd_dir_plane(:, 1) - time_vec) < 1e-8), ...
    'Platform behavior, heading and calcium timestamps must align frame by frame.');
delta_t = mean(diff(time_vec));
data = struct('time', time_vec, 'x', decode_para.behav_pos_plane(:, 2), ...
    'y', decode_para.behav_pos_plane(:, 3), 'dir', decode_para.hd_dir_plane(:, 2), ...
    'spike', zeros(n_frame, 1));
idx_valid = isfinite(data.time) & isfinite(data.x) & isfinite(data.y) & isfinite(data.dir);
assert(sum(idx_valid) >= 2, 'Too few valid platform position and heading samples.');
geometry = mlm_platform_estimate(data, n_dir, n_pos_x, n_pos_y, delta_t, max_iter, tol);
assert(geometry.valid, geometry.note);
bin_def = struct('x_edges', geometry.position.x_edges, 'y_edges', geometry.position.y_edges);
x_centers = geometry.position.x_centers;
y_centers = geometry.position.y_centers;
t_mat = geometry.joint.time_M;
position_occupancy = geometry.position.occupancy2d;
dir_angles_deg = geometry.direction.tc(:, 1);

num_cells = numel(cell_ids);
direction_factor = nan(num_cells, n_dir);
position_factor = nan(num_cells, n_pos_x, n_pos_y);
raw_tc = cell(num_cells, 1);
direction_tc = cell(num_cells, 1);
raw_tc_statistics = cell(num_cells, 1);
corrected_tc_statistics = cell(num_cells, 1);
n_mat_all = zeros(num_cells, n_dir, n_pos_x, n_pos_y);
joint_rate = nan(num_cells, n_dir, n_pos_x, n_pos_y);
event_frame_count = sum(decode_para.platform_calcium_event(:, cell_ids) > 0, 1)';
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
    icell = cell_ids(iselect);
    data.spike = decode_para.platform_calcium_event(:, icell) > 0;
    retained_event_frame_count(iselect) = sum(data.spike(idx_valid));
    if retained_event_frame_count(iselect) == 0
        warning('No events with valid platform position and heading for cell %d. Skipping.', icell);
        continue;
    end
    fitted = mlm_platform_estimate(data, n_dir, n_pos_x, n_pos_y, ...
        delta_t, max_iter, tol, bin_def);
    assert(fitted.valid, fitted.note);
    assert(isequal(fitted.joint.time_M, t_mat));
    dir_occupancy = fitted.direction.occupancy;
    direction_count = sum(sum(fitted.joint.nspike_M, 3), 2);
    raw_rate = zeros(n_dir, 1);
    occupied = dir_occupancy > 0;
    raw_rate(occupied) = direction_count(occupied) ./ (dir_occupancy(occupied) * delta_t);
    raw_tc{iselect} = [dir_angles_deg, raw_rate];
    direction_tc{iselect} = fitted.direction.tc;
    tc_stat_0 = analyses.tcStatistics(raw_tc{iselect}, 360 / n_dir, 50);
    tc_stat_sep = analyses.tcStatistics(direction_tc{iselect}, 360 / n_dir, 50);

    direction_factor(iselect, :) = fitted.direction.rate';
    position_factor(iselect, :, :) = fitted.position.rate2d;
    raw_tc_statistics{iselect} = tc_stat_0;
    corrected_tc_statistics{iselect} = tc_stat_sep;
    n_mat_all(iselect, :, :, :) = fitted.joint.nspike_M;
    joint_rate(iselect, :, :, :) = fitted.joint.rate_M;
    raw_MVL(iselect) = tc_stat_0.r;
    corrected_MVL(iselect) = tc_stat_sep.r;
    raw_pfd(iselect) = tc_stat_0.peakDirection;
    corrected_pfd(iselect) = tc_stat_sep.peakDirection;

    if export_figures
        fig = plot_platform_mlm_tuning(raw_tc{iselect}, direction_tc{iselect}, ...
            fitted.position.rate2d, x_centers, y_centers, tc_stat_0, tc_stat_sep, icell);
        saveas(fig, fullfile(figure_dir, sprintf('cell_%d_tuning.png', icell)));
        close(fig);
    end
    if mod(iselect, 25) == 0 || iselect == num_cells
        fprintf('Processed %d / %d selected platform HD cells.\n', iselect, num_cells);
    end
end

is_mlm_platform_hd = corrected_MVL > original_MVL_threshold;
mlm_platform_hd_cells = cell_ids(is_mlm_platform_hd);
screening_table = table(cell_ids, event_frame_count, retained_event_frame_count, ...
    original_MVL, original_stability, original_MVL_threshold, ...
    raw_MVL, corrected_MVL, raw_pfd, corrected_pfd, is_mlm_platform_hd, ...
    'VariableNames', {'CellID', 'EventFrames', 'RetainedEventFrames', ...
    'OriginalMVL', 'OriginalStability', 'OriginalShuffle95th', ...
    'RawMVL', 'CorrectedMVL', 'RawPFD_deg', 'CorrectedPFD_deg', 'IsMLMPlatformHD'});
parameters = struct('direction_bins', n_dir, 'position_bins_x', n_pos_x, ...
    'position_bins_y', n_pos_y, 'max_iter', max_iter, 'tol', tol, ...
    'sample_time_s', delta_t, 'position_edge_percentiles', [1, 99], ...
    'position_reference', 'platform_xy', 'heading_reference', 'platform_heading', ...
    'speed_filter_applied', false, 'raw_curve_smoothing_bins', 0, ...
    'valid_sample_count', sum(idx_valid), 'original_cell_count', n_original_cell, ...
    'initial_stability_threshold', 0.4, 'shuffle_percentile', 95, ...
    'post_mlm_threshold', 'original_cellwise_MVL_shuffle_95th', ...
    'screening_file', screening_file);
save(fullfile(output_dir, 'mlm_platform_hd_results.mat'), 'cell_ids', ...
    'mlm_platform_hd_cells', 'is_mlm_platform_hd', ...
    'original_MVL', 'original_stability', 'original_MVL_threshold', ...
    'direction_factor', 'position_factor', 'raw_tc', 'direction_tc', ...
    'raw_tc_statistics', 'corrected_tc_statistics', 'n_mat_all', 't_mat', ...
    'joint_rate', 'position_occupancy', 'x_centers', 'y_centers', 'bin_def', ...
    'screening_table', 'parameters', '-v7.3');
writetable(screening_table, fullfile(output_dir, 'mlm_platform_hd_summary.csv'));
fprintf('Retained %d / %d initially selected platform HD cells after MLM.\n', ...
    numel(mlm_platform_hd_cells), num_cells);
disp(mlm_platform_hd_cells');
