%% Dual-axis cell screening across the platform and cylindrical side
%
% Description:
%   Computes cross-surface MAE and tuning-curve correlation, and whole-session
%   split-half stability in the dual-axis reference frame. Selects cells with
%   MAE below the first percentile of 1,000 circular-shift shuffles and
%   stability > 0.3. Platform HD classification is not required.
%   A 2 cm/s speed cutoff is applied to platform tuning and whole-session
%   stability; the side MAE and correlation retain unfiltered side samples.
%
% Input:
%   Standardized Spilt_behave_calcium_data.mat, climb_position.mat,
%   head_direction.mat and NeuronActivity.mat.
%
% Output:
%   results/dual_axis_results.mat and results/dual_axis_screening.csv:
%   original cell IDs, tuning curves, MAE, correlation, stability and shuffles.
%
% Author: Xuan Zhang, Xin Yuan / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
data_dir = fullfile(project_root, 'example_data');
radius_cm = 11.5;
bin_width = 3;
smoothing_bins = 2;
minimum_speed_cm_s = 2;
stability_threshold = 0.3;
minimum_shift_seconds = 30;
n_shuffle = 1000;
seed = 1;

addpath(genpath(functions_dir));
rng(seed);

loaded = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), 'Spilt_behave_calcium_data');
data = loaded.Spilt_behave_calcium_data;
cell_ids = data.cell_filter_index(:);
assert(all(cell_ids == fix(cell_ids)) && numel(unique(cell_ids)) == numel(cell_ids));
assert(all(cell_ids >= 1 & cell_ids <= size(data.platform_calcium_event, 2)));
assert(size(data.platform_calcium_event, 2) == size(data.cylinder_calcium_event, 2));

plane_time = data.platform_calcium_time(:);
cylinder_time = data.cylinder_calcium_time(:);
assert(all(isfinite(plane_time)) && all(diff(plane_time) > 0));
assert(all(isfinite(cylinder_time)) && all(diff(cylinder_time) > 0));
assert(isequal(data.hd_dir_plane(:, 1), plane_time));
assert(isequal(data.behav_pos_plane(:, 1), plane_time));
assert(isequal(data.hd_cylinder(:, 1), cylinder_time));
assert(isequal(data.behav_pos_cylinder(:, 1), cylinder_time));
assert(size(data.platform_calcium_event, 1) == numel(plane_time));
assert(size(data.cylinder_calcium_event, 1) == numel(cylinder_time));
plane_dt = mean(diff(plane_time));
cylinder_dt = mean(diff(cylinder_time));
plane_events = data.platform_calcium_event(:, cell_ids) > 0;
cylinder_events = data.cylinder_calcium_event(:, cell_ids) > 0;

plane_direction = mod(data.hd_dir_plane(:, 2), 360);
plane_position = data.behav_pos_plane(:, 2:3);
plane_speed = speed2D(plane_position(:, 1), plane_position(:, 2), plane_time);
plane_direction(plane_speed < minimum_speed_cm_s) = NaN;
plane_position(plane_speed < minimum_speed_cm_s, :) = NaN;
plane_valid_rows = ~isnan(plane_direction) | any(~isnan(plane_position), 2);

cylinder_heading = mod(data.hd_cylinder(:, 3), 360);
cylinder_position_angle = mod(rad2deg(data.behav_pos_cylinder(:, 2) / radius_cm), 360);
cylinder_dualaxis_direction = mod(cylinder_heading + cylinder_position_angle + 90, 360);

position = load(fullfile(data_dir, 'climb_position.mat'), 'cylinder_position', 'plane_position', 'plane_index');
heading = load(fullfile(data_dir, 'head_direction.mat'), 'cylinder_hd', 'plane_hd');
whole_time = position.cylinder_position(:, 1);
assert(all(isfinite(whole_time)) && all(diff(whole_time) > 0));
assert(isequal(position.plane_position(:, 1), whole_time));
assert(isequal(heading.cylinder_hd(:, 1), whole_time));
assert(isequal(heading.plane_hd(:, 1), whole_time));
plane_index = position.plane_index(:);
assert(all(plane_index == fix(plane_index) & plane_index >= 1 & plane_index <= numel(whole_time)));
whole_direction = mod(heading.cylinder_hd(:, 2) + position.cylinder_position(:, 5) + 90, 360);
whole_direction(plane_index) = mod(heading.plane_hd(plane_index, 2), 360);
whole_cylinder_position = [deg2rad(position.cylinder_position(:, 5)) * radius_cm, position.cylinder_position(:, 4)];
whole_plane_position = position.plane_position(:, 2:3);
whole_speed = speed2D(whole_cylinder_position(:, 1), whole_cylinder_position(:, 2), whole_time);
whole_plane_speed = speed2D(whole_plane_position(:, 1), whole_plane_position(:, 2), whole_time);
whole_speed(plane_index) = whole_plane_speed(plane_index);
whole_speed_valid = isfinite(whole_speed) & whole_speed >= minimum_speed_cm_s;
whole_direction(~whole_speed_valid) = NaN;
whole_dt = mean(diff(whole_time));
whole_midpoint = (whole_time(1) + whole_time(end)) / 2;
whole_first_direction = whole_direction(whole_time <= whole_midpoint);
whole_second_direction = whole_direction(whole_time > whole_midpoint);
clear position heading whole_cylinder_position whole_plane_position whole_plane_speed

loaded = load(fullfile(data_dir, 'NeuronActivity.mat'), 'NeuronActivity');
whole_calcium_time = loaded.NeuronActivity.time(:);
assert(all(isfinite(whole_calcium_time)) && all(diff(whole_calcium_time) > 0));
assert(size(loaded.NeuronActivity.Event_filtered_exp2, 1) == numel(whole_calcium_time));
assert(size(loaded.NeuronActivity.Event_filtered_exp2, 2) == size(data.platform_calcium_event, 2));
assert(whole_calcium_time(1) >= whole_time(1) && whole_calcium_time(end) <= whole_time(end));
whole_events = loaded.NeuronActivity.Event_filtered_exp2(:, cell_ids) > 0;
clear loaded
whole_behavior_index = knnsearch(whole_time, whole_calcium_time);
whole_event_direction = whole_direction(whole_behavior_index);
whole_first_event = whole_calcium_time <= whole_midpoint;
whole_second_event = whole_calcium_time > whole_midpoint;

time_30s = find(cylinder_time >= cylinder_time(1) + minimum_shift_seconds, 1);
assert(~isempty(time_30s) && numel(cylinder_time) > 2 * time_30s, 'The side recording is too short for the requested shuffle interval.');
n_cell = numel(cell_ids);
cylinder_occupancy_seconds = general.circHist(cylinder_dualaxis_direction, bin_width)' * cylinder_dt + eps;
plane_tc = cell(n_cell, 1);
cylinder_tc = cell(n_cell, 1);
whole_tc_first = cell(n_cell, 1);
whole_tc_second = cell(n_cell, 1);
plane_pfd = nan(n_cell, 1);
dual_axis_stability = nan(n_cell, 1);
dual_axis_MAE = nan(n_cell, 1);
dual_axis_tc_correlation = nan(n_cell, 1);
MAE_shuffle = nan(n_cell, n_shuffle);
tc_correlation_shuffle = nan(n_cell, n_shuffle);
shuffle_shift_frames = nan(n_cell, n_shuffle);

for cell_index = 1:n_cell
    plane_event = plane_events(:, cell_index);
    cylinder_event = cylinder_events(:, cell_index);
    cylinder_tc{cell_index} = tuning_curve(cylinder_dualaxis_direction, cylinder_event, cylinder_dt, bin_width, smoothing_bins);

    whole_event = whole_events(:, cell_index);
    if sum(whole_event) >= 10 && sum(whole_event & isfinite(whole_event_direction)) >= 10
        [dual_axis_stability(cell_index), whole_tc_first{cell_index}, whole_tc_second{cell_index}] = ...
            whole_half_correlation(whole_first_direction, whole_second_direction, ...
            whole_event_direction(whole_event & whole_first_event), whole_event_direction(whole_event & whole_second_event), ...
            whole_dt, bin_width, smoothing_bins);
    end

    if sum(plane_event) >= 10 && sum(plane_event & plane_valid_rows) >= 10
        plane_tc{cell_index} = tuning_curve(plane_direction, plane_event, plane_dt, bin_width, smoothing_bins);
    end
    if isempty(plane_tc{cell_index})
        continue
    end
    plane_stats = analyses.tcStatistics(plane_tc{cell_index}, bin_width, 50);
    plane_pfd(cell_index) = plane_stats.peakDirection;
    predicted_heading = mod(plane_pfd(cell_index) - cylinder_position_angle + 270, 360);
    dual_axis_MAE(cell_index) = angular_mae(predicted_heading(cylinder_event), cylinder_heading(cylinder_event));
    if ~isempty(cylinder_tc{cell_index})
        dual_axis_tc_correlation(cell_index) = corr(plane_tc{cell_index}(:, 2), cylinder_tc{cell_index}(:, 2));
    end
    if sum(cylinder_event) < 10
        continue
    end

    shifts = time_30s + randi(numel(cylinder_time) - 2 * time_30s, n_shuffle, 1);
    shuffle_shift_frames(cell_index, :) = shifts';
    shuffled_event_counts = zeros(numel(cylinder_occupancy_seconds), n_shuffle);
    for shuffle_index = 1:n_shuffle
        shuffled_event = circshift(cylinder_event, shifts(shuffle_index));
        MAE_shuffle(cell_index, shuffle_index) = angular_mae(predicted_heading(shuffled_event), cylinder_heading(shuffled_event));
        shuffled_event_counts(:, shuffle_index) = general.circHist(cylinder_dualaxis_direction(shuffled_event), bin_width)';
    end
    shuffled_rate = shuffled_event_counts ./ cylinder_occupancy_seconds;
    if n_shuffle == 1
        shuffled_rate = general.circSmooth(shuffled_rate, smoothing_bins);
    else
        shuffled_rate = general.circSmooth(shuffled_rate, [smoothing_bins, 0]);
    end
    tc_correlation_shuffle(cell_index, :) = corr(plane_tc{cell_index}(:, 2), shuffled_rate);
    if mod(cell_index, 25) == 0 || cell_index == n_cell
        fprintf('Processed %d / %d cells\n', cell_index, n_cell);
    end
end

MAE_threshold = prctile(MAE_shuffle, 1, 2);
tc_correlation_threshold = prctile(tc_correlation_shuffle, 95, 2);
MAE_significant = dual_axis_MAE < MAE_threshold;
tc_significant = dual_axis_tc_correlation > tc_correlation_threshold;
is_dual_axis = MAE_significant & dual_axis_stability > stability_threshold;
dual_axis_candidate_cell = cell_ids(MAE_significant);
dual_axis_cell = cell_ids(is_dual_axis);
dual_axis_cell_tc_matched = cell_ids(is_dual_axis & tc_significant);

screening_table = table(cell_ids, dual_axis_MAE, dual_axis_tc_correlation, MAE_threshold, tc_correlation_threshold, dual_axis_stability, is_dual_axis, ...
    'VariableNames', {'CellID', 'MAE_deg', 'TCCorrelation', 'MAE_Shuffle1st', 'TCCorrelation_Shuffle95th', 'WholeSessionStability', 'IsDualAxis'});
parameters = struct('radius_cm', radius_cm, 'bin_width_deg', bin_width, 'smoothing_bins', smoothing_bins, ...
    'minimum_speed_cm_s', minimum_speed_cm_s, 'stability_threshold', stability_threshold, ...
    'whole_behavior_time_range_s', [whole_time(1), whole_time(end)], 'whole_calcium_time_range_s', [whole_calcium_time(1), whole_calcium_time(end)], ...
    'whole_behavior_sample_time_s', whole_dt, 'whole_split_time_s', whole_midpoint, ...
    'minimum_shift_seconds', minimum_shift_seconds, 'n_shuffle', n_shuffle, 'seed', seed);
result_dir = fullfile(project_root, 'results');
if ~isfolder(result_dir)
    mkdir(result_dir);
end
save(fullfile(result_dir, 'dual_axis_results.mat'), 'dual_axis_cell', 'dual_axis_candidate_cell', 'dual_axis_cell_tc_matched', ...
    'cell_ids', 'plane_tc', 'cylinder_tc', 'plane_pfd', 'dual_axis_MAE', 'dual_axis_tc_correlation', ...
    'MAE_shuffle', 'tc_correlation_shuffle', 'shuffle_shift_frames', 'MAE_threshold', 'tc_correlation_threshold', ...
    'dual_axis_stability', 'whole_tc_first', 'whole_tc_second', 'screening_table', 'parameters', '-v7.3');
writetable(screening_table, fullfile(result_dir, 'dual_axis_screening.csv'));
fprintf('Selected %d dual-axis cells from %d analyzed cells.\n', numel(dual_axis_cell), n_cell);
disp(dual_axis_cell');

function tc = tuning_curve(direction, event, sample_time, bin_width, smoothing_bins)
spike_direction = direction(event & isfinite(direction));
if isempty(spike_direction) || ~any(isfinite(direction))
    tc = [];
else
    tc = analyses.turningCurve(spike_direction, direction, sample_time, 'binWidth', bin_width, 'smooth', smoothing_bins);
end
end

function [correlation, first_tc, second_tc] = whole_half_correlation(first_direction, second_direction, first_spike_direction, second_spike_direction, sample_time, bin_width, smoothing_bins)
first_spike_direction = first_spike_direction(isfinite(first_spike_direction));
second_spike_direction = second_spike_direction(isfinite(second_spike_direction));
correlation = NaN;
first_tc = [];
second_tc = [];
if isempty(first_spike_direction) || isempty(second_spike_direction)
    return
end
first_tc = analyses.turningCurve(first_spike_direction, first_direction, sample_time, 'binWidth', bin_width, 'smooth', smoothing_bins);
second_tc = analyses.turningCurve(second_spike_direction, second_direction, sample_time, 'binWidth', bin_width, 'smooth', smoothing_bins);
correlation = corr(first_tc(:, 2), second_tc(:, 2));
end

function error_deg = angular_mae(predicted, observed)
valid = isfinite(predicted) & isfinite(observed);
error_deg = NaN;
if any(valid)
    error_deg = mean(abs(mod(predicted(valid) - observed(valid) + 180, 360) - 180));
end
end
