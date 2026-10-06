%% Platform head-direction cell screening
%
% Description:
%   Computes directional tuning, mean vector length (MVL) and split-half
%   tuning stability on the platform. Selects cells with MVL above the
%   95th percentile of 1,000 circular-shift shuffles and stability > 0.4,
%   using a 2 cm/s speed cutoff.
%
% Input:
%   Standardized Spilt_behave_calcium_data.mat, climb_position.mat,
%   head_direction.mat and NeuronActivity.mat.
%
% Output:
%   results/platform_hd_results.mat and results/platform_hd_screening.csv:
%   original cell IDs, tuning curves, MVL, preferred directions, stability
%   and shuffle thresholds.
%
% Author: Xuan Zhang, Xin Yuan / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
data_dir = fullfile(project_root, 'example_data');
bin_width = 3;
smoothing_bins = 2;
minimum_speed_cm_s = 2;
minimum_event_frames = 10;
stability_threshold = 0.4;
shuffle_percentile = 95;
minimum_shift_seconds = 30;
n_shuffle = 1000;
seed = 1;

addpath(genpath(functions_dir));
rng(seed);

loaded = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), 'Spilt_behave_calcium_data');
cell_ids = loaded.Spilt_behave_calcium_data.cell_filter_index(:);
clear loaded
assert(~isempty(cell_ids) && all(isfinite(cell_ids)) && all(cell_ids == fix(cell_ids)));
assert(numel(unique(cell_ids)) == numel(cell_ids));

position = load(fullfile(data_dir, 'climb_position.mat'), 'plane_position', 'plane_index');
heading = load(fullfile(data_dir, 'head_direction.mat'), 'plane_hd');
behavior_time = position.plane_position(:, 1);
assert(all(isfinite(behavior_time)) && all(diff(behavior_time) > 0));
assert(isequal(heading.plane_hd(:, 1), behavior_time));
plane_index = position.plane_index(:);
assert(~isempty(plane_index) && all(diff(plane_index) > 0));
assert(all(plane_index == fix(plane_index) & plane_index >= 1 & plane_index <= numel(behavior_time)));
is_plane = false(size(behavior_time));
is_plane(plane_index) = true;
plane_x = position.plane_position(:, 2);
plane_y = position.plane_position(:, 3);
speed = speed2D(plane_x, plane_y, behavior_time);
behavior_direction = mod(heading.plane_hd(:, 2), 360);
behavior_direction(~is_plane | ~isfinite(speed) | speed < minimum_speed_cm_s) = NaN;
behavior_dt = mean(diff(behavior_time));
plane_time_range = behavior_time(plane_index([1, end]))';
split_time = mean(plane_time_range);
first_direction = behavior_direction(behavior_time <= split_time);
second_direction = behavior_direction(behavior_time > split_time);
clear position heading plane_x plane_y speed

loaded = load(fullfile(data_dir, 'NeuronActivity.mat'), 'NeuronActivity');
calcium_time = loaded.NeuronActivity.time(:);
assert(all(isfinite(calcium_time)) && all(diff(calcium_time) > 0));
assert(size(loaded.NeuronActivity.Event_filtered_exp2, 1) == numel(calcium_time));
assert(all(cell_ids >= 1 & cell_ids <= size(loaded.NeuronActivity.Event_filtered_exp2, 2)));
assert(calcium_time(1) >= behavior_time(1) && calcium_time(end) <= behavior_time(end));
calcium_dt = mean(diff(calcium_time));
nearest_behavior = knnsearch(behavior_time, calcium_time);
plane_calcium_frame_index = find(is_plane(nearest_behavior));
plane_time = calcium_time(plane_calcium_frame_index);
event_direction = behavior_direction(nearest_behavior(plane_calcium_frame_index));
plane_events = loaded.NeuronActivity.Event_filtered_exp2(plane_calcium_frame_index, cell_ids) > 0;
clear loaded nearest_behavior
first_event = plane_time <= split_time;
second_event = plane_time > split_time;

n_frame = numel(plane_time);
sequence_time = (0:n_frame-1)' * calcium_dt;
time_30s = find(sequence_time >= minimum_shift_seconds, 1);
assert(~isempty(time_30s) && n_frame > 2 * time_30s, 'The platform recording is too short for the requested shuffle interval.');
occupancy_seconds = general.circHist(behavior_direction, bin_width)' * behavior_dt + eps;
n_cell = numel(cell_ids);
plane_tc = cell(n_cell, 1);
plane_tc_first = cell(n_cell, 1);
plane_tc_second = cell(n_cell, 1);
plane_MVL = nan(n_cell, 1);
plane_pfd = nan(n_cell, 1);
plane_stability = nan(n_cell, 1);
MVL_shuffle = nan(n_cell, n_shuffle);
shuffle_shift_frames = nan(n_cell, n_shuffle);
event_frame_count = sum(plane_events, 1)';
retained_event_frame_count = sum(plane_events & isfinite(event_direction), 1)';

for cell_index = 1:n_cell
    if event_frame_count(cell_index) < minimum_event_frames || retained_event_frame_count(cell_index) < minimum_event_frames
        continue
    end
    event = plane_events(:, cell_index);
    plane_tc{cell_index} = platform_tuning_curve(event_direction(event), behavior_direction, behavior_dt, bin_width, smoothing_bins);
    tc_stats = analyses.tcStatistics(plane_tc{cell_index}, bin_width, 50);
    plane_MVL(cell_index) = tc_stats.r;
    plane_pfd(cell_index) = tc_stats.peakDirection;

    spike_direction_first = event_direction(event & first_event);
    spike_direction_second = event_direction(event & second_event);
    if any(isfinite(spike_direction_first)) && any(isfinite(spike_direction_second))
        plane_tc_first{cell_index} = platform_tuning_curve(spike_direction_first, first_direction, behavior_dt, bin_width, smoothing_bins);
        plane_tc_second{cell_index} = platform_tuning_curve(spike_direction_second, second_direction, behavior_dt, bin_width, smoothing_bins);
        plane_stability(cell_index) = corr(plane_tc_first{cell_index}(:, 2), plane_tc_second{cell_index}(:, 2));
    end

    shifts = time_30s + randi(n_frame - 2 * time_30s, n_shuffle, 1);
    shuffle_shift_frames(cell_index, :) = shifts';
    shuffled_counts = zeros(numel(occupancy_seconds), n_shuffle);
    for shuffle_index = 1:n_shuffle
        shuffled_event = circshift(event, shifts(shuffle_index));
        shuffled_counts(:, shuffle_index) = general.circHist(event_direction(shuffled_event), bin_width)';
    end
    shuffled_rate = shuffled_counts ./ occupancy_seconds;
    if n_shuffle == 1
        shuffled_rate = general.circSmooth(shuffled_rate, smoothing_bins);
    else
        shuffled_rate = general.circSmooth(shuffled_rate, [smoothing_bins, 0]);
    end
    angle_radians = repmat(deg2rad(plane_tc{cell_index}(:, 1)), 1, n_shuffle);
    MVL_shuffle(cell_index, :) = circ_r(angle_radians, shuffled_rate);
    if mod(cell_index, 25) == 0 || cell_index == n_cell
        fprintf('Processed %d / %d cells\n', cell_index, n_cell);
    end
end

MVL_threshold = prctile(MVL_shuffle, shuffle_percentile, 2);
MVL_significant = plane_MVL > MVL_threshold;
is_platform_hd = MVL_significant & plane_stability > stability_threshold;
platform_hd_candidate_cells = cell_ids(MVL_significant);
platform_hd_cells = cell_ids(is_platform_hd);
screening_table = table(cell_ids, plane_MVL, MVL_threshold, plane_stability, plane_pfd, ...
    event_frame_count, retained_event_frame_count, is_platform_hd, ...
    'VariableNames', {'CellID', 'MVL', 'MVL_Shuffle95th', 'Stability', 'PFD_deg', 'EventFrames', 'RetainedEventFrames', 'IsPlatformHD'});
parameters = struct('bin_width_deg', bin_width, 'smoothing_bins', smoothing_bins, ...
    'minimum_speed_cm_s', minimum_speed_cm_s, 'minimum_event_frames', minimum_event_frames, ...
    'stability_threshold', stability_threshold, 'shuffle_percentile', shuffle_percentile, ...
    'minimum_shift_seconds', minimum_shift_seconds, 'minimum_shift_frames', time_30s, 'n_shuffle', n_shuffle, 'seed', seed, ...
    'heading_reference', 'platform_heading', 'behavior_sample_time_s', behavior_dt, 'calcium_sample_time_s', calcium_dt, ...
    'plane_time_range_s', plane_time_range, 'split_time_s', split_time, 'platform_calcium_frame_count', n_frame);
result_dir = fullfile(project_root, 'results');
if ~isfolder(result_dir)
    mkdir(result_dir);
end
save(fullfile(result_dir, 'platform_hd_results.mat'), ...
    'platform_hd_cells', 'platform_hd_candidate_cells', 'cell_ids', ...
    'plane_tc', 'plane_tc_first', 'plane_tc_second', 'plane_MVL', 'plane_pfd', 'plane_stability', ...
    'MVL_shuffle', 'MVL_threshold', 'shuffle_shift_frames', 'event_frame_count', 'retained_event_frame_count', ...
    'plane_calcium_frame_index', 'screening_table', 'parameters', '-v7.3');
writetable(screening_table, fullfile(result_dir, 'platform_hd_screening.csv'));
fprintf('Selected %d platform HD cells from %d analyzed cells.\n', numel(platform_hd_cells), n_cell);
disp(platform_hd_cells');

function tc = platform_tuning_curve(spike_direction, occupancy_direction, sample_time, bin_width, smoothing_bins)
spike_direction = spike_direction(isfinite(spike_direction));
tc = [];
if isempty(spike_direction) || ~any(isfinite(occupancy_direction))
    return
end
tc = analyses.turningCurve(spike_direction, occupancy_direction, sample_time, 'binWidth', bin_width, 'smooth', smoothing_bins);
end
