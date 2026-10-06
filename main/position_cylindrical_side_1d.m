%% One-dimensional spatial tuning of cylindrical-side position-related cells
%
% Description:
%   Computes horizontal circular and vertical height maps for previously
%   selected position-related cells, with a 2 cm/s speed cutoff. Calculates
%   information content, split-half stability, information z-scores from
%   1,000 circular-shift shuffles and horizontal circular STD.
%   No additional cell selection is performed.
%
% Input:
%   Standardized Spilt_behave_calcium_data.mat, climb_position.mat and
%   NeuronActivity.mat, plus results/cylindrical_side_position_results.mat
%   from position_cylindrical_side_screening.m.
%
% Output:
%   results/position_1d/cylindrical_side_position_1d_results.mat and
%   results/position_1d/cylindrical_side_position_1d_summary.csv: original
%   cell IDs, maps, bin centers, information, stability, z-scores and
%   horizontal circular STD in radians.
%
% Author: Xuan Zhang, Xin Yuan / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
data_dir = fullfile(project_root, 'example_data');
screening_file = fullfile(project_root, 'results', 'cylindrical_side_position_results.mat');
output_dir = fullfile(project_root, 'results', 'position_1d');
radius_cm = 11.5;
circular_bin_width_cm = 2*pi*radius_cm/60;
height_bin_width_cm = 1;
circular_limits_cm = [-pi*radius_cm, pi*radius_cm];
height_limits_cm = [0, 35];
smoothing_bins = 2;
minimum_speed_cm_s = 2;
minimum_shift_seconds = 30;
n_shuffle = 1000;
seed = 1;

addpath(genpath(functions_dir));
rng(seed);
assert(isfile(screening_file), 'Run position_cylindrical_side_screening.m first.');
selected = load(screening_file, 'position_related_cells', 'parameters', 'cylinder_calcium_frame_index');
cell_ids = selected.position_related_cells(:);
loaded = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), 'Spilt_behave_calcium_data');
cell_filter_index = loaded.Spilt_behave_calcium_data.cell_filter_index(:);
clear loaded
assert(all(ismember(cell_ids, cell_filter_index)));
assert(all(isfinite(cell_ids)) && all(cell_ids == fix(cell_ids)) && ...
    numel(unique(cell_ids)) == numel(cell_ids));

loaded = load(fullfile(data_dir, 'climb_position.mat'), 'cylinder_position', 'cylinder_index');
behavior_time = loaded.cylinder_position(:, 1);
cylinder_index = loaded.cylinder_index(:);
assert(all(isfinite(behavior_time)) && all(diff(behavior_time) > 0));
assert(~isempty(cylinder_index) && all(diff(cylinder_index) > 0));
assert(all(cylinder_index == fix(cylinder_index) & cylinder_index >= 1 & ...
    cylinder_index <= numel(behavior_time)));
is_cylinder = false(size(behavior_time));
is_cylinder(cylinder_index) = true;
arc = deg2rad(loaded.cylinder_position(:, 5)) * radius_cm;
height_cm = loaded.cylinder_position(:, 4);
speed = speed2D(arc, height_cm, behavior_time);
circumference_cm = 2*pi*radius_cm;
pos = [behavior_time, mod(arc + circumference_cm/2, circumference_cm) - circumference_cm/2, height_cm];
valid_position = is_cylinder & isfinite(speed) & speed >= minimum_speed_cm_s & ...
    all(isfinite(pos(:, 2:3)), 2);
pos(~valid_position, 2:3) = NaN;
cylinder_time_range = behavior_time(cylinder_index([1, end]))';
split_time = mean(cylinder_time_range);
assert(abs(split_time - selected.parameters.split_time_s) < 1e-8 && ...
    minimum_speed_cm_s == selected.parameters.minimum_speed_cm_s, ...
    'The position analysis and initial screening must use the same time split and speed cutoff.');
clear loaded arc height_cm speed

loaded = load(fullfile(data_dir, 'NeuronActivity.mat'), 'NeuronActivity');
calcium_time = loaded.NeuronActivity.time(:);
assert(all(isfinite(calcium_time)) && all(diff(calcium_time) > 0));
assert(size(loaded.NeuronActivity.Event_filtered_exp2, 1) == numel(calcium_time));
assert(all(cell_ids >= 1 & cell_ids <= size(loaded.NeuronActivity.Event_filtered_exp2, 2)));
assert(calcium_time(1) >= behavior_time(1) && calcium_time(end) <= behavior_time(end));
calcium_dt = mean(diff(calcium_time));
nearest_behavior = knnsearch(behavior_time, calcium_time);
cylinder_calcium_frame_index = find(is_cylinder(nearest_behavior));
assert(isequal(cylinder_calcium_frame_index, selected.cylinder_calcium_frame_index), ...
    'The initial screening result does not match the current side recording.');
cylinder_time = calcium_time(cylinder_calcium_frame_index);
events = loaded.NeuronActivity.Event_filtered_exp2(cylinder_calcium_frame_index, cell_ids) > 0;
clear loaded nearest_behavior selected
first_event = cylinder_time <= split_time;
second_event = cylinder_time > split_time;
first_behavior = behavior_time <= split_time;
second_behavior = behavior_time > split_time;

circular_options = {'datatime', 's', 'binWidth', circular_bin_width_cm, 'smooth', 0, ...
    'minTime', 0, 'maxGap', 0.3, 'blanks', 'off', 'limits', circular_limits_cm};
height_options = {'datatime', 's', 'binWidth', height_bin_width_cm, 'smooth', 0, ...
    'minTime', 0, 'maxGap', 0.3, 'blanks', 'off', 'limits', height_limits_cm};
circular_cache = prepare_map_1d(pos(:, 1:2), cylinder_time, circular_options, smoothing_bins, true);
height_cache = prepare_map_1d(pos(:, [1, 3]), cylinder_time, height_options, smoothing_bins, false);
circular_first_cache = prepare_map_1d(pos(first_behavior, 1:2), cylinder_time(first_event), ...
    circular_options, smoothing_bins, true);
circular_second_cache = prepare_map_1d(pos(second_behavior, 1:2), cylinder_time(second_event), ...
    circular_options, smoothing_bins, true);
height_first_cache = prepare_map_1d(pos(first_behavior, [1, 3]), cylinder_time(first_event), ...
    height_options, smoothing_bins, false);
height_second_cache = prepare_map_1d(pos(second_behavior, [1, 3]), cylinder_time(second_event), ...
    height_options, smoothing_bins, false);
assert(~isempty(circular_cache) && ~isempty(height_cache), 'No valid side spatial samples.');

n_frame = numel(cylinder_time);
minimum_shift_frames = find((0:n_frame-1)' * calcium_dt >= minimum_shift_seconds, 1);
assert(~isempty(minimum_shift_frames) && n_frame > 2*minimum_shift_frames, ...
    'The side recording is too short for the requested shuffle interval.');
shuffle_shift_frames = minimum_shift_frames + randi(n_frame - 2*minimum_shift_frames, n_shuffle, 1);
n_cell = numel(cell_ids);
circular_map = cell(n_cell, 1);
height_map = cell(n_cell, 1);
circular_map_first = cell(n_cell, 1);
circular_map_second = cell(n_cell, 1);
height_map_first = cell(n_cell, 1);
height_map_second = cell(n_cell, 1);
circular_information_content = nan(n_cell, 1);
height_information_content = nan(n_cell, 1);
circular_stability = nan(n_cell, 1);
height_stability = nan(n_cell, 1);
circular_std_rad = nan(n_cell, 1);
circular_resultant_length = nan(n_cell, 1);
circular_mean_rad = nan(n_cell, 1);
circular_information_content_shuffle = nan(n_cell, n_shuffle);
height_information_content_shuffle = nan(n_cell, n_shuffle);
event_frame_count = sum(events, 1)';
retained_event_frame_count = sum(events & isfinite(circular_cache.frame_bin) & ...
    isfinite(height_cache.frame_bin), 1)';

for cell_index = 1:n_cell
    event = events(:, cell_index);
    circular_map{cell_index} = event_map_1d(circular_cache, event);
    height_map{cell_index} = event_map_1d(height_cache, event);
    if ~isempty(circular_map{cell_index})
        information = analyses.mapStatsPDF(circular_map{cell_index});
        circular_information_content(cell_index) = information.content;
        [circular_std_rad(cell_index), circular_resultant_length(cell_index), ...
            circular_mean_rad(cell_index)] = ratemap_circ_std(circular_map{cell_index}.z);
    end
    if ~isempty(height_map{cell_index})
        information = analyses.mapStatsPDF(height_map{cell_index});
        height_information_content(cell_index) = information.content;
    end

    circular_map_first{cell_index} = event_map_1d(circular_first_cache, event(first_event));
    circular_map_second{cell_index} = event_map_1d(circular_second_cache, event(second_event));
    height_map_first{cell_index} = event_map_1d(height_first_cache, event(first_event));
    height_map_second{cell_index} = event_map_1d(height_second_cache, event(second_event));
    if ~isempty(circular_map_first{cell_index}) && ~isempty(circular_map_second{cell_index})
        circular_stability(cell_index) = zeroLagCorrelation( ...
            circular_map_first{cell_index}.z, circular_map_second{cell_index}.z);
    end
    if ~isempty(height_map_first{cell_index}) && ~isempty(height_map_second{cell_index})
        height_stability(cell_index) = zeroLagCorrelation( ...
            height_map_first{cell_index}.z, height_map_second{cell_index}.z);
    end

    circular_counts = zeros(numel(circular_cache.raw_time), n_shuffle);
    height_counts = zeros(numel(height_cache.raw_time), n_shuffle);
    for shuffle_index = 1:n_shuffle
        shuffled_event = circshift(event, shuffle_shift_frames(shuffle_index));
        circular_counts(:, shuffle_index) = event_counts_1d(circular_cache, shuffled_event);
        height_counts(:, shuffle_index) = event_counts_1d(height_cache, shuffled_event);
    end
    circular_information_content_shuffle(cell_index, :) = shuffle_information(circular_cache, circular_counts);
    height_information_content_shuffle(cell_index, :) = shuffle_information(height_cache, height_counts);
    if mod(cell_index, 10) == 0 || cell_index == n_cell
        fprintf('Processed %d / %d position-related cells.\n', cell_index, n_cell);
    end
end

circular_shuffle_mean = mean(circular_information_content_shuffle, 2, 'omitnan');
height_shuffle_mean = mean(height_information_content_shuffle, 2, 'omitnan');
circular_shuffle_std = std(circular_information_content_shuffle, 0, 2, 'omitnan');
height_shuffle_std = std(height_information_content_shuffle, 0, 2, 'omitnan');
circular_information_zscore = (circular_information_content - circular_shuffle_mean) ./ circular_shuffle_std;
height_information_zscore = (height_information_content - height_shuffle_mean) ./ height_shuffle_std;
circular_bin_centers_cm = (circular_cache.template.x(1:end-1) + circular_cache.template.x(2:end))/2;
height_bin_centers_cm = (height_cache.template.x(1:end-1) + height_cache.template.x(2:end))/2;
summary_table = table(cell_ids, circular_information_content, circular_stability, ...
    circular_information_zscore, circular_std_rad, height_information_content, ...
    height_stability, height_information_zscore, event_frame_count, retained_event_frame_count, ...
    'VariableNames', {'CellID', 'CircularInformationContent', 'CircularStability', ...
    'CircularInformationZscore', 'CircularSTD_rad', 'HeightInformationContent', ...
    'HeightStability', 'HeightInformationZscore', 'EventFrames', 'RetainedEventFrames'});
parameters = struct('radius_cm', radius_cm, 'circular_bin_width_cm', circular_bin_width_cm, ...
    'height_bin_width_cm', height_bin_width_cm, 'circular_limits_cm', circular_limits_cm, ...
    'height_limits_cm', height_limits_cm, 'smoothing_bins', smoothing_bins, ...
    'minimum_speed_cm_s', minimum_speed_cm_s, 'minimum_bin_time_s', 0.1, 'max_gap_s', 0.3, ...
    'blanks', 'on', 'horizontal_periodic', true, 'vertical_periodic', false, ...
    'n_shuffle', n_shuffle, 'seed', seed, 'minimum_shift_seconds', minimum_shift_seconds, ...
    'minimum_shift_frames', minimum_shift_frames, 'calcium_sample_time_s', calcium_dt, ...
    'split_time_s', split_time, 'cylinder_time_range_s', cylinder_time_range, ...
    'screening_file', screening_file, 'input_recording', 'full_synchronized_session', ...
    'zscore_std_normalization', 0, 'circular_std_unit', 'radian');
if ~isfolder(output_dir)
    mkdir(output_dir);
end
save(fullfile(output_dir, 'cylindrical_side_position_1d_results.mat'), 'cell_ids', ...
    'circular_map', 'height_map', 'circular_map_first', 'circular_map_second', ...
    'height_map_first', 'height_map_second', 'circular_bin_centers_cm', 'height_bin_centers_cm', ...
    'circular_information_content', 'height_information_content', ...
    'circular_stability', 'height_stability', 'circular_information_zscore', ...
    'height_information_zscore', 'circular_std_rad', 'circular_resultant_length', 'circular_mean_rad', ...
    'circular_information_content_shuffle', 'height_information_content_shuffle', ...
    'circular_shuffle_mean', 'height_shuffle_mean', 'circular_shuffle_std', 'height_shuffle_std', ...
    'shuffle_shift_frames', 'cylinder_calcium_frame_index', 'event_frame_count', ...
    'retained_event_frame_count', 'summary_table', 'parameters', '-v7.3');
writetable(summary_table, fullfile(output_dir, 'cylindrical_side_position_1d_summary.csv'));
fprintf('Saved horizontal/vertical maps and metrics for %d position-related cells.\n', n_cell);

function cache = prepare_map_1d(position, event_time, options, sigma, periodic)
cache = [];
valid = isfinite(position(:, 2));
if size(position, 1) < 2 || ~any(valid)
    return
end
if periodic
    template = analyses.map_cylinder(position, position(find(valid, 1), 1), options{:});
else
    template = analyses.map(position, position(find(valid, 1), 1), options{:});
end
n_bin = numel(template.z);
if periodic
    smoother = general.smooth_cylinder(eye(n_bin), [0, sigma])';
else
    smoother = general.smooth(eye(n_bin), [sigma, 0]);
end
raw_time = template.timeRaw(:);
template.time = (smoother * raw_time)';
template.time(raw_time == 0) = NaN;
template.timeRaw(raw_time == 0) = NaN;
template.timeRaw_original = raw_time';
template.p = rmfield(template.p, {'v', 'z', 'posGroups', 'spkGroups'});
template.p.smooth = sigma;
template.p.blanks = 'on';
coordinate = min(max(position(:, 2), template.p.limits(1)), template.p.limits(2));
coordinate(~isfinite(position(:, 2))) = NaN;
[~, ~, position_bin] = histcounts(coordinate, template.x);
behavior_row = knnsearch(position(:, 1), event_time);
frame_bin = position_bin(behavior_row);
frame_bin(frame_bin == 0) = NaN;
cache = struct('template', template, 'raw_time', raw_time, 'smoother', smoother, ...
    'behavior_row', behavior_row, 'frame_bin', frame_bin);
end

function counts = event_counts_1d(cache, event)
selected = find(event & isfinite(cache.frame_bin));
[~, unique_rows] = unique(cache.behavior_row(selected), 'stable');
counts = accumarray(cache.frame_bin(selected(unique_rows)), 1, [numel(cache.raw_time), 1]);
counts(cache.raw_time < 0.1) = 0;
end

function map = event_map_1d(cache, event)
map = [];
if isempty(cache) || ~any(event & isfinite(cache.frame_bin))
    return
end
map = cache.template;
counts = event_counts_1d(cache, event);
raw_rate = counts ./ (cache.raw_time + eps);
map.countRaw = counts';
map.Nspikes = counts';
map.zRaw = raw_rate';
map.z = (cache.smoother * raw_rate)';
map.count = (cache.smoother * counts)';
map.zRaw_original = map.zRaw;
map.countRaw_original = map.countRaw;
map.Nspikes_original = map.Nspikes;
blank = cache.raw_time == 0;
map.z(blank) = NaN;
map.zRaw(blank) = NaN;
map.countRaw(blank) = NaN;
map.Nspikes(blank) = NaN;
map.peakRate = max(map.z, [], 'omitnan');
map.meanRate = sum(counts) / (sum(cache.raw_time) + eps);
end

function content = shuffle_information(cache, counts)
rate = cache.smoother * (counts ./ (cache.raw_time + eps));
rate(cache.raw_time == 0, :) = NaN;
occupancy_pdf = cache.template.time(:) / (sum(cache.template.time, 'omitnan') + eps);
mean_rate = sum(rate .* occupancy_pdf, 1, 'omitnan');
rate_ratio = rate ./ mean_rate;
rate_ratio(rate_ratio < 1) = 1;
content = sum(occupancy_pdf .* rate .* log2(rate_ratio), 1, 'omitnan') ./ mean_rate;
end
