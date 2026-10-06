%% Cylindrical-side position-related cell screening
%
% Description:
%   Computes spatial maps with a periodic horizontal axis, information
%   content, place fields and split-half spatial tuning stability.
%   Selects cells with information content above the 95th percentile of
%   1,000 circular-shift shuffles and stability > 0.3, using a 2 cm/s cutoff.
%
% Input:
%   Standardized Spilt_behave_calcium_data.mat, climb_position.mat and
%   NeuronActivity.mat.
%
% Output:
%   results/cylindrical_side_position_results.mat and
%   results/cylindrical_side_position_screening.csv: original cell IDs,
%   spatial maps, place fields, information, stability and shuffle thresholds.
%
% Author: Xuan Zhang, Xin Yuan / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
data_dir = fullfile(project_root, 'example_data');
radius_cm = 11.5;
bin_width_cm = 2;
smoothing_bins = 2;
minimum_speed_cm_s = 2;
minimum_event_frames = 10;
stability_threshold = 0.3;
shuffle_percentile = 95;
minimum_shift_seconds = 30;
n_shuffle = 1000;
seed = 1;
field_threshold = 0.4;
field_min_peak = 0.2;
field_min_bins = 9;

addpath(genpath(functions_dir));
rng(seed);
loaded = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), 'Spilt_behave_calcium_data');
cell_ids = loaded.Spilt_behave_calcium_data.cell_filter_index(:);
clear loaded
assert(~isempty(cell_ids) && all(isfinite(cell_ids)) && all(cell_ids == fix(cell_ids)));
assert(numel(unique(cell_ids)) == numel(cell_ids));

loaded = load(fullfile(data_dir, 'climb_position.mat'), 'cylinder_position', 'cylinder_index');
behavior_time = loaded.cylinder_position(:, 1);
cylinder_index = loaded.cylinder_index(:);
assert(all(isfinite(behavior_time)) && all(diff(behavior_time) > 0));
assert(~isempty(cylinder_index) && all(diff(cylinder_index) > 0));
assert(all(cylinder_index == fix(cylinder_index) & cylinder_index >= 1 & cylinder_index <= numel(behavior_time)));
is_cylinder = false(size(behavior_time));
is_cylinder(cylinder_index) = true;
arc = deg2rad(loaded.cylinder_position(:, 5)) * radius_cm;
height_cm = loaded.cylinder_position(:, 4);
speed = speed2D(arc, height_cm, behavior_time);
circumference_cm = 2 * pi * radius_cm;
pos = [behavior_time, mod(arc + circumference_cm/2, circumference_cm) - circumference_cm/2, height_cm];
valid_position = is_cylinder & isfinite(speed) & speed >= minimum_speed_cm_s & all(isfinite(pos(:, 2:3)), 2);
pos(~valid_position, 2:3) = NaN;
side_height = height_cm(is_cylinder & isfinite(height_cm));
limits = [-circumference_cm/2, circumference_cm/2, min(side_height), max(side_height)];
cylinder_time_range = behavior_time(cylinder_index([1, end]))';
split_time = mean(cylinder_time_range);
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
cylinder_time = calcium_time(cylinder_calcium_frame_index);
events = loaded.NeuronActivity.Event_filtered_exp2(cylinder_calcium_frame_index, cell_ids) > 0;
clear loaded nearest_behavior
first_event = cylinder_time <= split_time;
second_event = cylinder_time > split_time;
map_options = {'datatime', 's', 'binWidth', bin_width_cm, 'smooth', smoothing_bins, ...
    'minTime', 0, 'maxGap', 0.3, 'blanks', 'off', 'limits', limits};
full_cache = prepare_map(pos, cylinder_time, map_options);
first_cache = prepare_map(pos(behavior_time <= split_time, :), cylinder_time(first_event), map_options);
second_cache = prepare_map(pos(behavior_time > split_time, :), cylinder_time(second_event), map_options);
map_size = size(full_cache.template.z);
vertical_smoother = general.smooth_cylinder(eye(map_size(1)), [smoothing_bins, 0]);
horizontal_smoother = general.smooth_cylinder(eye(map_size(2)), [0, smoothing_bins]);
occupancy_pdf = full_cache.template.time(:) / (nansum(nansum(full_cache.template.time)) + eps);

n_frame = numel(cylinder_time);
time_30s = find((0:n_frame-1)' * calcium_dt >= minimum_shift_seconds, 1);
assert(~isempty(time_30s) && n_frame > 2*time_30s, 'The side recording is too short for the requested shuffle interval.');
n_cell = numel(cell_ids);
cylinder_map = cell(n_cell, 1);
cylinder_map_first = cell(n_cell, 1);
cylinder_map_second = cell(n_cell, 1);
fields_map = cell(n_cell, 1);
place_fields = cell(n_cell, 1);
field_number = nan(n_cell, 1);
information_content = nan(n_cell, 1);
information_rate = nan(n_cell, 1);
spatial_stability = nan(n_cell, 1);
information_content_shuffle = nan(n_cell, n_shuffle);
shuffle_shift_frames = nan(n_cell, n_shuffle);
event_frame_count = sum(events, 1)';
retained_event_frame_count = sum(events & isfinite(full_cache.frame_bin), 1)';

for cell_index = 1:n_cell
    if event_frame_count(cell_index) < minimum_event_frames || retained_event_frame_count(cell_index) < minimum_event_frames
        continue
    end
    event = events(:, cell_index);
    cylinder_map{cell_index} = event_map(full_cache, event, smoothing_bins);
    information = analyses.mapStatsPDF(cylinder_map{cell_index});
    information_content(cell_index) = information.content;
    information_rate(cell_index) = information.rate;
    [fields_map{cell_index}, place_fields{cell_index}] = analyses.placefield_cylinder(cylinder_map{cell_index}, ...
        'pos', pos, 'binWidth', bin_width_cm, 'threshold', field_threshold, 'minPeak', field_min_peak, 'minBins', field_min_bins);
    field_number(cell_index) = numel(place_fields{cell_index});

    cylinder_map_first{cell_index} = event_map(first_cache, event(first_event), smoothing_bins);
    cylinder_map_second{cell_index} = event_map(second_cache, event(second_event), smoothing_bins);
    if ~isempty(cylinder_map_first{cell_index}) && ~isempty(cylinder_map_second{cell_index})
        spatial_stability(cell_index) = zeroLagCorrelation(cylinder_map_first{cell_index}.z, cylinder_map_second{cell_index}.z);
    end

    shifts = time_30s + randi(n_frame - 2*time_30s, n_shuffle, 1);
    shuffle_shift_frames(cell_index, :) = shifts';
    shuffled_counts = zeros(prod(map_size), n_shuffle);
    for shuffle_index = 1:n_shuffle
        shuffled_counts(:, shuffle_index) = event_counts(full_cache, circshift(event, shifts(shuffle_index)));
    end
    shuffled_counts(full_cache.template.timeRaw(:) < 0.1, :) = 0;
    shuffled_rate = shuffled_counts ./ (full_cache.template.timeRaw(:) + eps);
    shuffled_rate = smooth_maps(shuffled_rate, map_size, vertical_smoother, horizontal_smoother);
    mean_rate = nansum(shuffled_rate .* occupancy_pdf, 1);
    rate_ratio = shuffled_rate ./ mean_rate;
    rate_ratio(rate_ratio < 1) = 1;
    information_content_shuffle(cell_index, :) = nansum(occupancy_pdf .* shuffled_rate .* log2(rate_ratio), 1) ./ mean_rate;
    if mod(cell_index, 25) == 0 || cell_index == n_cell
        fprintf('Processed %d / %d cells\n', cell_index, n_cell);
    end
end

information_threshold = prctile(information_content_shuffle, shuffle_percentile, 2);
information_significant = information_content > information_threshold;
is_position_related = information_significant & spatial_stability > stability_threshold;
position_related_candidate_cells = cell_ids(information_significant);
position_related_cells = cell_ids(is_position_related);
screening_table = table(cell_ids, information_content, information_threshold, information_rate, spatial_stability, field_number, ...
    event_frame_count, retained_event_frame_count, is_position_related, ...
    'VariableNames', {'CellID', 'InformationContent', 'Information_Shuffle95th', 'InformationRate', 'SpatialStability', ...
    'FieldNumber', 'EventFrames', 'RetainedEventFrames', 'IsPositionRelated'});
parameters = struct('radius_cm', radius_cm, 'circumference_cm', circumference_cm, 'bin_width_cm', bin_width_cm, ...
    'smoothing_bins', smoothing_bins, 'minimum_speed_cm_s', minimum_speed_cm_s, 'minimum_event_frames', minimum_event_frames, ...
    'stability_threshold', stability_threshold, 'shuffle_percentile', shuffle_percentile, 'n_shuffle', n_shuffle, 'seed', seed, ...
    'minimum_shift_seconds', minimum_shift_seconds, 'minimum_shift_frames', time_30s, 'max_gap_s', 0.3, 'minimum_bin_time_s', 0.1, ...
    'map_limits', limits, 'cylinder_time_range_s', cylinder_time_range, 'split_time_s', split_time, 'calcium_sample_time_s', calcium_dt, ...
    'field_threshold', field_threshold, 'field_min_peak', field_min_peak, 'field_min_bins', field_min_bins, 'horizontal_periodic', true);
result_dir = fullfile(project_root, 'results');
if ~isfolder(result_dir)
    mkdir(result_dir);
end
save(fullfile(result_dir, 'cylindrical_side_position_results.mat'), 'position_related_cells', 'position_related_candidate_cells', 'cell_ids', ...
    'information_content', 'information_rate', 'spatial_stability', 'information_content_shuffle', 'information_threshold', ...
    'cylinder_map', 'cylinder_map_first', 'cylinder_map_second', 'fields_map', 'place_fields', 'field_number', ...
    'shuffle_shift_frames', 'event_frame_count', 'retained_event_frame_count', 'cylinder_calcium_frame_index', 'screening_table', 'parameters', '-v7.3');
writetable(screening_table, fullfile(result_dir, 'cylindrical_side_position_screening.csv'));
fprintf('Selected %d cylindrical-side position-related cells from %d analyzed cells.\n', numel(position_related_cells), n_cell);
disp(position_related_cells');

function cache = prepare_map(pos, event_time, options)
valid = all(isfinite(pos(:, 2:3)), 2);
assert(any(valid), 'No valid positions in a recording half.');
cache.template = analyses.map_cylinder(pos, pos(find(valid, 1), 1), options{:});
cache.template.p = rmfield(cache.template.p, {'v', 'z', 'posGroups', 'spkGroups'});
[~, ~, xbin] = histcounts(pos(:, 2), cache.template.x);
[~, ~, ybin] = histcounts(pos(:, 3), cache.template.y);
cache.behavior_row = knnsearch(pos(:, 1), event_time);
valid = xbin(cache.behavior_row) > 0 & ybin(cache.behavior_row) > 0;
cache.frame_bin = nan(size(event_time));
cache.frame_bin(valid) = sub2ind(size(cache.template.z), ybin(cache.behavior_row(valid)), xbin(cache.behavior_row(valid)));
end

function counts = event_counts(cache, event)
selected = find(event & isfinite(cache.frame_bin));
[~, unique_rows] = unique(cache.behavior_row(selected), 'stable');
counts = accumarray(cache.frame_bin(selected(unique_rows)), 1, [numel(cache.template.z), 1]);
end

function map = event_map(cache, event, smoothing_bins)
map = [];
if ~any(event & isfinite(cache.frame_bin))
    return
end
map = cache.template;
map.countRaw = reshape(event_counts(cache, event), size(map.z));
map.Nspikes = map.countRaw;
map.Nspikes(map.timeRaw < 0.1) = 0;
map.zRaw = map.Nspikes ./ (map.timeRaw + eps);
map.z = general.smooth_cylinder(map.zRaw, smoothing_bins);
map.count = general.smooth_cylinder(map.countRaw, smoothing_bins);
map.peakRate = nanmax(nanmax(map.z));
map.meanRate = nansum(nansum(map.Nspikes)) / (nansum(nansum(map.timeRaw)) + eps);
end

function smoothed = smooth_maps(rate, map_size, vertical_smoother, horizontal_smoother)
n_map = size(rate, 2);
vertical = vertical_smoother * reshape(rate, map_size(1), []);
horizontal_input = reshape(permute(reshape(vertical, [map_size, n_map]), [2, 1, 3]), map_size(2), []);
horizontal = horizontal_smoother' * horizontal_input;
smoothed = reshape(permute(reshape(horizontal, [map_size(2), map_size(1), n_map]), [2, 1, 3]), prod(map_size), n_map);
end
