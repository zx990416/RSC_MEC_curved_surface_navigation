% Locate place fields in a firing map.
%
% Identifies the place fields in 2D firing map. First map is converted to binary image
% by applying a global threshold. During the second step connected regions are labeled
% by function bwlabel. Final step includes filtering of identified regions and verification
% that they are valid place fields.
%
%  USAGE
%   [fieldsMap, fields] = analyses.placefield_cylinder(map, <options>)
%   map           Firing rate map either structure obtained using <a href="matlab:help analyses.map">analyses.map</a>
%                 or a matrix that represents firing map.
%   <options>     optional list of property-value pairs (see table below)
%
%   =========================================================================
%    Properties     Values
%   -------------------------------------------------------------------------
%    'threshold'    value above threshold*peak belong to a field (default = 0.2).
%    'minBins'      Minimum number of bins in a place field, i.e. total number
%                   of bins of place field or an area of place field. Fields with
%                   fewer bins are not considered as place fields. Remember to
%                   adjust this value when you change the bin width.
%                   (default = 9).
%    'minPeak'      peaks smaller than this value are considered spurious
%                   and ignored (default = 1). Peak is normally a rate, however
%                   it's units not necessary are Hz.
%    'binWidth'     Width of the bins in cm. It is used to calculate field size.
%                   (default = 1).
%    'pos'          Position samples. Used to calculate posInd. If not provided,
%                   then posInd will be an empty matrix.
%    'minMeanRate'  Fields with mean rate smaller that this value are ignored (default = 0).
%   =========================================================================
%
%  OUTPUT
%
%   fieldsMap       Matrix of the same size as firing map. The elements of fieldsMap
%                   are integer values greater or equal to 0. The elements labeled 0
%                   are the background (not a place field). The pixels labeled 1 make up
%                   first field; the pixels labeled 2 make up a second object; and so on.
%
%   fields          Structure with information about each field. Structure fields are:
%       row         vector of rows that constitute this field;
%       col         vector of columns that constitute this field;
%       area        area of field measured by regionprops function, i.e. number of bins in the field;
%       bbox        bounding box of the field;
%       peak        field peak value;
%       size        field size in cm; Calculated with help of binWidth argument.
%       x, y        field centre of mass point;
%       meanRate    mean firing rate;
%       PixelIdxList    linear list of indices that constitute the field.
%       peakX       X-coordinate of field's peak
%       peakY       Y-coordinate of field's peak
%       map         Binary map of the field. It has the same size as firing map. The elements of 'map'
%                   equal to one if they belong to the field and are zeros otherwise.
%       posInd      Indices of position samples that correspond to this field. In case position sample
%                   matrix is of size Nx5 ([t x y x1 y1]), posInd corresponds to the left most position
%                   columns ([x y]). If pos argument is not provided, then posInd will be an empty matrix.
%
%  EXAMPLES
%
%   1. To obtain positions of the first field run:
%       [~, fields] = analyses.placefield(map, 'threshold', 0.3, 'minBins', 5, 'minPeak', 0.1, 'pos', pos);
%       fieldPos = pos(fields(1).posInd, :);
%

% Note to developers:
% In order to improve speed and capture all fields the map is processed in two passes. First
% we extract fields based solely on global peak value. This should be enough for most of the cases.
% Secondly we check if extracted fields maxima is the same as results of imregionalmax. If they are
% different, then we perform another extraction, which is slower.
%

% This script is used for smoothing the rate map of cylinder arena;
% x must be circular;
% Must be a matrix.
%
% Modified by Xiang Zhang, November, 2023.
% Modified 2026-10-04: periodic component deduplication and one-based field indices.

function [fieldsMap, fields] = placefield_cylinder(map, varargin)
    inp = inputParser;
    defaultMinPeak = 1.;
    defaultThreshold = 0.2;
    defaultMinBins = 9;
    defaultBinWidth = 1;
    defaultPos = [];
    defaultMinMeanRate = 0;

    % input argument check functions
    checkThreshold = @(x) helpers.isdscalar(x, '>=0', '<=1');
    checkScalarZero = @(x) helpers.isiscalar(x, '>=0');
    checkDScalar = @(x) helpers.isdscalar(x, '>0');
    checkDScalarZero = @(x) helpers.isdscalar(x, '>=0');

    % fill input parser object
    addRequired(inp, 'map');
    addParameter(inp, 'threshold', defaultThreshold, checkThreshold);
    addParameter(inp, 'minBins', defaultMinBins, checkScalarZero);
    addParameter(inp, 'minPeak', defaultMinPeak, checkDScalarZero);
    addParameter(inp, 'binWidth', defaultBinWidth, checkDScalar);
    addParameter(inp, 'minMeanRate', defaultMinMeanRate, checkScalarZero);
    addParameter(inp, 'pos', defaultPos, @(x) ismatrix(x) && size(x, 2) >= 3);

    parse(inp, map, varargin{:});

    % get parsed arguments
    minPeak = inp.Results.minPeak;
    threshold = inp.Results.threshold;
    minBins = inp.Results.minBins;
    binWidth = inp.Results.binWidth;
    pos = inp.Results.pos;
    minMeanRate = inp.Results.minMeanRate;

    originalMap = [];
    if isstruct(map)
        originalMap = map;
        map = map.z;
    end
    map = [map, map, map];

    fieldsMap = zeros(size(map));
    fields = struct('row', {}, 'col', {}, ...
        'size', {}, 'peak', {}, 'peakX', {}, 'peakY', {}, ...
        'area', {}, 'bbox', {}, ...
        'x', {}, 'y', {}, ...
        'meanRate', {}, 'PixelIdxList', {}, ...
        'map', {}, 'posInd', {} ...
        ,'MinorAxisLength',{},'MajorAxisLength',{},'Orientation',{});

    if isempty(map)
        return
    end

    globalPeak = nanmax(nanmax(map));
    if isnan(globalPeak) || globalPeak == 0
        fieldsMap = zeros(size(map,1), size(map,2)/3);
        return;
    end

    mapNans = isnan(map);
    map(mapNans) = 0;
    regionalMaxMap = imregionalmax(map, 4); % obtain all local maxima
    testRegionalMx = zeros(size(map)); % map that will contain located fields maxima

    [ir, ic] = find(regionalMaxMap > 0);            % get peak locations
    foundPeaks = map(sub2ind(size(map), ir ,ic));   % obtain peaks value

    % remove peaks that have smaller rate than minPeak
    selected = foundPeaks < minPeak;
    regionalMaxMap(sub2ind(size(map), ir(selected), ic(selected))) = 0;

    % Counter for the number of fields
    nFields = 0;

    binMap = ones(size(map));

    binMap(map < globalPeak*threshold) = 0;
    binMap(mapNans) = 0;

    cc = bwconncomp(binMap, 4);
    stats = regionprops(cc, 'Area', 'BoundingBox', 'Centroid','MinorAxisLength','MajorAxisLength','Orientation');

    for i = 1:cc.NumObjects
        linInd = cc.PixelIdxList{i};
        [fieldPeak, peakLinearInd] = max(map(linInd));
        [r, c] = ind2sub(size(map), linInd);
        meanRate = nanmean(map(linInd));

        [pr, pc] = ind2sub(size(map), linInd(peakLinearInd));
        % mark this field as visited even if we reject if further down
        testRegionalMx(pr, pc) = 1;
        
        if ~any(c > size(map,2)/3 & c <= size(map,2)*2/3)
            continue;
        end
        if fieldPeak < minPeak
            continue;
        end
        if meanRate < minMeanRate
            continue;
        end

        [r, c, wrappedPixels] = periodic_pixels(linInd, size(map,1), size(map,2)/3);
        if any(arrayfun(@(field) isequal(field.PixelIdxList, wrappedPixels), fields))
            continue
        end

        if length(r) >= minBins
            nFields = nFields + 1;

            fields(nFields).row = r;
            fields(nFields).col = c;
            fields(nFields).size = length(r) * binWidth^2;
            fields(nFields).peak = fieldPeak;
            fields(nFields).peakX = mod(pc - 1, size(map,2)/3) + 1;
            fields(nFields).peakY = pr;
            fields(nFields).area = length(r);
            fields(nFields).bbox = stats(i).BoundingBox;
            fields(nFields).PixelIdxList =  wrappedPixels;

            fields(nFields).x = mod(stats(i).Centroid(1) - size(map,2)/3, size(map,2)/3);
            fields(nFields).y = stats(i).Centroid(2);

            fields(nFields).meanRate = meanRate;
            fields(nFields).map = zeros(size(map));
            fields(nFields).map(linInd) = 1;
            fields(nFields).map(mapNans) = nan;
            fields(nFields).map = fields(nFields).map(:, 1:size(map,2)/3) + ...
                fields(nFields).map(:, (size(map,2)/3+1):2*size(map,2)/3) + ...
                fields(nFields).map(:, (2*size(map,2)/3+1):size(map,2));
            fields(nFields).map = double(fields(nFields).map > 0);
            fields(nFields).map(mapNans(:, 1:size(map,2)/3)) = NaN;

            % fields(nFields).Extent = stats(i).Extent;
            fields(nFields).Orientation = stats(i).Orientation;
            % fields(nFields).Eccentricity = stats(i).Eccentricity;
            fields(nFields).MinorAxisLength = stats(i).MinorAxisLength;
            fields(nFields).MajorAxisLength = stats(i).MajorAxisLength;
            % fields(nFields).Perimeter = stats(i).Perimeter;

            fields(nFields).posInd = field_position_rows(fields(nFields).map > 0, pos, originalMap);

            fieldsMap(linInd) = nFields;
        end
    end
    fieldsMap_temp = fieldsMap;
    fieldsMap_temp(:, (size(map,2)/3+1):2*size(map,2)/3) = fieldsMap(:, 1:size(map,2)/3) + ...
        fieldsMap(:, (size(map,2)/3+1):2*size(map,2)/3) + ...
        fieldsMap(:, (2*size(map,2)/3+1):size(map,2));

    if ~isequaln(testRegionalMx, regionalMaxMap)
        map(fieldsMap_temp > 0) = 0; % turn off map values for known fields, prevent field duplicates

        % we have some uncounted fields
        leftPeaksMap = regionalMaxMap - testRegionalMx;
        [ir, ic] = find(leftPeaksMap > 0);                    % get peak locations
        foundPeaks = map(sub2ind(size(leftPeaksMap), ir ,ic));   % obtain peaks value
        mapThresholds = foundPeaks * threshold;
        if sum(foundPeaks) == 0
            fieldsMap = field_labels(fields, [size(map,1), size(map,2)/3]);
            return;
        end

        finalMap = zeros(size(map));
        for i = 1:length(foundPeaks)
            binMap = zeros(size(map));
            binMap( (map > mapThresholds(i)) & (map <= foundPeaks(i)) ) = 1;
            binMap(mapNans) = 0;

            [binMap, linInd] = bwselect(binMap, ic(i), ir(i), 4); % leave only region that relates to current field
            if isempty(linInd)
                continue;
            end

            [r, ~] = periodic_pixels(linInd, size(map,1), size(map,2)/3);
            % check for minimum number of bins
            if length(r) < minBins
                continue;
            end

            stats = regionprops(binMap, 'Centroid', 'EulerNumber'); % find statistics on field candidates
            distToPeak = sqrt((stats.Centroid(1) - ic(i))^2 + (stats.Centroid(2) - ir(i))^2);
    %         distToPeak = pdist2(stats.Centroid, [ic(i) ir(i)]);
            if distToPeak > 4 || stats.EulerNumber < 1 % we want object without holes (euler number)
                continue;
            end
            finalMap(linInd) = 1;
        end

        cc = bwconncomp(finalMap, 4); % somehow regionprops is not working if finalMap is passed directly
        stats = regionprops(cc, 'Area', 'BoundingBox', 'Centroid','MinorAxisLength','MajorAxisLength','Orientation');

        for fieldInd = 1:length(stats)
            linInd = cc.PixelIdxList{fieldInd};
            [r, c] = ind2sub(size(map), linInd);
            meanRate = nanmean(map(linInd));
            [peakRate, peakInd] = nanmax(map(linInd));
            [pr, pc] = ind2sub(size(map), linInd(peakInd));
            
            if ~any(c > size(map,2)/3 & c <= size(map,2)*2/3)
                continue;
            end
            if peakRate < minPeak
                continue;
            end
            if meanRate < minMeanRate
                continue;
            end

            [r, c, wrappedPixels] = periodic_pixels(linInd, size(map,1), size(map,2)/3);
            if length(r) < minBins || any(arrayfun(@(field) isequal(field.PixelIdxList, wrappedPixels), fields))
                continue
            end

            nFields = nFields + 1;

            fields(nFields).row = r;
            fields(nFields).col = c;
            fields(nFields).size = length(r) * binWidth^2;
            fields(nFields).peak = peakRate;
            fields(nFields).peakX = mod(pc - 1, size(map,2)/3) + 1;
            fields(nFields).peakY = pr;
            fields(nFields).area = length(r);
            fields(nFields).bbox = stats(fieldInd).BoundingBox;
            fields(nFields).PixelIdxList = wrappedPixels;

            fields(nFields).x = mod(stats(fieldInd).Centroid(1) - size(map,2)/3, size(map,2)/3);
            fields(nFields).y = stats(fieldInd).Centroid(2);

            fields(nFields).meanRate = meanRate;
            fields(nFields).map = zeros(size(map));
            fields(nFields).map(linInd) = 1;
            fields(nFields).map(mapNans) = nan;
            fields(nFields).map = fields(nFields).map(:, 1:size(map,2)/3) + ...
                fields(nFields).map(:, (size(map,2)/3+1):2*size(map,2)/3) + ...
                fields(nFields).map(:, (2*size(map,2)/3+1):size(map,2));
            fields(nFields).map = double(fields(nFields).map > 0);
            fields(nFields).map(mapNans(:, 1:size(map,2)/3)) = NaN;

            fields(nFields).MinorAxisLength = stats(fieldInd).MinorAxisLength;
            fields(nFields).MajorAxisLength = stats(fieldInd).MajorAxisLength;
            fields(nFields).Orientation = stats(fieldInd).Orientation;

            fields(nFields).posInd = field_position_rows(fields(nFields).map > 0, pos, originalMap);

            fieldsMap(linInd) = nFields;
        end
    end
    fieldsMap = field_labels(fields, [size(map,1), size(map,2)/3]);
end


function [row, col, pixels] = periodic_pixels(indices, nrow, ncol)
[rows, columns] = ind2sub([nrow, 3*ncol], indices);
columns = mod(columns - 1, ncol) + 1;
pixels = unique(sub2ind([nrow, ncol], rows, columns));
[row, col] = ind2sub([nrow, ncol], pixels);
end

function rows = field_position_rows(mask, pos, originalMap)
rows = [];
if isempty(pos)
    return
end
if isempty(originalMap)
    xedges = linspace(nanmin(pos(:,2)), nanmax(pos(:,2)), size(mask,2)+1);
    yedges = linspace(nanmin(pos(:,3)), nanmax(pos(:,3)), size(mask,1)+1);
else
    xedges = originalMap.x;
    yedges = originalMap.y;
end
[~, ~, xbin] = histcounts(pos(:,2), xedges);
[~, ~, ybin] = histcounts(pos(:,3), yedges);
valid = find(xbin > 0 & ybin > 0);
inside = mask(sub2ind(size(mask), ybin(valid), xbin(valid)));
rows = valid(inside);
end

function labels = field_labels(fields, map_size)
labels = zeros(map_size);
for field = 1:numel(fields)
    labels(fields(field).PixelIdxList) = field;
end
end
