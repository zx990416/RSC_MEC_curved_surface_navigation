function res = mlm_platform_estimate(data_in, n_dir, n_pos_x, n_pos_y, deltat, maxIter, tol, bin_def)

    if nargin < 8
        bin_def = [];
    end

    res = struct;
    res.valid = false;
    res.note = '';
    res.direction = struct;
    res.position = struct;
    res.joint = struct;
    res.meta = struct;

    required_fields = {'time','x','y','dir','spike'};
    for k = 1:numel(required_fields)
        if ~isfield(data_in, required_fields{k})
            res.note = ['Missing field: ', required_fields{k}];
            return;
        end
    end

    time = data_in.time(:);
    x_raw = data_in.x(:);
    y_raw = data_in.y(:);
    dir_deg = data_in.dir(:);
    spike = double(data_in.spike(:) > 0);

    idx_valid = isfinite(time) & isfinite(x_raw) & isfinite(y_raw) & isfinite(dir_deg) & isfinite(spike);
    time = time(idx_valid);
    x_raw = x_raw(idx_valid);
    y_raw = y_raw(idx_valid);
    dir_deg = dir_deg(idx_valid);
    spike = spike(idx_valid);

    if numel(time) < 2
        res.note = 'Too few valid samples.';
        return;
    end

    if isempty(bin_def)
        bin_def = prepare_top2d_bin_def(x_raw, y_raw, n_pos_x, n_pos_y);
    end

    x_edges = bin_def.x_edges;
    y_edges = bin_def.y_edges;

    x_clip = min(max(x_raw, x_edges(1) + eps), x_edges(end) - eps);
    y_clip = min(max(y_raw, y_edges(1) + eps), y_edges(end) - eps);

    [~,~,~,x_idx,y_idx] = histcounts2(x_clip, y_clip, x_edges, y_edges);

    dir_norm = mod(dir_deg, 360) / 360;
    dir_idx = floor(dir_norm * n_dir) + 1;
    dir_idx(dir_idx < 1) = 1;
    dir_idx(dir_idx > n_dir) = n_dir;

    good = x_idx > 0 & y_idx > 0 & dir_idx > 0;
    x_idx = x_idx(good);
    y_idx = y_idx(good);
    dir_idx = dir_idx(good);
    spike = spike(good);

    if isempty(x_idx)
        res.note = 'No valid binned samples.';
        return;
    end

    subs = [dir_idx(:), x_idx(:), y_idx(:)];
    time_M = accumarray(subs, 1, [n_dir, n_pos_x, n_pos_y], @sum, 0);
    nspike_M = accumarray(subs, spike, [n_dir, n_pos_x, n_pos_y], @sum, 0);

    rate_M = nan(n_dir, n_pos_x, n_pos_y);
    mask = time_M > 0;
    rate_M(mask) = nspike_M(mask) ./ (time_M(mask) * deltat);

    total_spikes = sum(nspike_M(:));

    dir_V = ones(n_dir, 1);
    pos_M = ones(n_pos_x, n_pos_y);

    for iter = 1:maxIter
        dir_old = dir_V;
        pos_old = pos_M;

        for idir = 1:n_dir
            occ_slice = squeeze(time_M(idir,:,:));
            denom = sum(sum(pos_M .* occ_slice)) * deltat;
            if denom > 0
                dir_V(idir) = sum(sum(squeeze(nspike_M(idir,:,:)))) / denom;
            else
                dir_V(idir) = 0;
            end
        end

        for ix = 1:n_pos_x
            for iy = 1:n_pos_y
                occ_col = squeeze(time_M(:,ix,iy));
                denom = sum(dir_V .* occ_col) * deltat;
                if denom > 0
                    pos_M(ix,iy) = sum(squeeze(nspike_M(:,ix,iy))) / denom;
                else
                    pos_M(ix,iy) = 0;
                end
            end
        end

        if sum(abs(dir_V(:) - dir_old(:))) < tol && sum(abs(pos_M(:) - pos_old(:))) < tol
            break
        end
    end

    dir_occ = squeeze(sum(sum(time_M,3),2));
    pos_occ2d = squeeze(sum(time_M,1));

    pred_spikes_from_dir = sum(dir_V .* dir_occ * deltat);
    pred_spikes_from_pos = sum(pos_M(:) .* pos_occ2d(:) * deltat);

    if pred_spikes_from_dir > 0
        dir_V = dir_V * (total_spikes / pred_spikes_from_dir);
    end

    if pred_spikes_from_pos > 0
        pos_M = pos_M * (total_spikes / pred_spikes_from_pos);
    end

    dir_V(dir_V < 0) = 0;
    pos_M(pos_M < 0) = 0;

    dir_centers_deg = ((0:n_dir-1)' + 0.5) * (360 / n_dir);
    x_centers = (x_edges(1:end-1) + x_edges(2:end)) / 2;
    y_centers = (y_edges(1:end-1) + y_edges(2:end)) / 2;

    dir_occ_prob = dir_occ / max(sum(dir_occ), eps);
    pos_occ_prob2d = pos_occ2d / max(sum(pos_occ2d(:)), eps);

    res.valid = true;
    res.note = 'ok';

    res.direction.rate = dir_V(:);
    res.direction.tc = [dir_centers_deg(:), dir_V(:)];
    res.direction.occupancy = dir_occ(:);
    res.direction.occupancy_prob = dir_occ_prob(:);

    res.position.rate2d = pos_M;
    res.position.rate = pos_M(:);
    res.position.occupancy2d = pos_occ2d;
    res.position.occupancy_prob2d = pos_occ_prob2d;
    res.position.occupancy_prob = pos_occ_prob2d(:);
    res.position.x_edges = x_edges;
    res.position.y_edges = y_edges;
    res.position.x_centers = x_centers(:);
    res.position.y_centers = y_centers(:);

    res.joint.time_M = time_M;
    res.joint.nspike_M = nspike_M;
    res.joint.rate_M = rate_M;

    res.meta.total_spikes = total_spikes;
    res.meta.n_samples = numel(time);
end

function bin_def = prepare_top2d_bin_def(x_raw, y_raw, n_pos_x, n_pos_y)

    valid = isfinite(x_raw) & isfinite(y_raw);
    x_raw = x_raw(valid);
    y_raw = y_raw(valid);

    x_low = prctile(x_raw, 1);
    x_high = prctile(x_raw, 99);
    y_low = prctile(y_raw, 1);
    y_high = prctile(y_raw, 99);

    if x_high <= x_low
        x_low = min(x_raw);
        x_high = max(x_raw) + eps;
    end

    if y_high <= y_low
        y_low = min(y_raw);
        y_high = max(y_raw) + eps;
    end

    bin_def = struct;
    bin_def.x_edges = linspace(x_low, x_high, n_pos_x + 1);
    bin_def.y_edges = linspace(y_low, y_high, n_pos_y + 1);
end
