function [rmap, centers] = compute_ratemap(x, y, rates, x_edges, y_edges)

    idx_x = discretize(x, x_edges);
    idx_y = discretize(y, y_edges);

    valid = ~isnan(idx_x) & ~isnan(idx_y);
    subs = [idx_y(valid), idx_x(valid)];
    val = rates(valid);

    sz = [length(y_edges) - 1, length(x_edges) - 1];
    
    count = accumarray(subs, val, sz, @sum, 0);
    occ   = accumarray(subs, 1,   sz, @sum, 0);

    rmap = count ./ occ;
    
    mask = occ > 0;
    rmap_filled = fillmissing(rmap, 'constant', 0);
    rmap_smooth = imgaussfilt(rmap_filled, 1, 'FilterSize', 5);
    
    rmap(~mask) = NaN;
    rmap(mask)  = rmap_smooth(mask);

    centers = {(x_edges(1:end-1) + x_edges(2:end)) / 2, ...
               (y_edges(1:end-1) + y_edges(2:end)) / 2};
end
