function [tuning, centers] = compute_tuning(hd, rates, edges)

    idx = discretize(hd, edges);
    valid = ~isnan(idx) & ~isnan(rates);

    sz = [length(edges) - 1, 1];
    tuning = accumarray(idx(valid), rates(valid), sz, @mean, 0);

    centers = (edges(1:end-1) + edges(2:end)) / 2;
    centers = centers(:);

    tuning = [tuning; tuning(1)];
    centers = [centers; centers(1)];
end
