function sparsity = get_sparsity(v)
if sum(~isnan(v))==0
    sparsity=nan;
else
    meanrate = nansum(nansum(v));
    meansquarerate = nansum(nansum( (v .^ 2)));
    L = length(v(:));
    sparsity = meanrate^2 / meansquarerate;
    sparsity = (1-sparsity/L)*L/(L-1);
end
end
