%% Rate-weighted circular dispersion of a one-dimensional spatial map
%
% Description:
%   Assigns equally spaced phases in [0, 2*pi) according to circular bin
%   order. Calculates the rate-weighted resultant length and mean phase,
%   and circular STD as sqrt(-2*log(R)).
%
% Input:
%   rate: Vector of spatial-bin rates in circular order. Only positive,
%   nonmissing rates contribute to the circular statistics.
%
% Output:
%   circStd: Circular standard deviation in radians.
%   R: Rate-weighted resultant length.
%   mu: Mean bin-order phase in [0, 2*pi), in radians.
%
% Author: Xuan Zhang, Xin Yuan / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

function [circStd, R, mu] = ratemap_circ_std(rate)
rate = rate(:);
N = numel(rate);
mask = ~isnan(rate) & rate > 0;
if ~any(mask)
    circStd = NaN;
    R = NaN;
    mu = NaN;
    warning('All bins are NaN or zero.');
    return;
end
w = rate(mask);
theta = linspace(0, 2*pi, N+1)';
theta(end) = [];
theta = theta(mask);
wsum = sum(w);
C = sum(w .* cos(theta)) / wsum;
S = sum(w .* sin(theta)) / wsum;
R = hypot(C, S);
mu = atan2(S, C);
if mu < 0
    mu = mu + 2*pi;
end
if R > 0
    circStd = sqrt(-2 * log(R));
else
    circStd = Inf;
end
end
