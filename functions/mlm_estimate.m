function [dir_factor, pos_factor] = mlm_estimate(n, t, delta_t, max_iter, tol)

    if nargin < 4 || isempty(max_iter)
        max_iter = 1000;
    end
    if nargin < 5 || isempty(tol)
        tol = 1e-6;
    end

    [nx, ny] = size(n);

    pos_factor = ones(1, ny);
    dir_factor = ones(nx, 1);

    for iter = 1:max_iter
        dir_old = dir_factor;
        pos_old = pos_factor;

        denom_dir = t * pos_factor'; 
        dir_factor = sum(n, 2) ./ max(denom_dir, eps);
        dir_factor(denom_dir == 0) = 0;

        denom_pos = dir_factor' * t;
        pos_factor = sum(n, 1) ./ max(denom_pos, eps);
        pos_factor(denom_pos == 0) = 0;

        if sum(abs(dir_factor - dir_old)) < tol && sum(abs(pos_factor - pos_old)) < tol
            break;
        end
    end

    temp_p = dir_factor .* sum(t, 2) * delta_t;
    temp_d = pos_factor .* sum(t, 1) * delta_t;

    scale_factor_p = nansum(n(:)) / max(nansum(temp_p), eps);
    scale_factor_d = nansum(n(:)) / max(nansum(temp_d), eps);

    dir_factor = dir_factor * scale_factor_p;
    pos_factor = pos_factor * scale_factor_d;
end
