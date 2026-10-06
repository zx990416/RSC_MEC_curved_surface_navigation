%% Zernike-Fourier GLM Encoding Model for Spatial & Head-Direction Cells
%
% Description:
%   Fits Poisson GLM using Zernike spatial basis and Fourier head direction
%   basis functions. Selects spatial order L using BIC and calculates MVL
%   from the fitted head-direction component with spatial terms held constant.
%
% Author: Xin Yuan, Xuan Zhang / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]

clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
addpath(functions_dir);
addpath(fullfile(functions_dir, 'externals'));

data_dir = fullfile(project_root, 'example_data');
save_dir = fullfile(project_root, 'results', 'Zernike_GLM');

loaded_data = load(fullfile(data_dir, 'Spilt_behave_calcium_data.mat'), ...
    'Spilt_behave_calcium_data');
decode_para_new_2508 = loaded_data.Spilt_behave_calcium_data;
required_fields = {'behav_pos_plane', 'hd_dir_plane', ...
    'platform_calcium_time', 'platform_calcium_event'};
assert(isstruct(decode_para_new_2508) && isscalar(decode_para_new_2508) && ...
    all(isfield(decode_para_new_2508, required_fields)), ...
    'Spilt_behave_calcium_data is missing required platform-recording fields.');

n_frame = size(decode_para_new_2508.platform_calcium_event, 1);
platform_time = decode_para_new_2508.platform_calcium_time(:);
assert(n_frame >= 2 && numel(platform_time) == n_frame && ...
    isequal(size(decode_para_new_2508.behav_pos_plane), [n_frame, 3]) && ...
    isequal(size(decode_para_new_2508.hd_dir_plane), [n_frame, 2]), ...
    'Platform positions, headings and calcium data must have matching frame counts.');
assert(all(isfinite(platform_time)) && all(diff(platform_time) > 0), ...
    'Platform calcium timestamps must be finite and strictly increasing.');
assert(all(abs(decode_para_new_2508.behav_pos_plane(:, 1) - platform_time) < 1e-8) && ...
    all(abs(decode_para_new_2508.hd_dir_plane(:, 1) - platform_time) < 1e-8), ...
    'Position and heading rows must align with platform calcium timestamps.');

required_functions = {'zernike_basis', 'compute_ratemap', 'compute_tuning', ...
    'get_sparsity', 'get_mvl', 'fitglm', 'imgaussfilt'};
missing_functions = required_functions(cellfun(@(name) isempty(which(name)), ...
    required_functions));
if ~isempty(missing_functions)
    error('ZernikeGLM:MissingDependencies', 'Missing analysis functions: %s.', ...
        strjoin(missing_functions, ', '));
end

if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

ncell = size(decode_para_new_2508.platform_calcium_event, 2);
num_grid = 23;        % Spatial rate map grid resolution (23x23)
num_dir_bins = 60;    % Direction bin count (6 degrees per bin)
dt = 0.02;            % Internal GLM bin width in seconds
J = 3;
nShuffle = 0;
hd_edges = linspace(-pi, pi, num_dir_bins + 1);
hd_centers = ((hd_edges(1:end-1) + hd_edges(2:end)) / 2)';
hd_centers_deg = rad2deg(hd_centers);
Hhd_curve = zeros(num_dir_bins, 2 * J);
for j = 1:J
    Hhd_curve(:, 2*j-1:2*j) = [sin(j * hd_centers), cos(j * hd_centers)];
end

p_hd_spatial = nan(ncell, 2);
L_kernel_spatial = nan(ncell, 1);
HD_tuning = nan(ncell, num_dir_bins);
HD_MVL = nan(ncell, 1);
Ori_MVL = nan(ncell, 1);
HD_coefficients = nan(ncell, 2 * J);
GLM_tuning  = nan(ncell, num_dir_bins);
Ori_tuning  = nan(ncell, num_dir_bins);
GLM_ratemap = nan(ncell, num_grid, num_grid);
Ori_ratemap = nan(ncell, num_grid, num_grid);

for icell = 1:ncell
    fprintf('\n================ Processing Cell %d/%d ================\n', icell, ncell);
    
    t_pos   = decode_para_new_2508.behav_pos_plane(:, 1);
    posx    = decode_para_new_2508.behav_pos_plane(:, 2);
    posy    = decode_para_new_2508.behav_pos_plane(:, 3);
    posdir  = decode_para_new_2508.hd_dir_plane(:, 2);
    
    spike_mask = decode_para_new_2508.platform_calcium_event(:, icell) > 0;
    t_spike    = t_pos(spike_mask);

    if isempty(t_spike)
        warning('No spikes detected for cell %d. Skipping.', icell);
        continue;
    end

    time_bins = (t_pos(1):dt:t_pos(end))';
    y_spikes  = histcounts(t_spike, [time_bins; time_bins(end) + dt])';

    x_raw  = interp1(t_pos, posx, time_bins, 'linear', 'extrap');
    y_raw  = interp1(t_pos, posy, time_bins, 'linear', 'extrap');
    hd_raw = deg2rad(interp1(t_pos, posdir, time_bins, 'linear', 'extrap'));

    x_center = (max(x_raw) + min(x_raw)) / 2;
    y_center = (max(y_raw) + min(y_raw)) / 2;
    x_norm   = x_raw - x_center;
    y_norm   = y_raw - y_center;

    max_r = max(sqrt(x_norm.^2 + y_norm.^2));
    [psi, rho] = cart2pol(x_norm / (max_r + 0.01), y_norm / (max_r + 0.01));

    Hhd = [];
    for j = 1:J
        Hhd = [Hhd, sin(j * hd_raw), cos(j * hd_raw)]; %#ok<AGROW>
    end

    L_range = 1:8;
    J_fixed = 5;
    bic_values = zeros(size(L_range));

    Hhd_bic = [];
    for j = 1:J_fixed
        Hhd_bic = [Hhd_bic, sin(j * hd_raw), cos(j * hd_raw)]; %#ok<AGROW>
    end

    for i = 1:length(L_range)
        L_curr = L_range(i);
        Hs_curr = [];
        for n = 0:L_curr
            for m = -n:2:n
                Hs_curr = [Hs_curr, zernike_basis(n, m, rho, psi)]; %#ok<AGROW>
            end
        end
        
        mdl_curr = fitglm([Hs_curr, Hhd_bic], y_spikes, 'Distribution', 'poisson', 'Link', 'log');
        k_params = mdl_curr.NumCoefficients;
        n_samples = mdl_curr.NumObservations;
        bic_values(i) = mdl_curr.Deviance + k_params * log(n_samples);
    end

    [~, best_idx] = min(bic_values);
    best_L = L_range(best_idx);
    L_kernel_spatial(icell) = best_L;
    fprintf('Optimal spatial order found: L = %d\n', best_L);

    Hs = [];
    for n = 0:best_L
        for m = -n:2:n
            Hs = [Hs, zernike_basis(n, m, rho, psi)]; %#ok<AGROW>
        end
    end

    X_final = [Hs, Hhd];
    mdl = fitglm(X_final, y_spikes, 'Distribution', 'poisson', 'Link', 'log');

    rate_pred = predict(mdl, X_final) / dt;
    rate_raw  = y_spikes / dt;
    rate_pred(isnan(rate_pred)) = 0;

    x_edges  = linspace(min(x_raw), max(x_raw), num_grid + 1);
    y_edges  = linspace(min(y_raw), max(y_raw), num_grid + 1);

    [raw_map, ~] = compute_ratemap(x_raw, y_raw, rate_raw, x_edges, y_edges);
    [glm_map, ~] = compute_ratemap(x_raw, y_raw, rate_pred, x_edges, y_edges);

    [raw_tuning, ~]          = compute_tuning(hd_raw, rate_raw, hd_edges);
    [glm_tuning, ~]          = compute_tuning(hd_raw, rate_pred, hd_edges);

    beta_hd = mdl.Coefficients.Estimate(end - 2*J + 1:end);
    log_hd_gain = Hhd_curve * beta_hd;
    hd_tuning = exp(log_hd_gain - max(log_hd_gain));
    HD_coefficients(icell, :) = beta_hd';
    HD_tuning(icell, :) = hd_tuning';
    HD_MVL(icell) = get_mvl(hd_tuning);
    Ori_MVL(icell) = get_mvl(raw_tuning(1:end-1));

    GLM_tuning(icell, :)   = glm_tuning(1:end-1);
    Ori_tuning(icell, :)   = raw_tuning(1:end-1);
    GLM_ratemap(icell,:,:) = glm_map;
    Ori_ratemap(icell,:,:) = raw_map;

    fprintf('Cell %d position-controlled HD MVL = %.4f\n', icell, HD_MVL(icell));

    if nShuffle > 0
        shuff_r_spatial = zeros(nShuffle, 1);
        shuff_r_hd      = zeros(nShuffle, 1);

        obs_sparsity = get_sparsity(raw_map);
        obs_mvl      = get_mvl(raw_tuning);

        min_shift = 10 * (1 / dt); % Minimum shift threshold (10 seconds)
        T = length(y_spikes);

        for s = 1:nShuffle
            shift_val = randi([min_shift, T - min_shift]);
            y_shuff   = circshift(y_spikes, shift_val);

            mdl_shuff = fitglm(X_final, y_shuff, 'Distribution', 'poisson');
            rate_shuff_pred = predict(mdl_shuff, X_final) / dt;

            [shuff_glm_map, ~]    = compute_ratemap(x_raw, y_raw, rate_shuff_pred, x_edges, y_edges);
            [shuff_glm_tuning, ~] = compute_tuning(hd_raw, rate_shuff_pred, hd_edges);

            shuff_r_spatial(s) = get_sparsity(shuff_glm_map(~isnan(shuff_glm_map)));
            shuff_r_hd(s)      = get_mvl(shuff_glm_tuning(~isnan(shuff_glm_tuning)));
        end

        p_spatial = sum(shuff_r_spatial >= obs_sparsity) / nShuffle;
        p_hd      = sum(shuff_r_hd >= obs_mvl) / nShuffle;

        p_hd_spatial(icell, :) = [p_hd, p_spatial];
        fprintf('Cell %d Significance -> Spatial p = %.4f, HD p = %.4f\n', icell, p_spatial, p_hd);
    end
end

save(fullfile(save_dir, 'GLM_top.mat'), ...
    'HD_tuning', 'HD_MVL', 'Ori_MVL', 'HD_coefficients', ...
    'hd_centers', 'hd_centers_deg', 'nShuffle', ...
    'p_hd_spatial', 'L_kernel_spatial', 'GLM_tuning', ...
    'Ori_tuning', 'GLM_ratemap', 'Ori_ratemap');

fprintf('Output saved to %s\n', fullfile(save_dir, 'GLM_top.mat'));
