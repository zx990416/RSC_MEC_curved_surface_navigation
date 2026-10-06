%% GLM model fitting and model comparison

% Description:
%   Fits 5 GLM variants to neural data, calculates AIC for model selection,
%   and simulates synthetic spikes to reconstruct tuning curves.
%
% Author: Xin Yuan, Xuan Zhang / Miao Lab
% Repository: [https://github.com/zx990416/RSC_MEC_curved_surface_navigation]
clear; clc; close all;

project_root = fileparts(fileparts(mfilename('fullpath')));
functions_dir = fullfile(project_root, 'functions');
addpath(functions_dir);
addpath(fullfile(functions_dir, 'externals'));

data_dir = fullfile(project_root, 'example_data');
session_dirs = {data_dir};
session_names = {'example_data'};
target_save_dir = fullfile(project_root, 'results', 'AIC');

assert(numel(session_dirs) == numel(session_names), ...
    'Each input directory must have one session name.');

for i = 1:numel(session_dirs)
    middle_dir = session_names{i};
    para_path = session_dirs{i};
    fprintf('Processing: %s\n', middle_dir);

    try
        input_file = fullfile(para_path, 'Spilt_behave_calcium_data.mat');
        loaded_data = load(input_file, 'Spilt_behave_calcium_data');
        data = loaded_data.Spilt_behave_calcium_data;
        required_fields = {'behav_pos_cylinder', 'hd_cylinder', ...
            'cylinder_calcium_time', 'cylinder_calcium_event', 'cell_filter_index'};
        assert(isstruct(data) && isscalar(data) && all(isfield(data, required_fields)), ...
            'Spilt_behave_calcium_data is missing required side-recording fields.');

        decode_para = struct;
        decode_para.behav_pos_cylinder = data.behav_pos_cylinder;
        decode_para.hd_cylinder = data.hd_cylinder;
        decode_para.calcium_time = data.cylinder_calcium_time(:)';
        decode_para.calcium_event = data.cylinder_calcium_event;
        cell_filter_index = data.cell_filter_index;

        n_frame = size(decode_para.calcium_event, 1);
        assert(n_frame >= 2 && numel(decode_para.calcium_time) == n_frame && ...
            isequal(size(decode_para.behav_pos_cylinder), [n_frame, 3]) && ...
            isequal(size(decode_para.hd_cylinder), [n_frame, 3]), ...
            'Side positions, headings and calcium data must have matching frame counts.');
        assert(all(isfinite(decode_para.calcium_time)) && ...
            all(diff(decode_para.calcium_time) > 0), ...
            'Side calcium timestamps must be finite and strictly increasing.');
        assert(all(abs(decode_para.behav_pos_cylinder(:, 1) - ...
            decode_para.calcium_time(:)) < 1e-8) && ...
            all(abs(decode_para.hd_cylinder(:, 1) - ...
            decode_para.calcium_time(:)) < 1e-8), ...
            'Position and heading rows must align with side calcium timestamps.');

        required_functions = {'fourier_bases', 'fourier_bases_cross'};
        missing_functions = required_functions(cellfun(@(name) isempty(which(name)), ...
            required_functions));
        if ~isempty(missing_functions)
            error('AIC:MissingDependencies', 'Missing original model functions: %s.', ...
                strjoin(missing_functions, ', '));
        end

        smooth_path = which('smooth');
        if isempty(smooth_path)
            error('AIC:MissingCurveFittingToolbox', ...
                'MATLAB smooth is unavailable. Install or enable Curve Fitting Toolbox before running AIC.');
        end
        if ~startsWith(smooth_path, [matlabroot, filesep])
            error('AIC:ShadowedSmooth', ...
                'MATLAB smooth is shadowed by %s. Remove the custom function from the search path.', ...
                smooth_path);
        end

        if ~exist(target_save_dir, 'dir')
            mkdir(target_save_dir);
        end

        sampleTime = 1;
        AngleSmooth = 2;
        AngleBinsize = 3;
        ncell = size(decode_para.calcium_event, 2);

        [inputA, nameA, feat] = extract_glm_features(decode_para);

        AIC_M = nan(5, ncell);
        struct_dualM = cell(ncell, 6);
        struct_bodyM = cell(ncell, 6);
        tc_dual = nan(ncell, 120, 6);
        tc_body = nan(ncell, 120, 6);

        valid_mask = all(~isnan(feat.full_matrix), 2);

        for icell = 1:ncell
            cal_event = decode_para.calcium_event(:, icell) > 0;
            event_time = decode_para.calcium_time(cal_event);
            behav_time = decode_para.behav_pos_cylinder(:, 1);

            matched_idx = knnsearch(behav_time, event_time');
            is_spike = false(size(behav_time));
            is_spike(matched_idx) = true;

            spike_poi0 = find(is_spike(valid_mask));
            if length(spike_poi0) < 20
                continue;
            end

            AIC_all = nan(5, 1);
            
            for imodel = 1:5
                inputM = inputA{imodel};

                mdl = fitglm(inputM(valid_mask, :), is_spike(valid_mask), ...
                             'linear', 'Distribution', 'poisson');
                B_t = mdl.predict(inputM(valid_mask, :));
                
                E = nanmean(B_t);
                mu = length(spike_poi0) / length(B_t);
                k = mu / max(E, eps);

                spike_ori = find(is_spike(valid_mask));
                no_spike_ori = find(~is_spike(valid_mask));

                p_ori = B_t(spike_ori) * k;
                no_p_ori = 1 - B_t(no_spike_ori) * k;

                AIC_all(imodel) = -2 * sum(log([p_ori; no_p_ori]));

                syn_spikes = simulate_inhom_poisson(B_t, k);

                if length(syn_spikes(syn_spikes > 0.5)) < 10
                    continue;
                end

                side_body0 = feat.side_body(valid_mask);
                side_dual0 = feat.side_dual(valid_mask);

                syn_spikes_valid = round(syn_spikes(syn_spikes > 0.5));
                
                hdtc_dual = analyses.turningCurve(side_dual0(syn_spikes_valid), side_dual0, ...
                    sampleTime, 'smooth', AngleSmooth, 'binWidth', AngleBinsize);
                hdtc_body = analyses.turningCurve(side_body0(syn_spikes_valid), side_body0, ...
                    sampleTime, 'smooth', AngleSmooth, 'binWidth', AngleBinsize);

                struct_dualM{icell, imodel} = analyses.tcStatistics(hdtc_dual, AngleBinsize, 50);
                struct_bodyM{icell, imodel} = analyses.tcStatistics(hdtc_body, AngleBinsize, 50);

                tc_dual(icell, :, imodel) = hdtc_dual(:, 2);
                tc_body(icell, :, imodel) = hdtc_body(:, 2);
            end

            AIC_M(:, icell) = AIC_all;

            real_spikes = find(is_spike(valid_mask));
            side_body0 = feat.side_body(valid_mask);
            side_dual0 = feat.side_dual(valid_mask);

            hdtc_dual_real = analyses.turningCurve(side_dual0(real_spikes), side_dual0, ...
                sampleTime, 'smooth', AngleSmooth, 'binWidth', AngleBinsize);
            hdtc_body_real = analyses.turningCurve(side_body0(real_spikes), side_body0, ...
                sampleTime, 'smooth', AngleSmooth, 'binWidth', AngleBinsize);

            struct_dualM{icell, 6} = analyses.tcStatistics(hdtc_dual_real, AngleBinsize, 50);
            struct_bodyM{icell, 6} = analyses.tcStatistics(hdtc_body_real, AngleBinsize, 50);

            tc_dual(icell, :, 6) = hdtc_dual_real(:, 2);
            tc_body(icell, :, 6) = hdtc_body_real(:, 2);
        end

        save_file_name = fullfile(target_save_dir, ['AIC_', middle_dir, '.mat']);
        save(save_file_name, 'struct_dualM', 'struct_bodyM', 'tc_dual', 'tc_body', 'AIC_M');
    catch ME
        fprintf('Error processing session %s: %s\n', middle_dir, ME.message);
        rethrow(ME);
    end
end

disp('All configured sessions completed.');
