function [inputA, model_names, features_struct] = extract_glm_features(decode_para, smooth_bin)

    if nargin < 2, smooth_bin = 10; end

    head_dirs = decode_para.hd_cylinder(:, 1:2);
    ahv = analyses.angularHeadVelocity(head_dirs);
    ahv = smooth_data_sep_nan(ahv, smooth_bin);
    
    x = decode_para.behav_pos_cylinder(:, 2);
    y = decode_para.behav_pos_cylinder(:, 3);
    t = decode_para.behav_pos_cylinder(:, 1);

    v   = smooth_data_sep_nan(speed2D_cy(x, y, t), smooth_bin);
    v_c = smooth_data_sep_nan(speed2D_cy(x, ones(size(y)), t), smooth_bin);
    v_h = smooth_data_sep_nan(speed2D_cy(ones(size(x)), y, t), smooth_bin);

    thresh_v = 18;
    v_c(v_c > thresh_v) = NaN;
    v_h(v_h > thresh_v) = NaN;

    thresh_ahv = 80;
    ahv(abs(ahv) > thresh_ahv) = NaN;

    height_cy = y;
    circ_pos = decode_para.behav_pos_cylinder(:, 2);
    R = (max(circ_pos) - min(circ_pos)) / 2;
    circ_ang = circ_pos / R * 180;

    side_body = decode_para.hd_cylinder(:, 3);
    x_dual = -circ_ang - 90;
    side_dual = decode_para.hd_cylinder(:, 3) - x_dual;

    fb_body  = fourier_bases(side_body);
    fb_circ  = fourier_bases(ang_mod_360(-circ_ang - 90));
    fb_cross = fourier_bases_cross(side_dual, ang_mod_360(-circ_ang - 90));

    model_names = {'full', 'no_self_motion', 'no_circ_pos', 'no_tp', 'no_cross_term'};
    inputA = {
        [ahv, v_c, v_h, height_cy, fb_body, fb_circ, fb_cross]; ...
        [height_cy, fb_body, fb_circ, fb_cross];               ...
        [ahv, v_c, v_h, height_cy, fb_body];                   ...
        [ahv, v_c, v_h, height_cy, fb_circ];                   ...
        [ahv, v_c, v_h, height_cy, fb_body, fb_circ]
    };

    features_struct.side_body = side_body;
    features_struct.side_dual = side_dual;
    features_struct.height_cy = height_cy;
    features_struct.full_matrix = inputA{1};
end
