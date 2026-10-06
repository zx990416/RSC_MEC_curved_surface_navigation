function fig = plot_tuning_curves(raw_tc, dir_factor, pos_factor, tc_stat_0, tc_stat_sep, cell_id)

    nx = length(dir_factor);
    ny = length(pos_factor);

    fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);

    subplot(1, 3, 1);
    raw_angles = deg2rad(raw_tc(:, 1));
    polarplot([raw_angles; raw_angles(1) + 2 * pi], [raw_tc(:, 2); raw_tc(1, 2)]);
    title({sprintf('Cell #%d Raw Tuning Curve', cell_id), ...
           sprintf('MVL = %.3f', tc_stat_0.r)});

    subplot(1, 3, 2);
    dir_angles = ((0:nx-1)' + 0.5) / nx * 2 * pi;
    polarplot([dir_angles; dir_angles(1) + 2 * pi], [dir_factor; dir_factor(1)]);
    title({sprintf('Direction Component'), ...
           sprintf('MVL = %.3f', tc_stat_sep.r)});

    subplot(1, 3, 3);
    pos_angles = ((0:ny-1) + 0.5) / ny * 2 * pi;
    polarplot([pos_angles, pos_angles(1) + 2 * pi], [pos_factor, pos_factor(1)]);
    title('Position Component');
end
