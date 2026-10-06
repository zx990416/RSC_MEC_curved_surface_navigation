function fig = plot_platform_mlm_tuning(raw_tc, direction_tc, position_rate, ...
    x_centers, y_centers, tc_stat_0, tc_stat_sep, cell_id)
fig = figure('Visible', 'off', 'Position', [100, 100, 1200, 400]);
subplot(1, 3, 1);
angles = deg2rad(raw_tc(:, 1));
polarplot([angles; angles(1) + 2 * pi], [raw_tc(:, 2); raw_tc(1, 2)]);
title({sprintf('Cell #%d Raw Tuning Curve', cell_id), sprintf('MVL = %.3f', tc_stat_0.r)});
subplot(1, 3, 2);
angles = deg2rad(direction_tc(:, 1));
polarplot([angles; angles(1) + 2 * pi], [direction_tc(:, 2); direction_tc(1, 2)]);
title({'Direction Component', sprintf('MVL = %.3f', tc_stat_sep.r)});
subplot(1, 3, 3);
imagesc(x_centers, y_centers, position_rate');
set(gca, 'YDir', 'normal');
axis image;
colorbar;
xlabel('x (cm)');
ylabel('y (cm)');
title('Position Component');
end
