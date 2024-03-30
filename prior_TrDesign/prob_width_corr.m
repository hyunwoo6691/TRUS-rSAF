h = findobj(gca,'Type','line');

x = get(h, 'Xdata');
y = get(h, 'Ydata');

%%
rsaf_cw = y{2};
rsaf_pw = y{1};

rsaf_ = zeros(2, numel(rsaf_pw));
rsaf_(1,:) = rsaf_cw;
rsaf_(2,:) = rsaf_pw;
%%
figure; 
plot(rsaf_cw, rsaf_pw,'LineStyle', 'None' ,'LineWidth', 4, 'color', 'k', 'Marker', '^', 'MarkerSize', 10);
xlabel('rSAF-CW','Interpreter','None'); ylabel('rSAF-PW', 'Interpreter','None');

%%
p = polyfit(rsaf_cw, rsaf_pw, 1);

fit_curve = polyval(p, rsaf_cw);

figure(1);
% plot(rsaf_cw, rsaf_pw,'LineStyle', 'None' ,'LineWidth', 2, 'color', 'k', 'Marker', '^', 'MarkerSize', 10); hold on;
scatter(rsaf_cw(1), rsaf_pw(1), 'MarkerEdgeColor', [0.12 0.12 0.12]*1, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100); hold on;
scatter(rsaf_cw(2), rsaf_pw(2), 'MarkerEdgeColor', [0.12 0.12 0.12]*2, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100);
scatter(rsaf_cw(3), rsaf_pw(3), 'MarkerEdgeColor', [0.12 0.12 0.12]*3, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100);
scatter(rsaf_cw(4), rsaf_pw(4), 'MarkerEdgeColor', [0.12 0.12 0.12]*4, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100);
scatter(rsaf_cw(5), rsaf_pw(5), 'MarkerEdgeColor', [0.12 0.12 0.12]*5, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100);
scatter(rsaf_cw(6), rsaf_pw(6), 'MarkerEdgeColor', [0.12 0.12 0.12]*6, 'LineWidth', 2, 'Marker', 'x', 'SizeData', 100);
plot(rsaf_cw, fit_curve, 'LineWidth', 2, 'color', 'b'); hold off;
xlabel('rSAF-CW [\circ]'); ylabel('rSAF-PW [\circ]');

legend('\sigma = 0.1\circ', '\sigma = 0.2\circ', '\sigma = 0.5\circ', ...
        '\sigma = 1\circ', '\sigma = 2\circ', '\sigma = 5\circ');

set(gca, 'FontName', 'Times New Roman', 'FontSize', 20, 'FontWeight', 'bold');
set(gcf, 'Position', [507 251 354 354]);