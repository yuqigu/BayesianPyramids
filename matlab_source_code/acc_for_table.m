raw = [91, 146, 160, 84, 134, 161]/1000;
dx = [54, 45, 72, 50, 35, 71]/1000;
indbin = [36, 45, 44, 35, 43, 35]/1000;
bp = [29, 25, 30, 16, 17, 24]/1000;

r_raw = 1- bp ./ raw
r_dx = 1- bp ./ dx
r_indbin = 1- bp ./ indbin

[min(r_raw(4:6)), max(r_raw(4:6))]
[min(r_dx(4:6)), max(r_dx(4:6))]
[min(r_indbin(4:6)), max(r_indbin(4:6))]

% bar plot 
err_mat = [raw', dx', indbin', bp'];

figure; bar(1-err_mat(4:6,:)); ylim([0.8,1])

figure
for j = 1:4
    plot(1-err_mat(4:6,j), ':p', 'MarkerSize', 16, 'LineWidth', 2)
    hold on
end
xticks([1 2 3 4])
xticklabels({'EI', 'IE', 'Nither'})
set(gca, 'FontSize', 20)
legend({'Raw DNA sequences', 'Univariate latent class', 'Single latent layer', ...
    'Bayesian Pyramid'})
xlabel('Three sequence types')
ylabel('classification accuracy')
title('Downstream classification accuracy on the test set')
print('-r300', 'splice_compare', '-dpng')
