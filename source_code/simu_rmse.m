n_vec = 250:250:2000;

beta_rmse = [1.2201    1.4314    1.5641;...
0.3556    0.3867    0.5986;...
0.2565    0.2723    0.2988;...
0.2137    0.2304    0.2450;...
0.1818    0.1933    0.2091;...
0.1587    0.1707    0.1905;...
0.1527    0.1620    0.1738;...
0.1424    0.1500    0.1626];

beta0_rmse = [1.1931    1.3210    1.4782;...
0.5914    0.6335    0.7024;...
0.4258    0.4515    0.4970;...
0.3450    0.3761    0.4220;...
0.2949    0.3335    0.3558;...
0.2707    0.3078    0.3795;...
0.2470    0.2821    0.3306;...
0.2236    0.2506    0.3172];


eta_rmse = [0.3500    0.3610    0.3991;...
0.0729    0.0838    0.2251;...
0.0427    0.0486    0.0535;...
0.0339    0.0380    0.0439;...
0.0263    0.0316    0.0375;...
0.0258    0.0303    0.0366;...
0.0235    0.0264    0.0296;...
0.0211    0.0246    0.0269];


% error bar
% 
beta_rmse_errbar = beta_rmse;
beta_rmse_errbar(:,1) = beta_rmse(:,2) - beta_rmse(:,1);
beta_rmse_errbar(:,3) = beta_rmse(:,3) - beta_rmse(:,2);


beta0_rmse_errbar = beta_rmse;
beta0_rmse_errbar(:,1) = beta0_rmse(:,2) - beta0_rmse(:,1);
beta0_rmse_errbar(:,3) = beta0_rmse(:,3) - beta0_rmse(:,2);


eta_rmse_errbar = eta_rmse;
eta_rmse_errbar(:,1) = eta_rmse(:,2) - eta_rmse(:,1);
eta_rmse_errbar(:,3) = eta_rmse(:,3) - eta_rmse(:,2);


%% plot aRMSE with error bars
figure
errorbar(n_vec, beta_rmse_errbar(:,2), beta_rmse_errbar(:,1), beta_rmse_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10, 'CapSize', 12)
% xticks(n_vec)
xticks(500:500:2000)
xlim([200, 2050]); ylim([0, 2])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('RMSE for main-effects $(\beta_{j,k,c})$', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'rmse_beta_mat', '-dpng')



%
figure
errorbar(n_vec, beta0_rmse_errbar(:,2), beta0_rmse_errbar(:,1), beta0_rmse_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10, 'CapSize', 12)
% xticks(n_vec)
xticks(500:500:2000)
xlim([200, 2050]); 
ylim([0, 2])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('RMSE for intercepts $(\beta_{j,0,c})$', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'rmse_beta0', '-dpng')



%
figure
errorbar(n_vec, eta_rmse_errbar(:,2), eta_rmse_errbar(:,1), eta_rmse_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10, 'CapSize', 12)
% xticks(n_vec)
xticks(500:500:2000)
xlim([200, 2050]); 
ylim([0, 0.6])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('RMSE for deeper tensor arms $(\eta_{k,b})$', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'rmse_eta', '-dpng')


%%
% -- look at the estimation accracy of the graphical matrix G -- %

Q_acc_mat_n = [1,0.4728;...
0.56,0.0792;...
0.24,0.0335;...
0.08,0.0085;...
0.04,0.0055;...
0.04,0.002;...
0.02,0.0003;...
0.02,0.0003];

n_vec = 250:250:2000;

figure
plot(n_vec, Q_acc_mat_n(:,1), ':o', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.9290, 0.6940, 0.1250])
hold on
plot(n_vec, Q_acc_mat_n(:,2), '-s', 'MarkerSize', 10, 'LineWidth', 2, 'Color', [0.4940, 0.1840, 0.5560])
legend({'Matrix error', 'Entry error'})
xticks(n_vec)
xlim([200, 2050]);
xlabel('sample size'); ylabel('average estimation error')
pbaspect([4 3 1])
title('Estimation error Graphical Matrix')
set(gca, 'FontSize', 18)

% print('-r300', 'Q_acc', '-dpng')





% -- median; look at the estimation accracy of the graphical matrix G -- %
Q_err_matrix = [ 1     1     1;...
0     1     1;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0];

Q_err_entry = [0.4000    0.4938    0.6375;...
0    0.0125    0.0625;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0];

Q_err_row = [0.8500    1.0000    1.0000;...
0    0.0500    0.2500;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0;...
0     0     0];

Q_err_matrix_errbar = Q_err_matrix;
Q_err_matrix_errbar(:,1) = Q_err_matrix(:,2) - Q_err_matrix(:,1);
Q_err_matrix_errbar(:,3) = Q_err_matrix(:,3) - Q_err_matrix(:,2);

Q_err_entry_errbar = Q_err_entry;
Q_err_entry_errbar(:,1) = Q_err_entry(:,2) - Q_err_entry(:,1);
Q_err_entry_errbar(:,3) = Q_err_entry(:,3) - Q_err_entry(:,2);

Q_err_row_errbar = Q_err_row;
Q_err_row_errbar(:,1) = Q_err_row(:,2) - Q_err_row(:,1);
Q_err_row_errbar(:,3) = Q_err_row(:,3) - Q_err_row(:,2);



n_vec = 250:250:2000;

figure
errorbar(n_vec, Q_err_matrix_errbar(:,2), Q_err_matrix_errbar(:,1), Q_err_matrix_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10)
% plot(n_vec, Q_err_matrix_errbar(:,2), ':o', 'MarkerSize', 10, 'LineWidth', 2)
xticks(n_vec); xlim([200, 2050]); ylim([0 1])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('Estimation error of $G^{(1)}$, matrix-wise', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'Q_acc_matrix_bar', '-dpng')


figure
errorbar(n_vec, Q_err_row_errbar(:,2), Q_err_row_errbar(:,1), Q_err_row_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10, 'Color', [0.8500, 0.3250, 0.0980])
xticks(n_vec); xlim([200, 2050]); ylim([0 1])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('Estimation error of $G^{(1)}$, row-wise', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'Q_acc_row', '-dpng')


figure
errorbar(n_vec, Q_err_entry_errbar(:,2), Q_err_entry_errbar(:,1), Q_err_entry_errbar(:,3), 'o', ...
    'LineWidth', 2, 'LineStyle', '-', 'MarkerSize', 10, 'Color', [0.4940, 0.1840, 0.5560])
xticks(n_vec); xlim([200, 2050]); ylim([0 1])
xlabel('sample size'); ylabel('average RMSE')
pbaspect([4 3 1])
title('Estimation error of $G^{(1)}$, entry-wise', 'interpreter','latex')
set(gca, 'FontSize', 22)
print('-r300', 'Q_acc_entry', '-dpng')