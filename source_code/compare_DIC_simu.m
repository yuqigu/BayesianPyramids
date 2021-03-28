% simulated data
K_vec = 2:6;
DIC_for_K = [42.4733 40.8229  39.8109  40.3975  40.9127];

% promoters data
K_vec = 2:6;
DIC_for_K = [153.3115,  150.8560,  149.7564, 149.9229, 150.0055];


% K_vec = 3:5;
% DIC_for_K = [161.6851, 161.2774, 161.3209];
K_vec = 2:6;
DIC_splice_K2 = 162.0127;
DIC_splice_K3 = 161.6851;
DIC_splice_K4 = 161.2774;
DIC_splice_K5 = 161.3209;
DIC_splice_K6 = 161.3598;
DIC_for_K = [DIC_splice_K2, DIC_splice_K3, DIC_splice_K4, DIC_splice_K5, DIC_splice_K6];



figure
plot(K_vec, DIC_for_K, '-o', 'LineWidth',2, 'MarkerSize',10)
xticks(K_vec)
set(gca,'FontSize',15)
title('Complete DIC, splice data')
xlabel('# binary latent attributes')
ylabel('DIC value')
grid on
% line([4,4], [39.5,42.5], 'Color',[1 0 0])
% legend({'DIC from MCMC', 'True K'})
% print('-r300', 'DIC_splice_K', '-dpng')

