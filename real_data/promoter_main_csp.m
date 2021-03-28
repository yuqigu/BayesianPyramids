% -- posterior means -- %
beta_mat_pomean = mean(beta_mat_arr(:,:,:,burn+1:thin:end), 4);
beta_mat_pomean


sig2_beta_pomean = mean(sig2_beta_arr(:,:,burn+1:thin:end), 3);
sig2_beta_pomean


beta0_pomean = mean(beta0_arr(:,:,burn+1:thin:end), 3);
beta0_pomean


% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);
Bern_K_pomean


% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:thin:end), 2);
fprintf('\ntau posterior mean is:\n')
tau_pomean


figure; 
plot(tau_arr(1,:))
hold on
plot(tau_arr(2,:))
title('tau traceplot')

figure; plot(squeeze(Bern_K_arr(2,2,:)))


% A_mat post. mean
A_post_mean_prob = mean(A_mat_arr(:,1:K,burn+1:nrun), 3);
A_mat_pomean = (A_post_mean_prob > 0.5);

% Q_mat post. mean
Q_post_mean_prob = mean(Q_mat_arr(:,1:K,burn+1:nrun), 3);
Q_mat_pomean = (Q_post_mean_prob > 0.5);


fprintf('\nlocation loading structure:\n')
[(1:p)', Q_mat_pomean]


% z_mat post. mean
z_tau_mat_pomean = mean(z_tau_mat_arr(:,:,burn+1:end), 3);
z_tau_pm_binary = (z_tau_mat_pomean>0.5);

acc_col1 = mean(is_promoter == z_tau_pm_binary(:,1));
acc_col2 = mean(is_promoter == z_tau_pm_binary(:,2));

acc_gene = max(acc_col1, acc_col2);

if acc_col1 > acc_col2
    choose_col = 1;
else
    choose_col = 2;
end

[(1:n)', is_promoter, z_tau_pm_binary(:,choose_col)]


acc_promo = mean(is_promoter(1:(n/2)) == z_tau_pm_binary(1:(n/2),choose_col));
acc_nonpromo = mean(is_promoter((n/2+1):end) == z_tau_pm_binary((n/2+1):end,choose_col));


% fprintf('\nOverall clustering accuracy is %1.4f\n', acc_gene);
% fprintf('\nPromoters clustering accuracy is %1.4f\n', acc_promo);
% fprintf('\nNon-promoters clustering accuracy is %1.4f\n', acc_nonpromo);



% sum(z_tau_pm_binary, 1)



% [is_promoter, (1:n)', A_mat_pomean]


figure
subplot(121)
imagesc(Q_mat_pomean)
xticks(1:K); yticks(1:p)
xlabel('binary latent vars')
ylabel('p=57 locations')
set(gca, 'FontSize', 12)
title('1st layer location loading structure')
%
subplot(122)
imagesc([A_mat_pomean])
xticks(1:K); yticks(1:3:n); 
xlabel('binary latent vars')
ylabel('n=106 genes')
set(gca, 'FontSize', 12)
title('1st layer gene grouping structure')
% print('-r300', 'promo_tauDir18_7925', '-dpng')



% % -- traceplots -- %
% figure
% plot(squeeze(Bern_K_arr(1,1,:)))
% hold on
% plot(squeeze(Bern_K_arr(1,2,:)))
% hold on
% plot(squeeze(Bern_K_arr(2,1,:)))
% hold on
% plot(squeeze(Bern_K_arr(2,2,:)))
% hold on
% plot(squeeze(Bern_K_arr(3,1,:)))
% hold on
% plot(squeeze(Bern_K_arr(3,2,:)))
% hold on
% plot(squeeze(Bern_K_arr(4,1,:)))
% hold on
% plot(squeeze(Bern_K_arr(4,2,:)))
% xlabel('Gibbs iterations');
% title('deep tensor arms'); set(gca, 'FontSize', 12)
% % print('-r300', 'trace_deep_arm', '-dpng')




% ------- Rule List preparation begins here ------- %
% conver z_tau and A_mat_pomean to multivariate binary rule lists

ZA = [z_tau_pm_binary, A_mat_pomean];

fprintf('\nOverall clustering accuracy is %1.4f\n', acc_gene);
fprintf('\nPromoters clustering accuracy is %1.4f\n', acc_promo);
fprintf('\nNon-promoters clustering accuracy is %1.4f\n', acc_nonpromo);

% -- write to csv file -- %
% csvwrite('Y_promoters_K6.csv', [ZA, is_promoter])

DIC

figure; 
subplot(121); imagesc(is_promoter)
subplot(122); imagesc(ZA)


%%
traits_select = find(z_beta_vec > (1:K)');

figure
subplot(121)
imagesc(Q_mat_pomean(:,traits_select))
xticks(1:K); yticks(1:p)
xlabel('binary latent vars')
ylabel('p=57 locations')
set(gca, 'FontSize', 12)
title('1st layer location loading structure')
%
subplot(122)
imagesc([A_mat_pomean(:,traits_select)])
xticks(1:K); yticks(1:3:n); 
xlabel('binary latent vars')
ylabel('n=106 genes')
set(gca, 'FontSize', 12)
title('1st layer gene grouping structure')
% print('-r300', 'promo_tauDir18_7925', '-dpng')

figure
imagesc(z_tau_mat_pomean); colorbar

figure
imagesc(Q_post_mean_prob(:,traits_select)); colorbar


%% pink figure
pink3 = pink(3);
reversepink = [pink3(3,:); pink3(1,:)];

figure
ax1 = subplot(131);
imagesc(Q_mat_pomean(:,traits_select));
xticks(1:length(traits_select))
xlabel('binary latent vars'); ylabel('p=57 locations'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old1 = colormap(ax1,reversepink); 
bar1 = colorbar; bar1.Ticks = [0.25 0.75]; bar1.TickLabels = {'0' '1'};
title('loci loading')
% second plot
ax2 = subplot(132);
imagesc(A_mat_pomean(:,traits_select));
xticks(1:length(traits_select))
xlabel('binary latent vars'); ylabel('n=106 sequences'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old2 = colormap(ax2,reversepink); colormap( flipud(old2) ); 
bar2 = colorbar; bar2.Ticks = [0.25 0.75]; bar2.TickLabels = {'0' '1'};
line([0 5], [53 53], 'Color', [0.5 0.5 0.5]); yticks([10:10:40, 53, 60:10:n])
title('1st layer traits')
% third plot
ax3 = subplot(133);
imagesc(z_tau_mat_pomean);
xlabel('z membership'); ylabel('n=106 sequences'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12); 
old3 = colormap(ax3,reversepink); colormap( flipud(old3) ); 
bar3 = colorbar; bar3.Ticks = [0.25 0.75]; bar3.TickLabels = {'0' '1'};
line([0 5], [53 53], 'Color', [0.5 0.5 0.5]); 
yticks([10:10:40, 53, 60:10:n])
title('2nd layer membership')

% print('-r300', 'promoter_csp_Kupper_7_Kfinal_4', '-dpng')


ZA_eff = [z_tau_mat_pomean, A_mat_pomean(:,traits_select)];
csvwrite('Y_promoter_K_upper7_alpha2_eff.csv', [ZA_eff, is_promoter])


A_mat_pomean(10:17,traits_select)

tabulate(sum(A_mat_pomean(:,traits_select), 2))

A_post_mean_prob(10:17,traits_select)

figure; imagesc(A_post_mean_prob(:, traits_select))


% -- plot confusion matrix -- %
t_hat = A_mat_pomean(:,4);
conf_mat = zeros(B, B);

for bb=1:B
    for cc=1:B
        conf_mat(bb,cc) = sum(is_promoter==bb-1 & t_hat==cc-1) / sum(is_promoter==bb-1);
    end
end


conf_mat





