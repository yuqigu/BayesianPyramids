% filename = strcat('mdata0_n', num2str(n), '_p', num2str(p), '_K', num2str(K), '_d', num2str(d), '.mat');
% save(filename)

beta_mat_tilde_arr = zeros(p, K, d, nrun);
for ii=1:nrun
    beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
end


% -- posterior means -- %
% beta_matcat1_pomean = mean(beta_mat_tilde_arr(:,1:K,1,burn+1:end), 4);
% beta_matcat2_pomean = mean(beta_mat_tilde_arr(:,1:K,2,burn+1:end), 4);
% beta_matcat3_pomean = mean(beta_mat_tilde_arr(:,1:K,3,burn+1:end), 4);
beta_mat_pomean = zeros(p,K,d-1);
for c=1:d-1
    beta_mat_pomean(:,:,c) = mean(beta_mat_tilde_arr(:,1:K,c,burn+1:end), 4);
end


% % beta posterior mean
% figure
% subplot(131); imagesc(beta_mat_pomean(:,:,1)); colorbar
% subplot(132); imagesc(beta_mat_pomean(:,:,2)); colorbar
% subplot(133); imagesc(beta_mat_pomean(:,:,3)); colorbar
% imagesc([beta_matcat1_pomean, beta_matcat2_pomean, beta_matcat3_pomean]); colorbar
% % print('-r300', 'beta_bindata_est', '-dpng')

beta0_pomean = mean(beta0_arr(:,:,burn+1:nrun), 3);
beta0_pomean

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(1:K,:,burn+1:nrun), 3);
Bern_K_pomean

% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:nrun), 2);
fprintf('\ntau posterior mean is:\n')
tau_pomean

figure; 
plot(tau_arr(1,:))
hold on
plot(tau_arr(2,:))
title('tau traceplot')


% Q post. mean
Q_post_mean = (mean(Q_mat_arr(:,1:K,burn+1:nrun), 3) > 0.5);


fprintf('\nlocation loading structure:\n')
[(1:p)', Q_post_mean]


% z_mat post. mean
z_mat_pomean = mean(z_mat_arr(:,:,burn+1:end), 3);
z_tau_pm_binary = (z_mat_pomean>0.5);

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

% A_mat post. mean
A_post_mean_prob = mean(A_mat_arr(:,1:K,burn+1:nrun), 3);
A_post_mean = (A_post_mean_prob > 0.5);

A_mat_arr(:,1:K,end) - A_mat_arr(:,1:K,end-1);
A_mat_arr(:,1:K,end-1) - A_mat_arr(:,1:K,end-2);


[is_promoter, (1:n)', A_post_mean]


figure
subplot(121)
imagesc(Q_post_mean)
xticks(1:K); yticks(1:p)
xlabel('binary latent vars')
ylabel('p=57 locations')
set(gca, 'FontSize', 12)
title('1st layer location loading structure')
%
subplot(122)
imagesc(A_post_mean)
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




mean(any(A_post_mean, 2) == is_promoter)
% [(1:n)', is_promoter, z_tau_pm_binary(:,choose_col), A_post_mean]


% ------- Rule List preparation begins here ------- %
% conver z_tau and A_post_mean to multivariate binary rule lists

ZA = [z_tau_pm_binary, A_post_mean];

fprintf('\nOverall clustering accuracy is %1.4f\n', acc_gene);
fprintf('\nPromoters clustering accuracy is %1.4f\n', acc_promo);
fprintf('\nNon-promoters clustering accuracy is %1.4f\n', acc_nonpromo);

% -- write to csv file -- %
% csvwrite('Y_promoters_K6.csv', [ZA, is_promoter])

DIC_rep



