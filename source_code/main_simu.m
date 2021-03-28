addpath('./from_cluster')

% load msimu_CSPbeta_n1000p20_Kini7_d4_allpos_rep50;

% look at Q in one trial
g = 4;
figure
subplot(121); imagesc(Q_true)
subplot(122); imagesc(Q_mat_pomean_rep(:,:,g));

figure
subplot(141); imagesc(Q_mat_pomean_rep(:,:,g)); colorbar
title('graphical matrix')
subplot(142); imagesc(beta_mat_pomean_rep(:,:,1,g)); colorbar
title('beta, category 1')
subplot(143); imagesc(beta_mat_pomean_rep(:,:,2,g)); colorbar
title('beta, category 2')
subplot(144); imagesc(beta_mat_pomean_rep(:,:,3,g)); colorbar
title('beta, category 3')

sig2_beta_rep(:,:,g)

beta_mat_pomean_rep

% %% CSP prior part begins
% %% traceplot
% sig2_beta_rep
% % beta variance
% figure
% plot(squeeze(sig2_beta_arr(1,1,:)))
% hold on
% plot(squeeze(sig2_beta_arr(1,2,:)))
% hold on
% plot(squeeze(sig2_beta_arr(1,3,:)))
% title('beta variance traceplot')
% 
% 
% % intercept beta0
% beta0_pomean_rep
% figure
% plot(squeeze(beta0_arr(1,1,:)))
% hold on
% plot(squeeze(beta0_arr(1,2,:)))
% hold on
% plot(squeeze(beta0_arr(1,3,:)))
% hold on
% 
% % beta_mat
% figure
% plot(squeeze(beta_mat_arr(1,1,1,:)))
% hold on
% plot(squeeze(beta_mat_arr(1,2,1,:)))
% hold on
% plot(squeeze(beta_mat_arr(1,3,1,:)))
% hold on
% 
% % K_star
% figure; plot(K_star_arr)
% 
% % tabulate(K_star_arr)
% % after burn-in
% tabulate(K_star_arr(burn+1:end))
% 
% % z_beta_vec
% % active_traits = cell(nrun, 1);
% active_traits = zeros(K, nrun);
% for ii=1:nrun
%     % active_traits{ii} = find(z_beta_vec_arr(:,ii) > (1:K)');
%     active_traits(:, ii) = (z_beta_vec_arr(:,ii) > (1:K)');
% end
% % -- CSP part ends here -- %



%% replications
% look at Q across replications
figure
for g=1:25
    subplot(5,5,g)
    imagesc(Q_mat_pomean_rep(:,:,g))
    colorbar
    title(strcat('replicate ', num2str(g)))
end

% look at beta across replications
beta_mat_ave = mean(beta_mat_pomean_rep(:,:,1:d-1,:), 4);
beta_mattil_ave = mean(beta_mattil_pomean_rep(:,:,1:d-1,:), 4);


% %
% figure
% for c=1:d-1
% subplot(3,2,2*c-1)
% imagesc(beta_mat_true(:,:,c)); colorbar; 
% title(strcat('true beta, category ', num2str(c)))
% subplot(3,2,2*c)
% imagesc(beta_mattil_ave(:,:,c)); colorbar; 
% title(strcat('estimated beta, category ', num2str(c)))
% % beta_mat_true(:,:,c) - beta_mattil_ave(:,:,c)
% end


% beta_mat_true - beta_mat_pomean_rep(:,[5 1 2 3],1:d-1,g)

rmse_beta = sqrt((beta_mat_true - beta_mattil_ave).^2);
figure
imagesc(rmse_beta(:,:,1))
hold on
imagesc(rmse_beta(:,:,2))



%% look at DIC
DIC_ave = mean(DIC_rep);




% load('msimu_DIC_Kini3_d4_n2000p20_rep25.mat')
% DIC_rep3 = DIC_rep;
% load('msimu_DIC_Kini4_d4_n2000p20_rep25.mat')
% DIC_rep4 = DIC_rep;
% load('msimu_DIC_Kini5_d4_n2000p20_rep25.mat')
% DIC_rep5 = DIC_rep;


% csvwrite('DICnew_vs_K_simu_n2000.csv', dic_matnew)


% %% DIC
% dic_mat = csvread('DIC_vs_K_simu_n2000.csv');
% K_vec = [3 4 5 6];
% 
% figure
% for g=1:size(dic_mat,1)
%     plot(K_vec, dic_mat(g,:), 'LineWidth', 1.5); hold on
% end
% 
% 
% overest = find(dic_mat(:,2) > dic_mat(:,3));
% 
% 
% figure
% for aa=1:length(overest)
%     plot(K_vec, dic_mat(overest(aa),:), 'LineWidth', 2); hold on
% end
% xticks(K_vec)


% % look at Q across replications
% figure
% for aa=1:length(overest)
%     subplot(1,length(overest),aa)
%     imagesc(Q_mat_pomean_rep(:,:,overest(aa)))
% end


% for aa=1:length(overest)
% g = overest(aa);
% figure
% subplot(141); imagesc(Q_mat_pomean_rep(:,:,g));
% subplot(142); imagesc(beta_mat_pomean_rep(:,:,1,g))
% subplot(143); imagesc(beta_mat_pomean_rep(:,:,2,g))
% subplot(144); imagesc(beta_mat_pomean_rep(:,:,3,g))
% end
% 
% 
% exact_est = setdiff(1:rep, overest);
% for aa=1:(rep-length(overest))
% g = exact_est(aa);
% figure
% subplot(141); imagesc(Q_mat_pomean_rep(:,:,g));
% subplot(142); imagesc(beta_mat_pomean_rep(:,:,1,g))
% subplot(143); imagesc(beta_mat_pomean_rep(:,:,2,g))
% subplot(144); imagesc(beta_mat_pomean_rep(:,:,3,g))
% end


% look at posterior variance of beta
for g=1:rep
    figure
    subplot(241); imagesc(Q_mat_pomean_rep(:,:,g));
    subplot(242); imagesc(beta_mat_pomean_rep(:,:,1,g))
    subplot(243); imagesc(beta_mat_pomean_rep(:,:,2,g))
    subplot(244); imagesc(beta_mat_pomean_rep(:,:,3,g))
    subplot(2,4,[6 7])
    imagesc(sig2_beta_rep(:,:,g)); colorbar
end


% %% cross validation
% figure
% for g=1:rep
%     plot(cvmat(g,:)); hold on
% end
% 
% [max_cv, which_max] = max(cvmat, [], 2);
% tabulate(which_max)




