addpath('./from_cluster')

% This script plots figures that visualize the simulation results
load msimu_weak_CSPbeta_n2000p20_Kini7_d4_allpos_rep50;


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

% beta_mat_true - beta_mat_pomean_rep(:,[5 1 2 3],1:d-1,g)

rmse_beta = sqrt((beta_mat_true - beta_mattil_ave(:,1:K0,:)).^2);
figure
imagesc(rmse_beta(:,:,1))
hold on
imagesc(rmse_beta(:,:,2))








