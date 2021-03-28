addpath('./from_cluster')
load msimu_weak_CSPbeta_n250p20_Kini7_d4_allpos_rep50;

active_traits = zeros(K, rep);

Q_active_all = zeros(p, K, rep);
beta_mat_active = zeros(p, K, d, rep);
beta_mattil_active = zeros(p, K, d, rep);

K_active = zeros(1, rep);

Bern_K_permute = zeros(K, B, rep);

% further permute Bern_K
Bern_K_permute2 = zeros(K, B, rep);
tau_permute2 = zeros(B, rep);

for g=1:rep
    Q_mat_g = Q_mat_pomean_rep(:,:,g);
    % active_traits(:, g) = z_beta_vec_rep(:,g) > (1:K)';
    
    % % select latent traits based on posterior uncertainty
    % active_traits(:, g) = (mean(Q_mat_g .* (1-Q_mat_g), 1) < 0.01)';
    
    % -- find the active traits -- %
    [sort_mean_var, rank_mean_var] = sort(mean(sig2_beta_rep(:,:,g), 2), 'descend');
    %
    [~, K_act] = max(sort_mean_var(1:end-1) - sort_mean_var(2:end));
    K_active(g) = K_act;
    %
    active_traits(rank_mean_var(1:K_act), g) = 1;
    
    find_active = find(active_traits(:, g));
    
    Q_binary_pomean = Q_mat_g(:, active_traits(:,g) == 1) > 0.5;
        
    % find the column permutation starts
    Q_firstK_binind = 2.^(0:K_act-1) * Q_binary_pomean(1:K_act, :);
    [~, permute_vec] = sort(Q_firstK_binind);
    permute_augment = [find_active(permute_vec); setdiff((1:K)', find_active)];
    
    Q_active_all(:, :, g) = Q_mat_pomean_rep(:, permute_augment, g);
    
    beta_mat_active(:,:,:, g) = beta_mat_pomean_rep(:, permute_augment, :, g);
    beta_mattil_active(:,:,:, g) = beta_mattil_pomean_rep(:, permute_augment, :, g);
    
    Bern_K_permute(:,:,g) = Bern_K_pomean_rep(permute_augment, :, g);
    
    % second permutation of tau and Bern_K (eta)
    if Bern_K_permute(1,1,g) > Bern_K_permute(1,2,g)
        Bern_K_permute2(:,:,g) = Bern_K_permute(:,:,g);
        tau_permute2(:,g) = tau_pomean_rep(:,g);
    else
        Bern_K_permute2(:,:,g) = [Bern_K_permute(:,2,g), Bern_K_permute(:,1,g)];
        tau_permute2(:,g) = [tau_pomean_rep(2,g); tau_pomean_rep(1,g)];
    end
    
    % permucation at the deepest layer for tau
    
end



% RMSE
se_beta = zeros(p,K0,d-1,rep);
se_beta0 = zeros(p,d-1,rep);
se_Bern = zeros(K0,B,rep);
se_tau = zeros(1,rep);


for g=1:rep
se_beta(:,:,:,g) = (beta_mat_true - beta_mat_active(:,1:K0,1:d-1,g)).^2;

se_beta0(:,:,g) = (beta0_true - beta0_pomean_rep(:,1:d-1,g)).^2;

% se_tau(:,g) = (tau_true - tau_pomean_rep(:,g)).^2;
% se_tau(:,g) = (tau_true - tau_permute2(:,g)).^2;
se_tau(g) = (tau_true(1) - tau_permute2(1,g)).^2;


% se_Bern(:,:,g) = (Bern_K_true - Bern_K_permute(1:K0,:,g)).^2;
se_Bern(:,:,g) = (Bern_K_true - Bern_K_permute2(1:K0,:,g)).^2;

end



%%%
RMSE_beta_mat = sqrt(squeeze(sum(squeeze(sum(squeeze(sum(se_beta, 1)), 1)), 1))/(p*K0*(d-1)));

RMSE_beta0 = sqrt(squeeze(sum(squeeze(sum(se_beta0, 1)), 1))/(p*(d-1)));

RMSE_Bern = sqrt(squeeze(sum(squeeze(sum(se_Bern, 1)), 1))/(K0*B));

RMSE_tau = sqrt(squeeze(sum(se_tau, 1)));


% figure; boxplot(RMSE_beta_mat); title('beta')
% 
% figure; boxplot(RMSE_beta0); title('beta0')
% 
% figure; boxplot(RMSE_Bern); title('eta')
% 
% figure; boxplot(RMSE_tau); title('tau')


% % look at Bern_K
% Bern_K_permute2(1:K0,:,:)


%
RMSE_mat = [RMSE_beta_mat; RMSE_beta0; RMSE_Bern; RMSE_tau];

[median(RMSE_mat, 2), iqr(RMSE_mat, 2)]

quantile(RMSE_mat, [0.25 0.5 0.75], 2)

% mean(RMSE_mat, 2)



% -- look at the estimation accracy of the graphical matrix G -- %
K_rep = sum(active_traits, 1);
sum(sum(active_traits,1) == K0)

Q_round_rep = zeros(p, K0, rep);

Q_acc_mat = zeros(rep,1);
Q_acc_row = zeros(rep,1);
Q_acc_entry = zeros(rep,1);

for g=1:rep
    Q_round_rep(:,:,g) = (Q_active_all(:,1:K0,g) > 0.5);
    Q_acc_mat(g) = all(all(Q_round_rep(:,:,g) == Q_true));
    Q_acc_entry(g) = sum(sum(Q_round_rep(:,:,g) == Q_true))/(p*K0);
    
    Q_acc_row(g) = sum(all(Q_round_rep(:,:,g) == Q_true, 2))/p;
end

Q_error2 = 1 - [Q_acc_mat, Q_acc_entry, Q_acc_row];

n

% 1 - [mean(Q_acc_mat), mean(Q_acc_entry)]

% [quantile(Q_error2, 0.25);...
% median(Q_error2);...
% quantile(Q_error2, 0.75)]


% matrix error
[quantile(Q_error2(:,1), 0.25),...
median(Q_error2(:,1)),...
quantile(Q_error2(:,1), 0.75)]


% entry error
[quantile(Q_error2(:,2), 0.25),...
median(Q_error2(:,2)),...
quantile(Q_error2(:,2), 0.75)]


% row error
[quantile(Q_error2(:,3), 0.25),...
median(Q_error2(:,3)),...
quantile(Q_error2(:,3), 0.75)]

