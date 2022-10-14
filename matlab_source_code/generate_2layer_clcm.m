clear; clc; close all
rng(0513)
% -- global parameters -- %
% This function generate data from Constrained Latent Class Model (CLCM)


n = 500;  
rep = 50;
K_trait = 2; p = 20; 
d = 3;
Q_true = repmat(eye(K_trait), [p/K_trait 1]);
p = size(Q_true, 1);


% % new beta_max
beta0_true1 = repmat([-4 -3], [p 1]);

beta_max = [repmat([2 3], [12,1]); ...
            repmat([4 5], [6,1]); ...
            repmat([5 6], [2,1])];

beta_mat_true_dense = zeros(p, K_trait, d-1);
for j=1:p
    for c=1:d-1
        beta_mat_true_dense(j,:,c) = (beta_max(j,c)-beta0_true1(j,c))/sum(Q_true(j,:));
    end
end
beta_mat_true = beta_mat_true_dense .* Q_true;

% Lambda, lower layer conditional prob. tables
% A_all = binary(0:(2^K_trait-1), K_trait);
A_all = get_I(K_trait, K_trait);
Lambda_arr_true = zeros(2^K_trait, p, d);
% for CLCM
beta_true1 = zeros(p,2^K_trait,d-1);
%
for j=1:p
    beta_true1(j,:,:) = beta0_true1(j,:) + A_all * squeeze(beta_mat_true(j,:,:));
    nume_temp = exp(beta_true1(j,:,:));
    denom_temp = 1 + sum(nume_temp, 2);
    Lambda_arr_true(:,j,1:end-1) = nume_temp ./ denom_temp;
    
end
Lambda_arr_true(:,:,end) = 1 - sum(Lambda_arr_true(:,:,1:end-1),3);


%%%% transform to CLCM form %%%%
lambda_tilde_true = permute(Lambda_arr_true, [2 3 1]);
S_mat_true = get_ideal_resp(Q_true, A_all);
beta_true = zeros(p,d,2^K_trait);
beta_true(:,1:d-1,:) = permute(beta_true1, [1 3 2]);
beta0_true = zeros(p,d);
beta0_true(:,1:d-1) = beta0_true1;


figure;
for j=1:10
subplot(10,2,(j-1)*2+1)
imagesc(S_mat_true(j,:)); colorbar; %title(num2str(j))
% title(strcat('true E-vector for item', num2str(j)))
% xlabel('extreme latent profiles');
% xticks(1:k); 
%
subplot(10,2,(j-1)*2+2)
imagesc(squeeze(lambda_tilde_true(j,:,:))); colorbar;
%title(num2str(j))
% title(strcat('true lambda matrix for item', num2str(j)))
% xlabel('extreme latent profiles'); ylabel('response categories')
% xticks(1:k); yticks(1:d); caxis([0 0.5])
end



% alpha_vec = 0.3 * ones(K_trait,1);
pi_vec = [0.2; 0.3; 0.3; 0.2];


Y_arr = zeros(n, p, d, rep);


K_all = 2^K_trait;
for g = 1:rep
    % latent Z_mat_true
    Z_mat_true = mnrnd(1, pi_vec, n); % n * K_trait
    Z_vec_true = Z_mat_true * (1:K_all)'; % n * 1
    
    lambda_g = squeeze(lambda_tilde_true(:,:,Z_vec_true)); % p * d * n
    
    for j = 1:p
        Y_arr(:,j,:,g) = reshape(mnrnd(1,squeeze(lambda_g(j,:,:))'), [n 1 d]);
    end
    
end

dataname = strcat('clcm_fromBP_rep50_n', num2str(n), '_p', num2str(p), '_K', num2str(K_trait), '.mat');

save(dataname)
