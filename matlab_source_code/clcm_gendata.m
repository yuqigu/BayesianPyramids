% generate data from constrained latent class models for multivariate categorical data
clear; close all;

rng(0513)
% -- global parameters -- %

n = 5000;  

d = 3; p = 20; G = 10; gsize = p/G;

K = 4;

S_mat_true = repmat([0 0 1 1; 0 1 0 1], p/2, 1); % p*k

beta0_true = repmat([-3 -1 0], [p 1]); % p * d
lambda0_nume = exp(beta0_true);
lambda0_true = lambda0_nume ./ sum(lambda0_nume,2);

beta_true = zeros(p,d,K);
for c = 1:d-1
beta_true(:,c,:) = repmat([1 2 4 7], [p 1]); % p * d * K
end
beta_tilde_true = beta_true .* reshape(S_mat_true,[p 1 K]); % p * K

% p * d * K
beta_array = reshape(beta0_true, [p d 1]) + beta_tilde_true;

% lambda_tilde_nume = ones(p,d,K);
lambda_tilde_nume = exp(beta_array);

lambda_tilde_true = lambda_tilde_nume ./ sum(lambda_tilde_nume,2);


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

% alpha_vec = 0.3 * ones(K,1);
pi_vec = [0.2; 0.3; 0.3; 0.2];

% consider binary responses
rep = 50;
% Y_mat = zeros(n, p, rep);
Y_arr = zeros(n, p, d, rep);



for g = 1:rep
    % latent Z_mat_true
    Z_mat_true = mnrnd(1, pi_vec, n); % n * K
    Z_vec_true = Z_mat_true * (1:K)'; % n * 1
    
    lambda_g = squeeze(lambda_tilde_true(:,:,Z_vec_true)); % p * d * n
    
    for j = 1:p
        Y_arr(:,j,:,g) = reshape(mnrnd(1,squeeze(lambda_g(j,:,:))'), [n 1 d]);
    end
    
end

dataname = strcat('clcm_new_rep50_n', num2str(n), '_p', num2str(p), '_K', num2str(K), '.mat');

save(dataname)


