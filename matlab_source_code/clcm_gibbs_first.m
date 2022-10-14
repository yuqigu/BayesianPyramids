function [lambda_arr, lambda0_arr, lambda_tilde_arr, pi_vec_arr, S_mat_arr] = ...
    clcm_gibbs_first(Y_arr_now, S_mat_ini, beta0, beta_mat, nrun)

% Use Stochastic Search Variable Selection to estimate S_mat

[n, p, d] = size(Y_arr_now);
K = size(beta_mat,2);


% prior for variance in SSVS
sigma2_spike = 0.01; 
sigma2_slab = 10;

% initialize rho
rho = betarnd(1,1, 1);


% initialize S_mat
S_mat = S_mat_ini;
% S_mat_ini = S_mat;


% initialize lambda0 and lambda_mat
alpha_lambda0 = 100;
alpha_lambda = 0.01;
% size p*d
lambda0_temp = gamrnd(alpha_lambda0, 1, [p d]); 
lambda0 = lambda0_temp ./ sum(lambda0_temp, 2);
% size p*d*K
lambda_temp = gamrnd(alpha_lambda, 1, [p d K]);
lambda = bsxfun(@rdivide, lambda_temp, sum(lambda_temp, 2));
lambda_tilde = bsxfun(@power, lambda, reshape(S_mat,[p 1 K])) ...
    .* bsxfun(@power, reshape(lambda0,[p d 1]), reshape(1-S_mat,[p 1 K]));


% initialize Z_mat and pi_vec
pi_temp = gamrnd(1,1, [K 1]); 
pi_vec = pi_temp / sum(pi_temp);
Z_mat = mnrnd(1, pi_vec, n); % n * K


% store MCMC output
lambda_arr = zeros(p,K,d,nrun);
lambda0_arr = zeros(p,d,nrun);
lambda_tilde_arr = zeros(p,K,d,nrun);
pi_vec_arr = zeros(K,nrun);
S_mat_arr = zeros(p,K,nrun);


for ii = 1:nrun
    
    % linear form inside the exp()
    % get Polya-Gamma parameters
    % [PG_param, ~, C_ij_minus_c] = get_linear_form(beta_mat, beta0, Q_mat, A_mat);
    [PG_param, ~, C_ij_minus_c] = get_linear_form_clcm(Z_mat, beta0, beta_mat, S_mat);
    
    keyboard
    
    % Gibbs sampling for PG varriables, beta_mat, and beta0    
    [beta_mat, beta0] = get_clcm_gibbs_beta_allpos(beta_mat, beta0, Y_arr, ...
        PG_param, C_ij_minus_c, Z_mat, S_mat);
    
    %%%%
    
    % K * p * d
    ZY_prod_Kpd = squeeze(sum(bsxfun(@times, reshape(Z_mat', [K n 1 1]), reshape(Y_arr_now, [1 n p d])), 2));
    
    % Sample S_mat  % p * d * K
    lambda_ratio_01 = reshape(lambda0, [p d 1]) ./ reshape(lambda_tilde, [p d K]);
    % % p * d-1 * K
    % lambda_ratio_01 = reshape(lambda0(:,1:d-1,:), [p d-1 1]) ./ reshape(lambda(:,1:d-1,:), [p d-1 K]);
    % p * K
    ratio_denom = squeeze(prod(bsxfun(@power, lambda_ratio_01, permute(ZY_prod_Kpd, [2 3 1])), 2));
    S_mat_prob = rho./(rho + (1-rho)*ratio_denom);
    S_mat = (rand(p,K) < S_mat_prob);
    
%     % try:
%     S_mat(:,1:2) = S_mat_true(:,1:2);
   
    % Sample lambda, size p * K * d
    lambda_temp = gamrnd(alpha_lambda0 + permute(ZY_prod_Kpd, [2 1 3]) .* reshape(S_mat, [p K 1]), 1,  [p K d]);
    lambda = lambda_temp ./ reshape(sum(lambda_temp, 3), [p K 1]);
    
    
    % Sample lambda0, size p * d
    lambda0_temp = gamrnd(alpha_lambda + squeeze(sum(permute(ZY_prod_Kpd,[2 1 3]) .* reshape(1-S_mat,[p K 1]), 2)), 1,  [p d]);
    lambda0 = lambda0_temp ./ sum(lambda0_temp, 2);    
    
    % Sample Z_mat
    % size p * K * d
    lambda_tilde = lambda .* reshape(S_mat,[p K 1]) + reshape(lambda0,[p 1 d]) .* reshape(1-S_mat,[p K 1]);
    % size n * K
    lambda_Y_nK = squeeze(prod(bsxfun(@power, reshape(lambda_tilde,[1 p K d]), reshape(Y_arr_now,[n p 1 d])), [2 4]));
    norm_const = max(log(lambda_Y_nK), [], 'all');
    lambda_Y_nK_ratio = exp(log(lambda_Y_nK) - norm_const);
    prob_Z_nume = lambda_Y_nK_ratio .* pi_vec';
    prob_Z = prob_Z_nume ./ sum(prob_Z_nume, 2);
    Z_mat = mnrnd(1, prob_Z);
    
    
    % Sample pi_vec
    pi_vec_temp = gamrnd(1 + sum(Z_mat,1),1, [1 K]);
    pi_vec = pi_vec_temp' / sum(pi_vec_temp);
    
     
    % Store data
    pi_vec_arr(:,ii) = pi_vec;
    % S_mat_arr(:,:,ii) = S_mat;
    lambda_arr(:,:,:,ii) = lambda;
    lambda0_arr(:,:,ii) = lambda0;
    lambda_tilde_arr(:,:,:,ii) = lambda_tilde;
    
    S_mat_arr(:,:,ii) = S_mat;
    
    fprintf('Gibbs iteration %d completed\n', ii);
end





end