function [lambda0_arr, lambda_tilde_arr, pi_vec_arr, S_mat_arr, beta_tilde_arr, beta0_arr, z_vec_arr] = ...
    clcm_gibbs_logistic(Y_arr_now, S_mat_ini, beta0, beta, nrun)

% Use Stochastic Search Variable Selection to estimate S_mat

[n, p, d] = size(Y_arr_now);
K = size(S_mat_ini,2);


% prior for variance in SSVS
sig2_spike = 0.1; 
sig2_slab = 10;

% initialize rho
rho = betarnd(1,1, 1);


% initialize S_mat
S_mat = S_mat_ini;
% S_mat_ini = S_mat;


% % initialize lambda0 and lambda_mat
% alpha_lambda0 = 100;
% alpha_lambda = 0.01;
% % size p*d
% lambda0_temp = gamrnd(alpha_lambda0, 1, [p d]); 
% lambda0 = lambda0_temp ./ sum(lambda0_temp, 2);
% % size p*d*K
% lambda_temp = gamrnd(alpha_lambda, 1, [p d K]);
% lambda = bsxfun(@rdivide, lambda_temp, sum(lambda_temp, 2));
% lambda_tilde = bsxfun(@power, lambda, reshape(S_mat,[p 1 K])) ...
%     .* bsxfun(@power, reshape(lambda0,[p d 1]), reshape(1-S_mat,[p 1 K]));


% initialize Z_mat and pi_vec
pi_temp = gamrnd(1,1, [K 1]); 
pi_vec = pi_temp / sum(pi_temp);
Z_mat = mnrnd(1, pi_vec, n); % n * K


% store MCMC output
lambda0_arr = zeros(p,d,nrun);
lambda_tilde_arr = zeros(p,d,K,nrun);
pi_vec_arr = zeros(K,nrun);
S_mat_arr = zeros(p,K,nrun);
z_vec_arr = zeros(n,nrun);

beta_tilde_arr = zeros(p,d,K,nrun);
beta0_arr = zeros(p,d,nrun);

for ii = 1:nrun
    
    % linear form inside the exp()
    % get Polya-Gamma parameters
    [PG_param, ~, C_ij_minus_c] = get_linear_form_clcm(Z_mat, beta0, beta, S_mat);
        
    % Gibbs sampling for PG varriables, beta, and beta0    
    [beta, beta0] = get_clcm_gibbs_logit_beta(beta, beta0, Y_arr_now, ...
        PG_param, C_ij_minus_c, Z_mat, S_mat);
    
    beta_tilde = beta .* reshape(S_mat,[p 1 K]);
    beta_array = reshape(beta0, [p d 1]) + beta_tilde;
    lambda_tilde_nume = exp(beta_array);
    % size p * d * K
    lambda_tilde = lambda_tilde_nume ./ sum(lambda_tilde_nume,2);
    % size p * d
    lambda0_nume = exp(beta0);
    lambda0 = lambda0_nume ./ sum(lambda0_nume,2);
    %%%%
   

    % Sample S_mat 
    for k = 1:K
        % size p*1
        S_mat_k1 = S_mat; S_mat_k1(:, k) = 1;
        S_mat_k0 = S_mat; S_mat_k0(:, k) = 0;

        % get linear forms, n * p * d 
        [~, linear_form_k1, ~] = get_linear_form_clcm(Z_mat, beta0, beta, S_mat_k1);
        [~, linear_form_k0, ~] = get_linear_form_clcm(Z_mat, beta0, beta, S_mat_k0);
    
        % normalize_lf_k1 = linear_form_k1 - log(sum(exp(linear_form_k1),3));
        normalize_lf_k1 = linear_form_k1 - max(linear_form_k1, [], 3);
        normalize_lf_k0 = linear_form_k0 - max(linear_form_k0, [], 3);

        % min(normalize_lf_k1,[],3)

        % size n_now * p
        denom_category_h1 = sum(exp(normalize_lf_k1), 3);
        denom_category_h0 = sum(exp(normalize_lf_k0), 3);

        % size n_now * p * d; sum(cond_prob_h1, 3) = 1
        cond_prob_h1 = bsxfun(@rdivide, exp(normalize_lf_k1), denom_category_h1);
        cond_prob_h0 = bsxfun(@rdivide, exp(normalize_lf_k0), denom_category_h0);

        % size n_now * p
        prod_part_h1 = prod( bsxfun(@power, cond_prob_h1, Y_arr_now), 3);
        prod_part_h0 = prod( bsxfun(@power, cond_prob_h0, Y_arr_now), 3);
        
%         prior_ratio_01 = sig2_spike^(-(d-1)/2) / sig2_slab^(-(d-1)/2) ...
%             * prod(exp(-1/2 * squeeze(beta(:,:,k).^2) * (sig2_spike^(-1) - sig2_slab^(-1))'), 2);
% 
%         likelihood_ratio_01 = exp(sum(log(prod_part_h0) - log(prod_part_h1), 1)');
%         
%         ratio_qk_01 = prior_ratio_01 .* likelihood_ratio_01;
        ratio_qk_01 = sig2_spike^(-(d-1)/2) / sig2_slab^(-(d-1)/2) ...
            * exp( sum(-1/2 * squeeze(beta(:,1:d-1,k).^2) * (sig2_spike^(-1) - sig2_slab^(-1))',2) ...
            + sum(log(prod_part_h0) - log(prod_part_h1), 1)' );

        prob_Sk = 1./(1 + (1-rho)/rho * ratio_qk_01);
        S_mat(:, k) = (rand(p,1) < prob_Sk);
    end
    
    
    % Sample Z_mat

    % size n * K
    lambda_Y_nK = squeeze(prod(bsxfun(@power, reshape(lambda_tilde,[1 p d K]), reshape(Y_arr_now,[n p d 1])), [2 3]));
    norm_const = max(log(lambda_Y_nK), [], 'all');
    lambda_Y_nK_ratio = exp(log(lambda_Y_nK) - norm_const);
    prob_Z_nume = lambda_Y_nK_ratio .* pi_vec';
    prob_Z = prob_Z_nume ./ sum(prob_Z_nume, 2);
    Z_mat = mnrnd(1, prob_Z);
    z_vec_arr(:,ii) = Z_mat * (1:K)';
    
    
    % Sample piSvec
    pi_vec_temp = gamrnd(1 + sum(Z_mat,1),1, [1 K]);
    pi_vec = pi_vec_temp' / sum(pi_vec_temp);
    
     
    % Store data
    pi_vec_arr(:,ii) = pi_vec;
    lambda0_arr(:,:,ii) = lambda0;
    lambda_tilde_arr(:,:,:,ii) = lambda_tilde;
    
    S_mat_arr(:,:,ii) = S_mat;
    
    beta_tilde_arr(:,:,:,ii) = beta_tilde;
    beta0_arr(:,:,ii) = beta0;
    
    %if mod(ii,100)==0
        fprintf('Gibbs iteration %d completed\n', ii);
    %end
end





end