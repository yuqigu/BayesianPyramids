function [Q_mat_arr, beta_mat_arr, beta0_arr, A_mat_arr, Bern_K_arr, tau_arr,...
    z_tau_mat_arr, sig2_beta_arr, K_star_arr, z_beta_vec_arr] ...
    = bp_csp_fun(Y_arr, K, B, nrun, alpha0)

% This function runs the Gibbs sampler for Bayesian Pyramid with CSP prior

[n, p, d] = size(Y_arr);


%% specify the priors
% prior for intercept
mu_beta0 = -3 * ones(3,1);
var_beta0 = 10 * ones(d-1,1);

% inverse Gamma prior for hypervariances
a_sig = 2; b_sig = 2;

% sig2_pseudo = v0^2;
sig2_pseudo = 0.07;

% -- initialize deep layer parameters -- %
tau = ones(1,B)/B;
Bern_K = rand(K, B); 
% -- initialize local deep Z -- %
z_tau_mat = mnrnd(1, tau, n);

% -- initialize local binary patterns A -- %
A_mat = (rand(K, n) < Bern_K(:, z_tau_mat*(1:B)'))';


% % Option 2: give gamma_q a prior
gamma_q = betarnd(1,1, 1);

% -- initialize Q matrix -- %
Q_mat = (rand(p,K) < gamma_q');
% make sure every row contains at least one entry of 1
% for any all-zero row, uniformly sample one entry to make it =1 
if any(sum(Q_mat,2)==0)
    Q_mat(sum(Q_mat,2)==0, :) = mnrnd(1, 1/K*ones(K,1), sum(sum(Q_mat,2)==0));
end

% -- CSP prior begins -- %
sig2_inf = 0.1;

% the last component of v_beta equals 1; v_beta has size K * 1
v_beta = [betarnd(1, alpha0, [K-1 1]); 1];
w_beta = v_beta .* cumprod([1; 1-v_beta(1:end-1)]);
% size K * K
z_beta_mat = mnrnd(1, w_beta, K);
% -- CSP prior ends -- %


% -- initialize sig2_beta -- %
% sig2_beta, size K * d-1
sig2_beta = zeros(K, d-1);
for k=1:K
    if find(z_beta_mat(k,:)) > k
        % slab
        sig2_beta(k, :) = 1./(gamrnd(a_sig, 1/b_sig, [1, d-1]));
    else
        % spike
        sig2_beta(k, :) = sig2_inf;
    end
end


% -- initialize beta-mat -- %
% the p variables for each (c,k) are a priori independent, so it has a
% diagonal covariance structure as follows
beta_mat = zeros(p, K, d);
for c=1:d-1
    for k=1:K
        pd1 = makedist('Normal', 'mu', 0, 'sigma', sqrt(sig2_beta(k,c)));
        beta_mat(:, k, c) = random(truncate(pd1,0,Inf), p, 1);
    end
end


% initialize intercept beta0, size p * d -- %
beta0 = zeros(p, d);
beta0(:, 1:d-1) = repmat(mu_beta0', [p,1]);

% -- book-keeping of the data through MCMC iterations -- %
beta_mat_arr = zeros(p,K,d,nrun);
sig2_beta_arr = zeros(K,d-1,nrun);
beta0_arr = zeros(p,d,nrun);
Q_mat_arr = zeros(p,K,nrun); gamma_q_arr = zeros(1,nrun);
Bern_K_arr = zeros(K, B, nrun);
tau_arr = zeros(B, nrun);
% subject-specific local parameters
A_mat_arr = zeros(n, K, nrun);
z_tau_mat_arr = zeros(n, B, nrun);
% CSP specific quantities begin
K_star_arr = zeros(1, nrun);
z_beta_vec_arr = zeros(K, nrun);


% Gibbs sampler begins here
for ii=1:nrun   

    % linear form inside the exp()
    % get Polya-Gamma parameters
    [PG_param, ~, C_ij_minus_c] = get_linear_form(beta_mat, beta0, Q_mat, A_mat);

    % (1) Gibbs sampling for PG varriables, beta_mat, and beta0    
    [beta_mat, beta0] = get_gibbs_beta(beta_mat, beta0, Y_arr, PG_param, C_ij_minus_c,...
        A_mat, Q_mat, sig2_beta, sig2_pseudo, var_beta0, mu_beta0);
    % --Alternative: if we want all beta to be positive, run the following:
    % [beta_mat, beta0] = get_gibbs_beta_allpos(beta_mat, beta0, Y_arr, PG_param, C_ij_minus_c,...
    %     A_mat, Q_mat, sig2_beta, sig2_pseudo, var_beta0, mu_beta0);

    % (2) Gibbs sampling for Q_mat
    % Q_mat = get_gibbs_Q(Q_mat, Y_arr, beta_mat, beta0, A_mat, gamma_q, sig2_beta, sig2_pseudo);
    Q_mat = get_gibbs_Q_noall0(Q_mat, Y_arr, beta_mat, beta0, A_mat, gamma_q, sig2_beta, sig2_pseudo);

    % (3) Gibbs sampling for gamma_q
    gamma_q = betarnd(1+sum(sum(Q_mat)), 1+p*K-sum(sum(Q_mat)), 1);
    
    % (4) Gibbs sampling for sig2_beta, size K * d-1
    for k=1:K
        if find(z_beta_mat(k,:)) > k 
            % slab
            for c=1:d-1 
                sig2_beta(k,c) = 1/gamrnd(a_sig + 0.5 * sum(Q_mat(:,k)), ...
                    1/(b_sig + 0.5 * sum(Q_mat(:,k) .* beta_mat(:,k,c).^2)) );
            end
        else
            % spike
            sig2_beta(k, :) = sig2_inf;
        end
    end
    
    % (5) Gibbs sampling for v_beta
    v_new_a = (sum(z_beta_mat(:,1:end-1), 1))';
    v_new_b = flip(cumsum([0; flip(v_new_a(2:end))]));
    v_beta(1:end-1) = betarnd(1 + v_new_a, alpha0 + v_new_b);
    % update w_beta according to v_beta
    w_beta = v_beta .* cumprod([1; 1-v_beta(1:end-1)]);
    
    % (6) Gibbs sampling for z_beta_mat
    for k=1:K
        
        log_normal_k = - 0.5/sig2_inf * sum(sum(Q_mat(:,k) .* squeeze(beta_mat(:,k,1:d-1)).^2)) ...
            - 0.5 * (sum(Q_mat(:,k)) * (d-1))* log(sig2_inf);
        log_t_k = -(a_sig+0.5) * sum(Q_mat(:,k)' * log( 1 + squeeze(beta_mat(:,k,1:d-1)).^2 / (2*b_sig) ));
        
        log_combine = [log_normal_k, log_t_k];
        log_combine_normalize = log_combine - max(log_combine);
        
        prob_kl = zeros(1, K);
        prob_kl(1:k) =  w_beta(1:k) * exp(log_combine_normalize(1));
        prob_kl((k+1):K) = w_beta((k+1):K) * exp(log_combine_normalize(2));
        prob_kl1 = prob_kl / sum(prob_kl);
        
        z_beta_mat(k,:) = mnrnd(1, prob_kl1);
    end
    
    
    % (7) Gibbs sampling for A_mat
    A_mat = get_gibbs_A(A_mat, Y_arr, beta_mat, beta0, Q_mat, Bern_K, z_tau_mat);
    
    
    % (8) Gibbs sampling for the deeper latent layer and Z
    % -- update Bern_K, deep tensor arms -- %
    Bern_K = betarnd(1 + A_mat' * z_tau_mat, 1 + (1-A_mat)' * z_tau_mat);
    % -- update tau, deep tensor core -- %
    tau_temp = gamrnd(1 + sum(z_tau_mat, 1), 1); tau = tau_temp/sum(tau_temp);
    % -- update z_tau_mat, deep tensor core membership -- %
    z_tau_mat = get_gibbs_z(tau, Bern_K, A_mat);
    

    %% store output data
    Q_mat_arr(:,:,ii) = Q_mat;
    gamma_q_arr(ii) = gamma_q;
    beta_mat_arr(:,:,:,ii) = beta_mat;
    beta0_arr(:,:,ii) = beta0;
    Bern_K_arr(:,:,ii) = Bern_K;
    tau_arr(:,ii) = tau';
    % subject-specific local parameters
    A_mat_arr(:,:,ii) = A_mat;
    z_tau_mat_arr(:,:,ii) = z_tau_mat;
    % beta variances
    sig2_beta_arr(:,:,ii) = sig2_beta;
    % CSP-specific quantities
    z_beta_vec_arr(:,ii) = z_beta_mat * (1:K)';
    K_star_arr(ii) = sum(z_beta_vec_arr(:,ii) > (1:K)');

    fprintf('Gibbs iteration %d completed\n', ii);
    

end






end