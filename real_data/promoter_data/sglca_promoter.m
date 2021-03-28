%function [] = sglca_promoter(K)
% clear; clc; close all;

% K=4; sglca_promoter

tic;
addpath('../../')

rng(0514)

% MCMC setup
nrun = 10000; burn = 4000; thin = 5; 
% nrun = 100; burn = 40; thin = 1; 


Y_with_type = csvread('promoters.csv');
Y_mat_d4 = Y_with_type(:, 2:end);
is_promoter = Y_with_type(:, 1);
Y_mat = Y_mat_d4;


[n, p] = size(Y_mat);
d = max(max(Y_mat));
Y_arr = zeros(n,p,d);
for c=1:d
    Y_arr(:,:,c) = (Y_mat==c);
end

% global
B = 2; 
tau_dir = 1;

% prior for intercept
% mu_beta0 = -3 * ones(3,1);
mu_beta0 = - rand(d-1,1);
var_beta0 = 10 * ones(d-1,1);


% inverse Gamma prior for hypervariances
a_sig = 2; b_sig = 2;

% sig2_pseudo = v0^2;
sig2_pseudo = 0.07;

tic;
addpath('../../')

rng(0514)


Y_with_type = csvread('promoters.csv');
Y_mat_d4 = Y_with_type(:, 2:end);
is_promoter = Y_with_type(:, 1);
Y_mat = Y_mat_d4;

[n, p] = size(Y_mat);
d = max(max(Y_mat));
Y_arr = zeros(n,p,d);
for c=1:d
    Y_arr(:,:,c) = (Y_mat==c);
end


        

    % -- initialize parameters -- %
    % -- initialize deep layer parameters -- %
    % tau_temp = gamrnd(1, 1, [1 B]); % ones(B,1)/B; % C * 1
    % tau = tau_temp/sum(tau_temp);
    tau = ones(1,B)/B;
    Bern_K = rand(K, B); 
    Bern_K_ini = Bern_K; tau_ini = tau;
    % -- initialize local deep Z -- %
    z_tau_mat = mnrnd(1, tau, n);
    
    % -- initialize local binary patterns A -- %
    A_mat = (rand(K, n) < Bern_K(:, z_tau_mat*(1:B)'))';

    % % Option 1: fix gamma_q
    gamma_q = 0.5;
    % % Option 2: give gamma_q a prior
    % gamma_q = betarnd(1,1, 1);
    % gamma_q_ini = gamma_q;

    % -- initialize Q matrix -- %
    Q_mat = (rand(p,K) < gamma_q');
    % make sure every row contains at least one entry of 1
    % for any all-zero row, uniformly sample one entry to make it =1 
    if any(sum(Q_mat,2)==0)
        Q_mat(sum(Q_mat,2)==0, :) = mnrnd(1, 1/K*ones(K,1), sum(sum(Q_mat,2)==0));
    end
    Q_ini = Q_mat;
    

    % -- initialize sig2_beta -- %
    % sig2_beta, size K * d-1
    sig2_beta = 1./(gamrnd(a_sig, 1/b_sig, [K, d-1]));
    sig2_beta_ini = sig2_beta;

    % -- initialize beta-mat -- %
    % the p variables for each (c,k) are a priori independent, so it has a
    % diagonal covariance structure as follows
    beta_mat = zeros(p, K, d);
    for c=1:d-1
        for k=1:K
            num_spike = sum(Q_mat(:,k)==0);
            
            if any(Q_mat(:,k)==1)
                pd1 = makedist('Normal', 'mu', 0, 'sigma', sqrt(sig2_beta(k,c)));
                beta_mat(Q_mat(:,k)==1, k, c) = random(truncate(pd1,0,Inf), p-num_spike, 1);
                % normrnd(0, sqrt(sig2_beta(k,c)), [p-num_spike 1]);
            end
            
            if any(Q_mat(:,k)==0)
                pd0 = makedist('Normal', 'mu', 0, 'sigma', sqrt(sig2_pseudo));
                beta_mat(Q_mat(:,k)==0, k, c) = random(truncate(pd0,0,Inf), num_spike, 1);
                % normrnd(0, sqrt(sig2_pseudo), [num_spike 1]);
            end
            
        end
    end
    beta_mat_ini = beta_mat;
    
    
    % initialize intercept beta0, size p * d -- %
    beta0 = zeros(p, d);
    beta0(:, 1:d-1) = repmat(mu_beta0', [p,1]);
        %rmvnrnd(mu_beta0, diag(var_beta0)/10, p, [-eye(d-1);eye(d-1)], [Inf*ones(d-1,1); zeros(d-1,1)]);
    beta0_ini = beta0;
    

    % -- Polya-Gamma random variables; no need to initialize
    W = zeros(n, p, d-1);

    

    % -- book-keeping of the data through MCMC iterations -- %
    beta_mat_arr = zeros(p,K,d,nrun);
    beta0_arr = zeros(p,d,nrun);
    Q_mat_arr = zeros(p,K,nrun);
    Bern_K_arr = zeros(K, B, nrun);
    tau_arr = zeros(B, nrun);
    % subject-specific local parameters
    A_mat_arr = zeros(n, K, nrun);
    z_mat_arr = zeros(n, B, nrun);
    sig2_beta_arr = zeros(K,d-1,nrun);
    
    % gibbs_acrate = zeros(p,K,nrun);
    
    fprintf('Initialization completed\n');
    

    %% Gibbs sampler begins here
    for ii=1:nrun   
        
        % linear form inside the exp()
        % get Polya-Gamma parameters
        [PG_param, ~, C_ij_minus_c] = get_linear_form(beta_mat, beta0, Q_mat, A_mat);
  
        % (1)-(2) Gibbs sampling for PG varriables, beta_mat, and beta0    
        for c=1:d-1
            % -- update beta0 and beta_mat -- %
            % get a binary matrix for category c
            Y_binc = Y_arr(:,:,c);
            
            for j=1:p
                W(:,j,c) = pgdraw(PG_param(:,j,c));
                
                
                %%
                if any(Q_mat(j,:)==1)
                % size n * |K_j| + 1, |K_j| = sum(Q_mat(j,:)==1);
                A_mat_jeff = A_mat(:, Q_mat(j,:)==1);

                % |K_j| * |K_j|
                Sigma0_jc = diag(sig2_beta(Q_mat(j,:)==1, c));
                % mu0_jc = zeros(sum(Q_mat(j,:)),1);

                % |K_j| * |K_j|
                Sigma_jc = inv(inv(Sigma0_jc) + A_mat_jeff' * (W(:,j,c) .* A_mat_jeff));
                Sigma_jc = 0.5 * (Sigma_jc + Sigma_jc');

                % |K_j| * 1
                mu_jc = Sigma_jc * (A_mat_jeff' * (Y_binc(:,j) - 0.5 + W(:,j,c) .* (C_ij_minus_c(:,j,c) - beta0(j,c)) ));

                % Option 1: no truncation
                beta_mat(j, Q_mat(j,:)==1, c) = mvnrnd(mu_jc, Sigma_jc, 1);
                
                                                
                end % end of "if any(Q_mat(j,:)==1)"
                
                if any(Q_mat(j,:)==0)
                    % sample from the pseudo prior
                    
                    % % option 2: unconstrained
                    beta_mat(j, Q_mat(j,:)==0, c) = normrnd(0, sqrt(sig2_pseudo), [1 sum(Q_mat(j,:)==0)]);
                end
                
                %% beta0
                var_beta00 = inv( 1/var_beta0(c) + sum(W(:,j,c)) );

                mu_beta00 = var_beta00 * ( var_beta0(c)\mu_beta0(c) + ...
                    sum( Y_binc(:,j)-1/2 - W(:,j,c) .* (A_mat * (Q_mat(j,:).*beta_mat(j,:,c))' - C_ij_minus_c(:,j,c)) ));
                
                % % option 1
                % pd_beta0 = makedist('Normal', 'mu', mu_beta00, 'sigma', sqrt(var_beta00));
                % beta0(j,c) = random(truncate(pd_beta0,-Inf,0), 1);
                
                % option 2: unconstrained
                beta0(j,c) = normrnd(mu_beta00, sqrt(var_beta00), 1);
            end
                   
        end

        % (3) Gibbs sampling for Q_mat
        for k=1:K

            % size p*1
            Q_mat_k1 = Q_mat; Q_mat_k1(:, k) = 1;
            Q_mat_k0 = Q_mat; Q_mat_k0(:, k) = 0;
            
            % get linear forms, n * p * d
            [~, linear_form_k1, ~] = get_linear_form(beta_mat, beta0, Q_mat_k1, A_mat);
            [~, linear_form_k0, ~] = get_linear_form(beta_mat, beta0, Q_mat_k0, A_mat);
            
            %% collapse out the PG variables

            normalize_lf_k1 = linear_form_k1 - log(sum(exp(linear_form_k1),3));
            normalize_lf_k0 = linear_form_k0 - log(sum(exp(linear_form_k0),3));
            
            % min(normalize_lf_k1,[],3)
            
            % size n * p
            denom_category_h1 = sum(exp(normalize_lf_k1), 3);
            denom_category_h0 = sum(exp(normalize_lf_k0), 3);

            % size n * p * d; sum(cond_prob_h1, 3) = 1
            cond_prob_h1 = bsxfun(@rdivide, exp(normalize_lf_k1), denom_category_h1);
            cond_prob_h0 = bsxfun(@rdivide, exp(normalize_lf_k0), denom_category_h0);

            % size n * p
            prod_part_h1 = prod( bsxfun(@power, cond_prob_h1, Y_arr), 3);
            prod_part_h0 = prod( bsxfun(@power, cond_prob_h0, Y_arr), 3);
           
            ratio_qk_01 = sig2_pseudo^(-(d-1)/2) / prod(sig2_beta(k,:).^(-1/2)) ...
                * exp(-1/2 * squeeze(beta_mat(:,k,1:d-1).^2) * (sig2_pseudo^(-1) - sig2_beta(k,:).^(-1))' ...
                + sum(log(prod_part_h0) - log(prod_part_h1), 1)');

            prob_Qk = 1./(1 + (1-gamma_q)/gamma_q * ratio_qk_01);

            Q_mat(:, k) = (rand(p,1) < prob_Qk);

        end


        % (6) Gibbs sampling for A_mat
        % -- update A_mat -- %
        for k=1:K
            % size n*1; conditional prob. of a_{ik}=1 given z_{i}=1:B
            z_tau_vec = z_tau_mat*(1:B)';
            Bern_part_k1 = (Bern_K(k, z_tau_vec))';
            Bern_part_k0 = 1 - Bern_part_k1;

            A_mat_k1 = A_mat; A_mat_k1(:, k) = 1;
            A_mat_k0 = A_mat; A_mat_k0(:, k) = 0;

            % get linear forms, n * p * d
            [~, linear_form_ak1, ~] = get_linear_form(beta_mat, beta0, Q_mat, A_mat_k1);
            [~, linear_form_ak0, ~] = get_linear_form(beta_mat, beta0, Q_mat, A_mat_k0);
            
            %% collapse out the PG variables
            normalize_lf_ak1 = linear_form_ak1;
            normalize_lf_ak0 = linear_form_ak0;
            
            normalize_lf_ak1 = linear_form_ak1 - log(sum(exp(linear_form_ak1),3));
            normalize_lf_ak0 = linear_form_ak0 - log(sum(exp(linear_form_ak0),3));
            
            % size n * p
            denom_category_ah1 = sum(exp(normalize_lf_ak1), 3);
            denom_category_ah0 = sum(exp(normalize_lf_ak0), 3);

            % size n * p * d; sum(cond_prob_ah1, 3) = 1
            cond_prob_ah1 = bsxfun(@rdivide, exp(normalize_lf_ak1), denom_category_ah1);
            cond_prob_ah0 = bsxfun(@rdivide, exp(normalize_lf_ak0), denom_category_ah0);

            % size n * p
            prod_part_ah1 = prod( bsxfun(@power, cond_prob_ah1, Y_arr), 3);
            prod_part_ah0 = prod( bsxfun(@power, cond_prob_ah0, Y_arr), 3);
            
            ratio_alikelihood = exp(sum(log(prod_part_ah0) - log(prod_part_ah1), 2));
            
            prob_Ak = 1./( 1 + (Bern_part_k0./Bern_part_k1) .*  ratio_alikelihood);

            A_mat(:, k) = (rand(n,1) < prob_Ak);
        end


        
        % (4) Gibbs sampling for sig2_beta, size K * d-1
        % -- update sig2_beta -- %
        for c=1:d-1 
            for k=1:K
                sig2_beta(k,c) = 1/gamrnd(a_sig + 0.5 * sum(Q_mat(:,k)), ...
                    1/(b_sig + 0.5 * sum(Q_mat(:,k) .* beta_mat(:,k,c).^2)) );
            end
        end
        
        % (5) Gibbs sampling for the deeper latent layer and Z
        % -- update Bern_K, deep tensor arms -- %
        Bern_K = betarnd(1 + A_mat' * z_tau_mat, 1 + (1-A_mat)' * z_tau_mat);
        % -- update tau, deep tensor core -- %
        tau_temp = gamrnd(1 + sum(z_tau_mat, 1), 1); 
        tau = tau_temp/sum(tau_temp);
        % -- update z_tau_mat, deep tensor core membership -- %
        z_tau_prob = zeros(n, B);
        for b=1:B
            z_tau_prob(:, b) = tau(b) * ...
                prod( bsxfun(@power, (Bern_K(:,b))', A_mat), 2 ) .* ...
                prod( bsxfun(@power, (1-Bern_K(:,b))', 1-A_mat), 2 );
                % prod(beta(1 + A_mat .* z_tau_mat(:,b), 1 + (1-A_mat) .* z_tau_mat(:,b)), 2);
        end
        z_tau_prob = z_tau_prob ./ sum(z_tau_prob,2);
        z_tau_mat = mnrnd(1, z_tau_prob);
        
        
        %% store output data
        Q_mat_arr(:,:,ii) = Q_mat;
        beta_mat_arr(:,:,:,ii) = beta_mat;
        beta0_arr(:,:,ii) = beta0;
        Bern_K_arr(:,:,ii) = Bern_K;
        tau_arr(:,ii) = tau';
        % subject-specific local parameters
        A_mat_arr(:,:,ii) = A_mat;
        z_mat_arr(:,:,ii) = z_tau_mat;
        
        sig2_beta_arr(:,:,ii) = sig2_beta;

        fprintf('Gibbs iteration %d completed\n', ii);

    end
        
    % get beta_tilde
    beta_mat_tilde_arr = zeros(p, K, d, nrun);
    for ii=1:nrun
        beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
    end
    
    % Q matrix posterior mean
    Q_pm_raw = (mean(Q_mat_arr(:,:,burn+1:thin:end), 3) > 0.5);
    
    % find the column permutation starts
    Q_firstK_binind = 2.^(0:K-1) * Q_pm_raw(1:K, :);
    [~, permute_vec] = sort(Q_firstK_binind);
    
    % permute everything
    Q_mat_arr = Q_mat_arr(:,permute_vec,:);
    beta_mat_arr = beta_mat_arr(:,permute_vec,:,:);
    beta_mat_tilde_arr = beta_mat_tilde_arr(:,permute_vec,:,:);
    A_mat_arr = A_mat_arr(:,permute_vec,:);
    Bern_K_arr = Bern_K_arr(permute_vec,:,:);
    
    
    % store other posterior means
    % Q-matrix, permuted
    Q_mat_rep = Q_pm_raw(:, permute_vec);
    
    
    % beta_mat, permuted
    beta_mat_pomean_rep = mean(beta_mat_arr, 4);
    
    beta_mattil_pomean_rep = mean(beta_mat_tilde_arr, 4);
    
    
    % beta0 posterior mean
    beta0_pomean_rep = mean(beta0_arr(:,:,burn+1:thin:end), 3);
    
    % Bern_K posterior mean
    Bern_K_pomean_rep = mean(Bern_K_arr(:,:,burn+1:thin:end), 3);

    % tau posterior mean
    tau_pomean_rep = mean(tau_arr(:,burn+1:thin:end), 2);
    
    % A_mat and z posterior mean; A_mat permuted
    A_mat_pomean_rep = mean(A_mat_arr(:,:,burn+1:thin:end), 3);
    z_mat_pomean_rep = mean(z_mat_arr(:,:,burn+1:thin:end), 3);
    
    permute_vec_rep = permute_vec;
    
    
%     %% compute the DIC value
%     eff_sample = (nrun-burn)/thin;
% 
%     DIC_part1 = zeros(nrun, 1);
%     DIC_part2 = zeros(nrun, 1);
% 
%     for rr=1:nrun
%         DIC_part1(rr) = loglik(Y_arr, Q_mat_arr(:,:,rr), beta_mat_arr(:,:,:,rr), beta0_arr(:,:,rr), ...
%             Bern_K_arr(:,:,rr), tau_arr(:,rr), A_mat_arr(:,:,rr), z_mat_arr(:,:,rr));
%     end
%  
%  
%     for rr=1:nrun
%         DIC_part2(rr) = loglik(Y_arr, Q_mat_rep, beta_mat_pomean_rep, beta0_pomean_rep, ...
%              Bern_K_pomean_rep, tau_pomean_rep, A_mat_arr(:,:,rr), z_mat_arr(:,:,rr));
%     end
% 
%     DIC_vec = - 4/eff_sample/n * DIC_part1 + 2/eff_sample/n * DIC_part2;
%     DIC_rep = mean(DIC_vec(burn+1:thin:nrun));
    
%% compute the DIC value
eff_nrun = (nrun-burn)/thin;
loglik_data_eff = zeros(n, eff_nrun);
penalty_eff = zeros(n, eff_nrun);

for ii=burn+1:thin:nrun
    ind = ceil((ii-burn)/thin);
    loglik_data_eff(:,ind) = loglik(Y_arr, Q_mat_arr(:,:,ind), beta_mat_arr(:,:,:,ind), beta0_arr(:,:,ind), ...
        Bern_K_arr(:,:,ind), tau_arr(:,ind), A_mat_arr(:,:,ind), z_mat_arr(:,:,ind));

    Q_rep_g = (Q_mat_rep > 0.5);
    penalty_eff(:,ind) = loglik(Y_arr, Q_mat_rep, beta_mat_pomean_rep, beta0_pomean_rep, ...
        Bern_K_pomean_rep, tau_pomean_rep, A_mat_arr(:,:,ind), z_mat_arr(:,:,ind));
end

DIC = - 4/eff_nrun * sum(sum(loglik_data_eff)) ...
    + 2/eff_nrun * sum(sum(penalty_eff));


    
    
%     % new DIC
%     DIC_part3 = zeros(nrun, 1);
%     for rr=1:nrun
%         DIC_part3(rr) = loglik(Y_arr, Q_mat_rep, beta_mat_pomean_rep, beta0_pomean_rep, ...
%              Bern_K_pomean_rep, tau_pomean_rep, A_mat_arr(:,:,rr), z_mat_arr(:,:,rr));
%     end
%     DIC_vecnew = - 4/eff_sample/n * DIC_part1 + 2/eff_sample/n * DIC_part3;
%     DIC_repnew = mean(DIC_vecnew(burn+1:thin:nrun));
%     DIC_repnew
    
    
mytime = toc;


filename = strcat('mdata_promoter_K', num2str(K), 'd', num2str(d), ...
    '_n', num2str(n), 'p', num2str(p), '.mat');

save(filename);



%end