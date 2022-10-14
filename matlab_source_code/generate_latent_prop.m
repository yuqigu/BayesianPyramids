function [prop_vec_true, Bern_K_true, tau_true] = generate_latent_prop(B, K0)

% This function generates the probability tensor for the joint distribution
% of the K binary latent variables

% rng(0513); B = 2; K0 = 4;
tau_true = ones(B,1)/B; % C * 1

% Bern_K_true = repmat([0.2 0.8; 0.8 0.2; 0.8 0.2], [K/3, 1]); % K * C; Sparsity matrix
% Bern_K_true = rand(K0, B);

Bern_K_square = 0.6 * eye(B) + 0.2 * ones(B,B);
Bern_K_true_temp = repmat(Bern_K_square, [ceil(K0/B), 1]);
Bern_K_true = Bern_K_true_temp(1:K0, :);

I_full_K = get_I(K0, K0);
prob_cond_K = get_resp_prob_cond(Bern_K_true, I_full_K);

prop_vec_true = prob_cond_K * tau_true;

% get_MI(prob_cond_K, tau_true)

end