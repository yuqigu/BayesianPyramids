function [A_mat] = get_gibbs_A(A_mat, Y_arr, beta_mat, beta0, Q_mat, Bern_K, z_tau_mat)

[n, K] = size(A_mat);
B = size(z_tau_mat,2);

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
    % normalize_lf_ak1 = linear_form_ak1 - log(sum(exp(linear_form_ak1),3));
    normalize_lf_ak1 = linear_form_ak1 - max(linear_form_ak1,[],3);
    normalize_lf_ak0 = linear_form_ak0 - max(linear_form_ak0,[],3);

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
end