function [Q_mat] = get_gibbs_Q(Q_mat, Y_arr, beta_mat, beta0, A_mat, gamma_q, sig2_beta, sig2_pseudo)

K = size(Q_mat,2);
[~, p, d] = size(Y_arr);


for k=1:K

    % size p*1
    Q_mat_k1 = Q_mat; Q_mat_k1(:, k) = 1;
    Q_mat_k0 = Q_mat; Q_mat_k0(:, k) = 0;

    % get linear forms, n_now * p * d
    [~, linear_form_k1, ~] = get_linear_form(beta_mat, beta0, Q_mat_k1, A_mat);
    [~, linear_form_k0, ~] = get_linear_form(beta_mat, beta0, Q_mat_k0, A_mat);

    %% collapse out the PG variables
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
    prod_part_h1 = prod( bsxfun(@power, cond_prob_h1, Y_arr), 3);
    prod_part_h0 = prod( bsxfun(@power, cond_prob_h0, Y_arr), 3);

    ratio_qk_01 = sig2_pseudo^(-(d-1)/2) / prod(sig2_beta(k,:).^(-1/2)) ...
        * exp(-1/2 * squeeze(beta_mat(:,k,1:d-1).^2) * (sig2_pseudo^(-1) - sig2_beta(k,:).^(-1))' ...
        + sum(log(prod_part_h0) - log(prod_part_h1), 1)');

    prob_Qk = 1./(1 + (1-gamma_q)/gamma_q * ratio_qk_01);

    Q_mat(:, k) = (rand(p,1) < prob_Qk);

end
    
end