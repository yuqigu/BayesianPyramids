function [vec] = loglik_marginal(Y_i, Q_mat, beta_mat, beta0, Bern_K, tau)

% Y_i = squeeze(Y_arr(1, :, :)); A_i = A_mat(1, :); 
% z_tau_i = z_tau_mat(1, :);

[n, p, d] = size(Y_i);
[K, B] = size(Bern_K);

A_all = binary(0:(2^K-1), K);

linear_part = zeros(2^K, p, d);

% deep layer
deep_posi = prod(bsxfun(@power, reshape(Bern_K, [1 B K]), reshape(A_all, [2^K 1 K])), 3);
deep_zero = prod(bsxfun(@power, reshape(1-Bern_K, [1 B K]), reshape(1-A_all, [2^K 1 K])), 3);

% 2^K * 1
long_prop = (deep_posi .* deep_zero) * tau;

% linear part
for c=1:d
    linear_part(:, :, c) = beta0(:, c)' + A_all * (Q_mat .* beta_mat(:,:,c))';
end

linear_part_norm = linear_part - max(linear_part, [], 3);

% 2^K * p * d
prob_response = exp(linear_part_norm) ./ sum(exp(linear_part_norm), 3);

% n * 2^K
prob_cond = prod(prod(bsxfun(@power, reshape(prob_response,[1 2^K p d]), reshape(Y_i,[n 1 p d])), 4), 3);

vec = log(prob_cond * long_prop);



end