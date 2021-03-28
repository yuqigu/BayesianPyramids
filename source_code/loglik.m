function [vec] = ...
    loglik(Y_i, Q_mat, beta_mat_cat, beta0_cat, Bern_K, tau, A_i, z_tau_i)

% Y_i = squeeze(Y_arr(1, :, :)); A_i = A_mat(1, :); 
% z_tau_i = z_tau_mat(1, :);

[n, p, d] = size(Y_i);
linear_part = zeros(n, p, d);

% keyboard

for c=1:d
    linear_part(:, :, c) = beta0_cat(:, c)' + A_i * (Q_mat .* beta_mat_cat(:,:,c))';
end

% loglik_part1 = sum(sum(sum(Y_i .* linear_part))) - sum(sum(log( sum(exp(linear_part),3) )));

loglik_part1 = sum(sum(Y_i .* linear_part, 3), 2) - sum(log( sum(exp(linear_part), 3)), 2);

log_deep = log(tau') + (A_i * log(Bern_K) + (1-A_i) * log(1-Bern_K));

% loglik_part2 = sum(sum(z_tau_i .* log_deep));

loglik_part2 = sum(z_tau_i .* log_deep, 2);

vec = loglik_part1 + loglik_part2;

end