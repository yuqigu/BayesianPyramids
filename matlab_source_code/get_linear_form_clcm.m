function [PG_param, linear_form, C_ij_minus_c] = get_linear_form_clcm(Z_mat, beta0, beta, S_mat)

n = size(Z_mat,1); [p, K] = size(S_mat); d = size(beta0,2);

% beta0 has size p * d
% beta has size p * d * K

% linear_form = zeros(n, p, d);

% size p * d * K
beta_arr_pdK = beta .* reshape(S_mat, [p 1 K]) + reshape(beta0, [p d 1]);

% size n * p * d, linear form
linear_form = sum(bsxfun(@times, reshape(beta_arr_pdK,[1 p d K]), reshape(Z_mat,[n 1 1 K])), 4);

% % n * p * (d-1)
% linear_form(:,:,1:d-1) = reshape(beta0(:,1:d-1), [1 p d-1]) + ...
%     reshape(beta_tilde_np, [n p 1]);

C_ij_minus_c = zeros(n, p, d);

for c=1:d
    C_ij_minus_c(:,:,c) = log(sum(exp(linear_form(:,:,[1:c-1,c+1:d])), 3));
end

% linear_form = linear_form_long(:,:,1:d-1);

PG_param = linear_form - C_ij_minus_c;


end