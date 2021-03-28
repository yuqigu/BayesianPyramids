function [PG_param, linear_form, C_ij_minus_c] = get_linear_form(beta_mat, beta0, Q_mat, A_mat)

n = size(A_mat,1); p = size(Q_mat,1); d = size(beta_mat,3);

linear_form = zeros(n, p, d);
C_ij_minus_c = zeros(n, p, d);

for c=1:d
    linear_form(:,:,c) = beta0(:, c)' + A_mat * (Q_mat .* beta_mat(:,:,c))';
end

for c=1:d
    C_ij_minus_c(:,:,c) = log(sum(exp(linear_form(:,:,[1:c-1,c+1:d])), 3));
end

% linear_form = linear_form_long(:,:,1:d-1);

PG_param = linear_form - C_ij_minus_c;

end