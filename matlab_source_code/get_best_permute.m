function [z_tau_perm, best_perm] = get_best_permute(z_tau_vec, base_vec)

% calculates the best permutation that matches a clustering to truth

B_tau = max(base_vec);
n = length(z_tau_vec);

rows = (1:n)'; 

cols = z_tau_vec;
lin_idx = sub2ind([n,B_tau],rows,cols);
z_tau_mat = zeros(n, B_tau);
z_tau_mat(lin_idx) = 1;

cols2 = base_vec; 
lin_idx2 = sub2ind([n,B_tau],rows,cols2);
base_mat = zeros(n, B_tau);
base_mat(lin_idx2) = 1;

all_perms = perms(1:B_tau);
err_vec = zeros(size(all_perms,1), 1);

for jj=1:size(all_perms,1)
    err_vec(jj) = sum(sum( ( base_mat - z_tau_mat(:,all_perms(jj,:)) ).^2 )) / (n*B_tau);
end

best_perm0 = all_perms(err_vec == min(err_vec), :);
best_perm = best_perm0(1, :);
z_tau_perm = z_tau_mat(:, best_perm) * (1:B_tau)';

end