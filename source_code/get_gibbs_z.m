function [z_tau_mat] = get_gibbs_z(tau, Bern_K, A_mat)

n = size(A_mat,1);
B = length(tau);

z_tau_prob = zeros(n, B);
for b=1:B
    z_tau_prob(:, b) = tau(b) * ...
        prod( bsxfun(@power, (Bern_K(:,b))', A_mat), 2 ) .* ...
        prod( bsxfun(@power, (1-Bern_K(:,b))', 1-A_mat), 2 );
end
z_tau_prob = z_tau_prob ./ sum(z_tau_prob,2);
z_tau_mat = mnrnd(1, z_tau_prob);

end