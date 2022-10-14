function [beta, beta0] = get_clcm_gibbs_logit_beta(beta, beta0, Y_arr, ...
        PG_param, C_ij_minus_c, Z_mat, S_mat)

K = size(S_mat,2);
[n, p, d] = size(Y_arr);
W = zeros(n, p, d-1);

sig2_spike = 0.01;
sig2_slab = 10;

var_beta0 = 10 * ones(d-1,1);
    
mu_beta0 = [-5; -3; -1];

for c = 1:d-1
    % -- update beta0 and beta -- %
    % get a binary matrix for category c
    Y_binc = Y_arr(:,:,c);

    for j=1:p
        W(:,j,c) = pgdraw(PG_param(:,j,c));

        %% beta
        if any(S_mat(j,:)==0)
            beta(j, c, S_mat(j,:)==0) = normrnd(0, sqrt(sig2_spike), [1 sum(S_mat(j,:)==0)]);
        end
                
        for k = 1:K
            if S_mat(j,k) == 1
                sig2_jc = inv(inv(sig2_slab) + Z_mat(:,k)' * (W(:,j,c) .* Z_mat(:,k)));
                x_jc = C_ij_minus_c(:,j,c) + diag(1./W(:,j,c)) * (Y_binc(:,j)-1/2);
                % size n * 1
                z_beta_prod =  Z_mat(:,[1:k-1 k+1:K]) * (squeeze(beta(j,c,[1:k-1 k+1:K])) .* S_mat(j,[1:k-1 k+1:K])');
                mu_jc = sig2_jc * ( (Z_mat(:,k) .* W(:,j,c))' * (x_jc - beta0(j,c) - z_beta_prod) );
                
%                 pd_beta =  makedist('Normal', 'mu', mu_jc, 'sigma', sqrt(sig2_jc));
%                 beta(j, c, k) = random(truncate(pd_beta,0,Inf), 1);
                beta(j, c, k) = normrnd(mu_jc, sqrt(sig2_jc), 1);
            end
        end

%         A_mat_jeff = A_mat(:, Q_mat(j,:)==1);
% 
%         % |K_j| * |K_j|
%         Sigma0_jc = diag(sig2_beta(Q_mat(j,:)==1, c));
%         % mu0_jc = zeros(sum(Q_mat(j,:)),1);
% 
%         % |K_j| * |K_j|
%         Sigma_jc = inv(inv(Sigma0_jc) + A_mat_jeff' * (W(:,j,c) .* A_mat_jeff));
%         Sigma_jc = 0.5 * (Sigma_jc + Sigma_jc');
% 
%         % |K_j| * 1
%         mu_jc = Sigma_jc * (A_mat_jeff' * (Y_binc(:,j) - 0.5 + W(:,j,c) .* (C_ij_minus_c(:,j,c) - beta0(j,c)) ));
                

        %% beta0, intercept
        var_beta00 = inv( 1/var_beta0(c) + sum(W(:,j,c)) );

        linear_part = C_ij_minus_c(:,j,c) - Z_mat * (squeeze(beta(j,c,:)) .* S_mat(j,:)'); % n * 1
        mu_beta00 = var_beta00 * ( var_beta0(c)\mu_beta0(c) + sum( Y_binc(:,j)-1/2 + W(:,j,c) .*  linear_part) );
        
%         %% beta0, intercept
%         var_beta00 = inv( 1/var_beta0(c) + sum(W(:,j,c)) );
% 
%         mu_beta00 = var_beta00 * ( var_beta0(c)\mu_beta0(c) + ...
%             sum( Y_binc(:,j)-1/2 - W(:,j,c) .* (A_mat * (Q_mat(j,:).*beta(j,:,c))' - C_ij_minus_c(:,j,c)) ));

        beta0(j,c) = normrnd(mu_beta00, sqrt(var_beta00), 1);

     end
end
    
end