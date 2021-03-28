function [beta_mat, beta0] = get_gibbs_beta_uncons(beta_mat, beta0, Y_arr, PG_param, C_ij_minus_c,...
    A_mat, Q_mat, sig2_beta, sig2_pseudo, var_beta0, mu_beta0)

% unconstrained Gibbs sampling for beta

[n, p, d] = size(Y_arr);
W = zeros(n, p, d-1);

for c=1:d-1
    % -- update beta0 and beta_mat -- %
    % get a binary matrix for category c
    Y_binc = Y_arr(:,:,c);

    for j=1:p
        W(:,j,c) = pgdraw(PG_param(:,j,c));

        %%
        if any(Q_mat(j,:)==1)
        % size n_now * |K_j| + 1, |K_j| = sum(Q_mat(j,:)==1);
        A_mat_jeff = A_mat(:, Q_mat(j,:)==1);

        % |K_j| * |K_j|
        Sigma0_jc = diag(sig2_beta(Q_mat(j,:)==1, c));
        % mu0_jc = zeros(sum(Q_mat(j,:)),1);

        % |K_j| * |K_j|
        Sigma_jc = inv(inv(Sigma0_jc) + A_mat_jeff' * (W(:,j,c) .* A_mat_jeff));
        Sigma_jc = 0.5 * (Sigma_jc + Sigma_jc');

        % |K_j| * 1
        mu_jc = Sigma_jc * (A_mat_jeff' * (Y_binc(:,j) - 0.5 + W(:,j,c) .* (C_ij_minus_c(:,j,c) - beta0(j,c)) ));

        % Option 1: no truncation
        beta_mat(j, Q_mat(j,:)==1, c) = mvnrnd(mu_jc, Sigma_jc, 1);

%         % Option 2: truncate normal to be positive
%         num_k = sum(Q_mat(j,:)==1);
%         if num_k==1
%             % % option 2.1: sampling from truncated normal
%             % pd_beta =  makedist('Normal', 'mu', mu_jc, 'sigma', sqrt(Sigma_jc));
%             % beta_mat(j, Q_mat(j,:)==1, c) = random(truncate(pd_beta,0,Inf), 1);
% 
%             % option 2.2: more robust sampling from truncated normal
%             rv_unif = unifrnd(normcdf(0,mu_jc,sqrt(Sigma_jc)),1, 1);
%             beta_candidate = norminv(rv_unif, mu_jc, sqrt(Sigma_jc));
%             if ~isnan(beta_candidate)
%                 beta_mat(j, Q_mat(j,:)==1, c) = norminv(rv_unif, mu_jc, sqrt(Sigma_jc));
%             else
%                 beta_mat(j, Q_mat(j,:)==1, c) = 0;
%             end
% 
%         else
%             try
%                 beta_mat(j, Q_mat(j,:)==1, c) = mvnrnd(mu_jc, Sigma_jc, 1);
%             catch
%                 beta_mat(j, Q_mat(j,:)==1, c) = mvnrnd(mu_jc, diag(max(diag(Sigma_jc),0)), 1);
%             end
%         end
%         % Option 2 ends here


        if any(Q_mat(j,:)==0)
            % option 1: no truncation
            beta_mat(j, Q_mat(j,:)==0, c) = normrnd(0, sqrt(sig2_pseudo), [1 sum(Q_mat(j,:)==0)]);
            %
            % % option 2: sample from truncated distribution
            % pd_pseudo =  makedist('Normal', 'mu', 0, 'sigma', sqrt(sig2_pseudo));
            % beta_mat(j, Q_mat(j,:)==0, c) = random(truncate(pd_pseudo,0,Inf), [1 sum(Q_mat(j,:)==0)]);
        end

        %% beta0, intercept
        var_beta00 = inv( 1/var_beta0(c) + sum(W(:,j,c)) );

        mu_beta00 = var_beta00 * ( var_beta0(c)\mu_beta0(c) + ...
            sum( Y_binc(:,j)-1/2 - W(:,j,c) .* (A_mat * (Q_mat(j,:).*beta_mat(j,:,c))' - C_ij_minus_c(:,j,c)) ));

        beta0(j,c) = normrnd(mu_beta00, sqrt(var_beta00), 1);

        end
    end
end
    
end