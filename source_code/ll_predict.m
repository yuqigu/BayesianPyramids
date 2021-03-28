function [ll_pred] = ll_predict(Y_arr_test, Q_mat_arr, beta_mat_arr, beta0_arr, ...
    Bern_K_arr, tau_arr, nrun, burn, thin)

% only retain the effective samples
Q_mat_arr = Q_mat_arr(:, :, burn+1:thin:nrun);
beta0_arr = beta0_arr(:, :, burn+1:thin:nrun);
beta_mat_arr = beta_mat_arr(:, :, :, burn+1:thin:nrun);
Bern_K_arr = Bern_K_arr(:, :, burn+1:thin:nrun);
tau_arr = tau_arr(:, burn+1:thin:nrun);

n = size(Y_arr_test, 1);

%
nrun_eff = (nrun-burn)/thin;
ll_test_mat = zeros(n, nrun_eff);

for ii=1:nrun_eff
    ll_test_mat(:, ii) = loglik_marginal(Y_arr_test, Q_mat_arr(:,:,ii), beta_mat_arr(:,:,:,ii), beta0_arr(:,:,ii), ...
        Bern_K_arr(:,:,ii), tau_arr(:,ii));
end

ll_pred = sum(log(sum(exp(ll_test_mat), 2)/nrun_eff));

end