g = 1;
Bern_K_permute(1:4,:,g)

A_mat_pomean_rep(:,:,g)

tau_pomean_rep(:,g)

I_full_K0 = get_I(K0, K0);
nu_g = get_resp_prob_cond(Bern_K_permute(1:4,:,1), I_full_K0) * tau_pomean_rep(:,g);

% nu_true_new = get_resp_prob_cond(Bern_K_true(1:4,:,1), I_full_K0) * tau_true;

[nu_true, nu_g]

[I_full_K0, nu_true, nu_g]

% %% look at the response pattern proportions
% I_full_p = get_I(p, p);
% 
% prop_est_resp_pat = get_resp_prob_cond(Lambda_est_g, I_full_p) * nu_g;
% % ----


figure
plot(nu_true);
hold on
plot(nu_g)
legend({'true', 'est'})


figure
plot(sort(nu_true));
hold on
plot(sort(nu_g));
legend({'true', 'est'})