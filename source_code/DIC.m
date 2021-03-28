
Q_mat_pomean = (mean(Q_mat_arr(:,:,burn+1:end), 3) > 0.5);

beta0_cat_pomean = mean(beta0_catarr(:,:,burn+1:end), 3);

beta_mat_cat_pomean = zeros(p,K_ini+1,d);
for c=1:d-1
    beta_mat_cat_pomean(:,:,c) = mean(beta_mat_catarr(:,:,c,burn+1:end), 4);
end

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(:,:,burn+1:end), 3);
 
% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:end), 2);


% -- compute complete DIC -- %
% tau_arr
% Bern_K_arr
% beta0_catarr
% beta_mat_catarr
% Q_mat_arr
% 
% 
% tau_pomean
% Bern_K_pomean
% beta0_cat_pomean
% beta_mat_cat_pomean
% Q_mat_pomean

end_iter = nrun;

eff_sample = (end_iter-burn);

DIC_part1 = zeros(nrun, 1);
DIC_part2 = zeros(nrun, 1);

for rr=1:nrun
    DIC_part1(rr) = loglik(Y_arr, Q_mat_arr(:,:,rr), beta_mat_catarr(:,:,:,rr), beta0_catarr(:,:,rr), ...
        Bern_K_arr(:,:,rr), tau_arr(:,rr), A_mat_arr(:,:,rr), z_tau_arr(:,:,rr));
end


for rr=1:nrun
    DIC_part2(rr) = loglik(Y_arr, Q_mat_pomean, beta_mat_cat_pomean, beta0_cat_pomean, ...
        Bern_K_pomean, tau_pomean, A_mat_arr(:,:,rr), z_tau_arr(:,:,rr));
end

DIC_vec = - 4/eff_sample/n * DIC_part1 + 2/eff_sample/n * DIC_part2;

figure
plot(DIC_part1)
hold on
plot(DIC_part2)

figure
plot(DIC_vec)

DIC_complete = sum(DIC_vec(burn+1:end_iter));

DIC_complete


% save('splicedata_DIC.mat')


