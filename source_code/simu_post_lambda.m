A_all = get_I(K0, K0);

A_all_bin = binary(0:(2^K0-1), K0);

eff_sample = (nrun-burn)/thin;
Lambda_est_effsample = zeros(2^K0, p, d, rep, eff_sample);
Lambda_est_rep = zeros(2^K0, p, d, rep);

for g=1:rep
%     for ss = 1:eff_sample
    for j=1:p
        nume_temp = exp(beta0_pomean_rep(j,:,g) + A_all_bin * squeeze(beta_mattil_pomean_rep(j,1:K0,:,g)));
        denom_temp = sum(nume_temp, 2);
        Lambda_est_rep(:,j,:,g) = nume_temp ./ denom_temp;
    end
%     end

end


sum(Lambda_est_rep, 3)

figure
imagesc(Lambda_est_rep(:,:,1,1)); colorbar