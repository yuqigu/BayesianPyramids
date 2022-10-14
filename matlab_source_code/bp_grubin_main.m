%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('mdata_rev_bp_grubin_new.mat')



fprintf('Median Gelman-Rubin statistic for beta-intercepts:\n')
med_beta0 = quantile(squeeze(grubin_beta0), 0.5, 3)

fprintf('Median Gelman-Rubin statistic for beta-main-effects:\n')
med_betamat = quantile(squeeze(grubin_betamat), 0.5, 4)

for c = 1:d-1
    [med_beta0(:,c), med_betamat(:,:,c)]
end

% fprintf('Median Gelman-Rubin statistic for beta-coefficients:\n')
% med_beta_all = [quantile(squeeze(grubin_beta0), 0.5, 2:3), ...
%         quantile(squeeze(grubin_betamat), 0.5, 3:4)]

fprintf('Median Gelman-Rubin statistic for deep conditional prob. tables:\n')
med_bern = quantile(squeeze(grubin_Bern_K), 0.5, 3)

fprintf('Median Gelman-Rubin statistic for deep latent proportions tau:\n')
med_tau = quantile(squeeze(grubin_tau), 0.5, 2)


% %%%
% quantile(squeeze(grubin_tau(1,1,:)), 0.75)
% 
% quantile(squeeze(grubin_beta0), 0.75, 3)
% 
% quantile(squeeze(grubin_betamat), 0.75, 3:4)
% 
% quantile(squeeze(grubin_Bern_K), 0.75, 3)