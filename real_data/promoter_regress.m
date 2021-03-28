A_mat_pomean;

% covariate_alt = [z_tau_perm, A_mat_pomean];
covariate_alt = [A_mat_pomean];
covariate_alt_std = (covariate_alt - mean(covariate_alt,1)) ./ sqrt(var(covariate_alt));

[conf_mat_reg, prob_3cat_mat, model_coef] = get_latent_reg(covariate_alt_std, is_promoter);

% supervised confusion matrix
conf_mat_reg;

acc_clust_reg = sum(diag(conf_mat_reg)) / n;
acc_mat_reg = conf_mat_reg ./ sum(conf_mat_reg, 2);
acc_mat_reg

fprintf('\nOverall prediction accuracy is %1.4f\n', acc_clust_reg);

fprintf('\nPredction accuracies for the 3 classes are %1.4f, %1.4f, %1.4f\n', ...
    conf_mat_reg(1,1)/sum(conf_mat_reg(1,:)), conf_mat_reg(2,2)/sum(conf_mat_reg(1,:)),...
    conf_mat_reg(3,3)/sum(conf_mat_reg(3,:)));


1 - [acc_clust, acc_clust_reg]