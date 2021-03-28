% load splice_K3_iter8000
% load splice_K3_new_iter8000

A_mat_pomean;

% covariate_alt = [z_tau_perm, A_mat_pomean];
covariate_alt = [A_mat_pomean];
covariate_alt_std = (covariate_alt - mean(covariate_alt,1)) ./ sqrt(var(covariate_alt));

[conf_mat_reg, prob_3cat_mat, model_coef] = get_latent_reg(covariate_alt_std, gene_type);
% acc_clust_reg = mean(ind_3cat_mat == gene_type);
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




% --- figures ---- %
%
figure;
subplot(151); imagesc(z_mat_pomean); colorbar
title('tau posterior mean'); set(gca, 'FontSize', 12)
%
subplot(152); imagesc(z_mat_pomean(:,best_perm))
colormap(parula(3)); cbh2 = colorbar; set(cbh2,'YTick', 1:3)
title('rounded tau'); set(gca, 'FontSize', 12)
%
subplot(153); imagesc(z_tau_perm)
colormap(parula(3)); cbh2 = colorbar; set(cbh2,'YTick', 1:3)
title('rounded tau permuted'); set(gca, 'FontSize', 12)
%
subplot(154); imagesc(prob_3cat_mat * (1:B)')
title('regression prediction'); set(gca, 'FontSize', 12)
colormap(parula(3)); cbh3 = colorbar; set(cbh3,'YTick', 1:3)
%
subplot(155); imagesc(gene_type); 
colormap(parula(3)); cbh3 = colorbar; set(cbh3,'YTick', 1:3)
title('gene types'); set(gca, 'FontSize', 12)
%
% print('-r300', 'splice_K4_compare_ind', '-dpng')



% ---- how about replace prob. with labels ---- %
% ---- how about replace prob. with labels ---- %
covariates_prob = [z_tau_perm, A_mat_pomean];
covariates_prob_std = (covariates_prob - mean(covariates_prob,1)) ./ sqrt(var(covariates_prob));


[conf_mat_regprob, ~] = get_latent_reg(covariates_prob_std, gene_type);
acc_mat_regprob = conf_mat_regprob ./ sum(conf_mat_regprob, 2);
acc_clust_regprob = sum(diag(conf_mat_regprob))/n

% covariates_prob = [z_tau_mat_pomean(:,best_perm), A_post_mean_prob];
% conf_mat_regprob = get_latent_reg(covariates_prob, gene_type);
% acc_clust_regprob = sum(diag(conf_mat_regprob))/n

[acc_clust, acc_clust_reg, acc_clust_regprob]


% conf_mat ./ sum(conf_mat,2)
% 
% conf_mat_reg ./ sum(conf_mat_reg,2)
% 
% conf_mat_regprob ./ sum(conf_mat_regprob,2)

% heatmap
x = repmat(1:B,B,1); % generate x-coordinates
y = x'; % generate y-coordinates
% Generate Labels
t = num2cell(acc_mat); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
%
t_reg = num2cell(acc_mat_reg); % extact values into cells
t_reg = cellfun(@num2str, t_reg, 'UniformOutput', false); % convert to string
%
t_regprob = num2cell(acc_mat_regprob); % extact values into cells
t_regprob = cellfun(@num2str, t_regprob, 'UniformOutput', false); % convert to string


% Draw Image and Label Pixels
figure; 
subplot(131); imagesc(acc_mat); colorbar %colormap('pink'); colorbar
pbaspect([1 1 1])
text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 12)
set(gca, 'FontSize', 12)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','N'}); yticklabels({'EI','IE','N'})
title(strcat('clustering prob., overall acc. ', num2str(acc_clust)))
%
subplot(132); imagesc(acc_mat_reg); colorbar %colormap('pink'); colorbar
pbaspect([1 1 1])
text(x(:), y(:), t_reg, 'HorizontalAlignment', 'Center', 'FontSize', 12)
set(gca, 'FontSize', 12)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','N'}); yticklabels({'EI','IE','N'})
title(strcat('prediction prob., overall acc. ', num2str(acc_clust_reg)))
%
subplot(133); imagesc(acc_mat_regprob); colorbar %colormap('pink'); colorbar
pbaspect([1 1 1])
text(x(:), y(:), t_regprob, 'HorizontalAlignment', 'Center', 'FontSize', 12)
set(gca, 'FontSize', 12)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','N'}); yticklabels({'EI','IE','N'})
title(strcat('prediction prob., overall acc. ', num2str(acc_clust_regprob)))
% print('-r300', 'splice_K4_new_compare3_acc', '-dpng')

% save('splice_K5_iter8000.mat')

figure; imagesc(prob_3cat_mat); colorbar
