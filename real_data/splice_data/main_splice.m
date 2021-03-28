% --- traceplots --- %
figure; 
plot(tau_arr(1,1:end))
hold on
plot(tau_arr(2,1:end))
hold on
plot(tau_arr(3,1:end))
title('tau traceplot'); set(gca, 'FontSize', 12)
% print('-r300', 'splice_tau_trace', '-dpng')

figure
plot(squeeze(Bern_K_arr(1,1,:)))
hold on
plot(squeeze(Bern_K_arr(1,2,:)))
hold on
plot(squeeze(Bern_K_arr(2,1,:)))
hold on
plot(squeeze(Bern_K_arr(2,2,:)))%
hold on
plot(squeeze(Bern_K_arr(3,1,:)))
hold on
plot(squeeze(Bern_K_arr(3,2,:)))
hold on
plot(squeeze(Bern_K_arr(4,1,:)))
hold on
plot(squeeze(Bern_K_arr(4,2,:)))
xlabel('Gibbs iterations');
title('deep tensor arms traceplot'); set(gca, 'FontSize', 12)
% print('-r300', 'splice_trace_deep_arm', '-dpng')


% filename = strcat('mdata0_n', num2str(n), '_p', num2str(p), '_K', num2str(K_ini), '_d', num2str(d), '.mat');
% save(filename)

beta_mat_tilde_arr = zeros(p, K, d, nrun);
for ii=1:nrun
    beta_mat_tilde_arr(:,:,:,ii) = beta_mat_arr(:,:,:,ii) .* Q_mat_arr(:,:,ii);
end


% -- posterior means -- %
% beta_matcat1_pomean = mean(beta_mat_tilde_catarr(:,1:K_ini,1,burn+1:end), 4);
% beta_matcat2_pomean = mean(beta_mat_tilde_catarr(:,1:K_ini,2,burn+1:end), 4);
% beta_matcat3_pomean = mean(beta_mat_tilde_catarr(:,1:K_ini,3,burn+1:end), 4);

beta_mat_pomean = zeros(p,K,d-1);
for c=1:d-1
    beta_mat_pomean(:,:,c) = mean(beta_mat_tilde_arr(:,1:K,c,burn+1:end), 4);
end

% % beta posterior mean
% figure
% subplot(131); imagesc(beta_mat_cat_pomean(:,:,1)); colorbar
% subplot(132); imagesc(beta_mat_cat_pomean(:,:,2)); colorbar
% subplot(133); imagesc(beta_mat_cat_pomean(:,:,3)); colorbar
% % print('-r300', 'beta_bindata_catest', '-dpng')

beta0_pomean = mean(beta0_arr(:,:,burn+1:end), 3);
beta0_pomean;

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(1:K,:,burn+1:end), 3);
Bern_K_pomean;

% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:end), 2);
fprintf('\ntau posterior mean is:\n')
tau_pomean




% % Q post. mean
% Q_mat_pomean = (mean(Q_mat_arr(:,1:K,burn+1:end), 3) > 0.5);

% Q_mat_arr(:,:,end) - Q_mat_arr(:,:,end-1)

fprintf('\nlocation loading structure:\n')
[(1:p)', Q_mat_pomean];


% z_tau_mat post. mean
z_mat_pomean = mean(z_tau_mat_arr(:,:,burn+1:end), 3);

[~, z_tau_mat_int] = max(z_mat_pomean, [], 2);

% calculates the best permutation
[z_tau_perm, best_perm] = get_best_permute(z_tau_mat_int, gene_type);


% figure for comparison
figure; 
subplot(141); imagesc(z_mat_pomean); colorbar
title('tau posterior mean'); set(gca, 'FontSize', 12)
%
subplot(142); imagesc(z_mat_pomean(:,best_perm))
colormap(parula(3)); cbh2 = colorbar; set(cbh2,'YTick', 1:3)
title('permuted tau post. mean'); set(gca, 'FontSize', 12)
%
subplot(143); imagesc(z_tau_perm)
colormap(parula(3)); cbh2 = colorbar; set(cbh2,'YTick', 1:3)
title('rounded tau'); set(gca, 'FontSize', 12)
%
subplot(144); imagesc(gene_type); 
colormap(parula(3)); cbh3 = colorbar; set(cbh3,'YTick', 1:3)
title('gene types'); set(gca, 'FontSize', 12)
% print('-r300', 'splice_cluster_ind_9310', '-dpng')


% clustering accuracy
z_and_gene = [z_mat_pomean(:, best_perm), gene_type];
conf_mat = zeros(B,B);

for bb=1:B
    for cc=1:B
        conf_mat(bb,cc) = sum(gene_type==bb & z_tau_perm==cc);
    end
end

acc_clust = mean(z_tau_perm == gene_type);

% confusion matrix
conf_mat
acc_mat = conf_mat ./ sum(conf_mat, 2)

fprintf('\nOverall clustering accuracy is %1.4f\n', acc_clust);

fprintf('\nClustering accuracies for the 3 classes are %1.4f, %1.4f, %1.4f\n', ...
    conf_mat(1,1)/sum(conf_mat(1,:)), conf_mat(2,2)/sum(conf_mat(1,:)),...
    conf_mat(3,3)/sum(conf_mat(3,:)));

% figure
% imagesc(acc_mat)
% colormap('pink'); colorbar
% title('clustering prob. matrix')
% set(gca, 'FontSize', 12)

% heatmap
x = repmat(1:B,B,1); % generate x-coordinates
y = x'; % generate y-coordinates
% Generate Labels
t = num2cell(acc_mat); % extact values into cells
t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% Draw Image and Label Pixels
figure; imagesc(acc_mat); colorbar %colormap('pink'); colorbar
pbaspect([1 1 1])
text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 15)
set(gca, 'FontSize', 15)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','n'}); yticklabels({'EI','IE','n'})
title('clustering prob. mat, splice data')
% print('-r300', 'splice_conf_mat', '-dpng')



% % A_mat post. mean
% A_post_mean_prob = mean(A_mat_arr(:,1:K,burn+1:end), 3);
% A_mat_pomean = (A_post_mean_prob > 0.5);
% 
% [gene_type, (1:n)', A_mat_pomean];


figure
subplot(121)
imagesc(Q_mat_pomean)
% xticks(1:K_ini); yticks(1:p)
xlabel('binary latent vars')
ylabel('p=60 locations')
set(gca, 'FontSize', 12)
title('1st layer location loading structure')
%
subplot(122)
imagesc(A_mat_pomean)
% xticks(1:K_ini); yticks(1:3:n);
xlabel('binary latent vars')
ylabel('n=3175 genes')
set(gca, 'FontSize', 12)
title('1st layer gene grouping structure')
% print('-r300', 'splice_K4_QA', '-dpng')


% ------- Rule List preparation begins here ------- %
% conver z_tau_perm and A_mat_pomean to multivariate binary rule lists
z_tau_binmat = zeros(n, B);
rows = (1:n)'; 
cols = z_tau_perm;
lin_idx = sub2ind([n, B],rows,cols);
z_tau_binmat(lin_idx) = 1;

ZA = [z_tau_binmat, A_mat_pomean];


% % -- write to csv file -- %
% csvwrite('Y1_splice_K4.csv', [ZA, gene_type==1])
% csvwrite('Y2_splice_K4.csv', [ZA, gene_type==2])
% csvwrite('Y3_splice_K4.csv', [ZA, gene_type==3])
