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
for k=1:K
    plot(squeeze(Bern_K_arr(k,1,:)))
    hold on
    plot(squeeze(Bern_K_arr(k,2,:)))
end
xlabel('Gibbs iterations');
title('deep tensor arms traceplot'); set(gca, 'FontSize', 12)
% print('-r300', 'splice_trace_deep_arm', '-dpng')




beta0_pomean = mean(beta0_arr(:,:,burn+1:end), 3);
beta0_pomean

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr(1:K,:,burn+1:end), 3);
Bern_K_pomean

% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:end), 2);
fprintf('\ntau posterior mean is:\n')
tau_pomean



% % Q post. mean
Q_mat_pomean = (mean(Q_mat_arr(:,1:K,burn+1:end), 3) > 0.5);


fprintf('\nlocation loading structure:\n')
[(1:p)', Q_mat_pomean]


% z_tau_mat post. mean
% z_mat_arr = z_mat_arr;
z_mat_pomean = mean(z_mat_arr(:,:,burn+1:end), 3);

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
A_post_mean_prob = A_mat_pomean_rep;
A_mat_pomean = (A_post_mean_prob > 0.5);
% 
% [gene_type, (1:n)', A_mat_pomean];


% ------- Rule List preparation begins here ------- %
% conver z_tau_perm and A_mat_pomean to multivariate binary rule lists
z_tau_binmat = zeros(n, B);
rows = (1:n)';
cols = z_tau_perm;
lin_idx = sub2ind([n, B],rows,cols);
z_tau_binmat(lin_idx) = 1;

ZA = [z_tau_binmat, A_mat_pomean];




% 
% figure
% subplot(131)
% imagesc(Q_mat_pomean)
% % xticks(1:K_ini); yticks(1:p)
% xlabel('binary latent vars')
% ylabel('p=60 locations')
% set(gca, 'FontSize', 12)
% title('1st layer location loading structure')
% %
% subplot(132)
% imagesc(A_mat_pomean)
% % xticks(1:K_ini); yticks(1:3:n);
% xlabel('binary latent vars')
% ylabel('n=3175 genes')
% set(gca, 'FontSize', 12)
% title('1st layer gene grouping structure')
% subplot(133)
% imagesc(ZA)
% % xticks(1:K_ini); yticks(1:3:n);
% xlabel('binary latent vars')
% ylabel('n=3175 genes')
% set(gca, 'FontSize', 12)
% title('2 layers gene groupings combined')
% % print('-r300', 'splice_K4_QA', '-dpng')



%%%%%%
pink3 = pink(3);
reversepink = [pink3(3,:); pink3(1,:)];

figure
ax1 = subplot(141);
imagesc(Q_mat_pomean);
xlabel('binary latent vars'); ylabel('p=60 locations'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old1 = colormap(ax1,reversepink); 
bar1 = colorbar; bar1.Ticks = [0.25 0.75]; bar1.TickLabels = {'0' '1'};
title('loci loading')
% second plot
ax2 = subplot(142);
imagesc(A_mat_pomean); 
xlabel('binary latent vars'); ylabel('n=3175 genes'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old2 = colormap(ax2,reversepink); colormap( flipud(old2) ); 
bar2 = colorbar; bar2.Ticks = [0.25 0.75]; bar2.TickLabels = {'0' '1'};
title('1st layer gene traits')
% third plot
ax3 = subplot(143);
imagesc(z_tau_binmat);
xlabel('z membership'); ylabel('n=3175 genes'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12); 
old3 = colormap(ax3,reversepink); colormap( flipud(old3) ); 
bar3 = colorbar; bar3.Ticks = [0.25 0.75]; bar3.TickLabels = {'0' '1'};
title('2nd layer gene membership')
% forth plot
ax4 = subplot(144); 
imagesc(gene_type); xticks(1); xticklabels({''})
xlabel('held out');  ylabel('n=3175 genes')
title('gene types'); set(gca, 'FontSize', 12); pbaspect([1 7 1])
%
colormap(ax4, parula(3));
% bar = colorbar; bar.Ticks = [1.3 2 2.7]; bar.TickLabels = {'EI' 'IE' 'N'};
% print('-r300', 'splice_sglca_K3', '-dpng')


% figure
% subplot(141)
% imagesc(Q_mat_pomean_binary)
% % xticks(1:K_ini); yticks(1:p)
% xlabel('binary latent vars')
% ylabel('p=60 locations')
% set(gca, 'FontSize', 12)
% colorbar; colormap(parula)
% title('1st layer loci loading')
% %
% subplot(142)
% imagesc(A_mat_pomean_bin)
% % xticks(1:K_ini); yticks(1:3:n);
% xlabel('binary latent vars')
% ylabel('n=3175 genes')
% set(gca, 'FontSize', 12)
% colorbar; colormap(parula)
% title('1st layer gene clustering')
% %
% subplot(143)
% imagesc(z_tau_binmat)
% % xticks(1:K_ini); yticks(1:3:n);
% xlabel('z membership'); ylabel('n=3175 genes')
% set(gca, 'FontSize', 12); colorbar; colormap(parula)
% title('2nd layer gene membership')
% subplot(144); imagesc(gene_type); xticks(1)
% colormap(pink(3)); cbh3 = colorbar; set(cbh3,'YTick', 1:3)
% title('gene types'); set(gca, 'FontSize', 12)


% print('-r300', 'splice_csp_Kupper_4', '-dpng')




% % -- write to csv file -- %
% csvwrite('Y1_splice_K4.csv', [ZA, gene_type==1])
% csvwrite('Y2_splice_K4.csv', [ZA, gene_type==2])
% csvwrite('Y3_splice_K4.csv', [ZA, gene_type==3])
