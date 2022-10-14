addpath('./from_cluster')
% load splice_CSPbeta_alpha2_n3175p60_Kini7_d4;
% load splice1_overfitmix_n3175p60_Kini7_B5_d4;
load splice_alpha2_overfitmix_n3175p60_Kini7_B5_d4;

% % Option 1: permute
% [tau_arr_permu, B_permu_ind] = sort(tau_arr, 1);
% % note that tau_arr(B_permu_ind(:,ii), ii) is sorted from smallest to large
% Bern_K_arr_permu = zeros(K, B, nrun);
% for ii = 1:nrun
%     Bern_K_arr_permu(:,:,ii) = Bern_K_arr(:, B_permu_ind(:,ii), ii);
% end

% Option 2: do not permute
tau_arr_permu = tau_arr;
Bern_K_arr_permu = Bern_K_arr;


% --- traceplots --- %
figure; 
for b = 1:B
plot(tau_arr_permu(b,1:nrun))
hold on
end
% legend({'1', '2', '3'})
legend({'1', '2', '3', '4', '5'})
xlabel('Gibbs iterations');
% title('tau traceplot');
title('traceplot for deep tensor core $(\tau_{b})$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 18)
% print('-r300', 'splice_trace_tau_overfitmix', '-dpng')


figure
for k = 2:6
    plot(squeeze(Bern_K_arr_permu(k,1,:)))
    hold on
    plot(squeeze(Bern_K_arr_permu(k,2,:)))
    hold on
    plot(squeeze(Bern_K_arr_permu(k,3,:)))
    hold on
end
xlabel('Gibbs iterations');
title('traceplot for deep tensor arms $(\eta_{kb})$', 'Interpreter', 'Latex'); 
set(gca, 'FontSize', 18)
% print('-r300', 'splice_trace_deep_arm3', '-dpng')

beta0_pomean = mean(beta0_arr(:,:,burn+1:thin:nrun), 3);
beta0_pomean

% Bern_K posterior mean
Bern_K_pomean = mean(Bern_K_arr_permu(1:K,:,burn+1:thin:nrun), 3);
Bern_K_pomean

figure; imagesc(Bern_K_pomean); colorbar

% tau posterior mean
tau_pomean = mean(tau_arr(:,burn+1:thin:nrun), 2);
fprintf('\ntau posterior mean is:\n')
tau_pomean



% % Q post. mean
Q_mat_pomean_binary = (mean(Q_mat_arr(:,1:K-1,burn+1:thin:nrun), 3) > 0.5);


% z_tau_mat post. mean
z_mat_pomean = mean(z_tau_mat_arr(:,:,burn+1:thin:nrun), 3);

[~, z_tau_mat_int] = max(z_mat_pomean, [], 2);

% option 1: if B = 3
% % calculates the best permutation
% [z_tau_perm, best_perm] = get_best_permute(z_tau_mat_int, gene_type);
% option 2: if B > 3
z_tau_perm = z_tau_mat_int;


% 
% heatmap
x = repmat(1:B,B,1); % generate x-coordinates
y = x'; % generate y-coordinates
% % Generate Labels
% t = num2cell(acc_mat); % extact values into cells
% t = cellfun(@num2str, t, 'UniformOutput', false); % convert to string
% % Draw Image and Label Pixels
% figure; imagesc(acc_mat); colorbar
% pbaspect([1 1 1])
% text(x(:), y(:), t, 'HorizontalAlignment', 'Center', 'FontSize', 15)
% set(gca, 'FontSize', 15);
% xticks(1:B); yticks(1:B);
% xticklabels({'EI','IE','n'}); yticklabels({'EI','IE','n'})
% title('CSP prior, matching z to gene type')
% % print('-r300', 'splice_confmat_csp', '-dpng')


% % A_mat post. mean
% A_post_mean_prob = mean(A_mat_arr(:,1:K,burn+1:end), 3);
A_post_mean_prob = A_mat_pomean(:,1:K-1);
A_mat_pomean_bin = (A_post_mean_prob > 0.5);
% 
% [gene_type, (1:n)', A_mat_pomean];


% ------- Rule List preparation begins here ------- %
% conver z_tau_perm and A_mat_pomean to multivariate binary rule lists
z_tau_binmat = zeros(n, B);
rows = (1:n)';
cols = z_tau_perm;
lin_idx = sub2ind([n, B],rows,cols);
z_tau_binmat(lin_idx) = 1;

ZA = [z_tau_binmat, A_mat_pomean_bin];

% % % % -- write to csv file -- %
% csvwrite('Y1_splice_K_upper7_alpha3.csv', [ZA, gene_type==1])
% csvwrite('Y2_splice_K_upper7_alpha3.csv', [ZA, gene_type==2])
% csvwrite('Y3_splice_K_upper7_alpha3.csv', [ZA, gene_type==3])


K_star_arrmat = zeros(K, nrun);
K_star_arrmat(K_star_arr, :) = 1;
% figure; imagesc(K_star_arrmat); colorbar
tabulate(K_star_arr)
tabulate(K_star_arr(burn+1:thin:nrun))

K_select = mode(K_star_arr(burn+1:thin:nrun));
traits_select = find(z_beta_vec > (1:K)');


% --- only retain those effective latent traits --- %
ZA_eff = [z_tau_binmat, A_mat_pomean_bin(:, traits_select)];
csvwrite('Y1_splice_B5.csv', [ZA_eff, gene_type==1])
csvwrite('Y2_splice_B5.csv', [ZA_eff, gene_type==2])
csvwrite('Y3_splice_B5.csv', [ZA_eff, gene_type==3])

% % --- only retain those effective latent traits --- %

pink3 = pink(3);
% reversepink = [pink3(3,:); pink3(1,:)];
reversepink = [pink3(1,:); pink3(3,:)];


figure
ax1 = subplot(131);
imagesc(Q_mat_pomean_binary(:,traits_select));
xticks(1:length(traits_select))
xlabel('binary latent vars'); ylabel('p=60 locations'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old1 = colormap(ax1,reversepink); 
bar1 = colorbar; bar1.Ticks = [0.25 0.75]; bar1.TickLabels = {'0' '1'};
title('loci loading')
% second plot
ax2 = subplot(132);
imagesc(A_mat_pomean_bin(:,traits_select));
xticks(1:length(traits_select))
xlabel('binary latent vars'); ylabel('n=3175 sequences'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12)
old2 = colormap(ax2,reversepink); colormap( flipud(old2) ); 
bar2 = colorbar; bar2.Ticks = [0.25 0.75]; bar2.TickLabels = {'0' '1'};
title('1st layer traits')
% third plot
ax3 = subplot(133);
imagesc(z_tau_binmat);
xticks(1:B)
xlabel('z membership'); ylabel('n=3175 sequences'); pbaspect([1 1.5 1])
set(gca, 'FontSize', 12); 
old3 = colormap(ax3,reversepink); colormap( flipud(old3) );
bar3 = colorbar; bar3.Ticks = [0.25 0.75]; bar3.TickLabels = {'0' '1'};
title('2nd layer membership')

% print('-r300', 'splice_csp_Kupper_7_Kfinal_5', '-dpng')
% print('-r300', 'splice_csp_overfit_B', '-dpng')

% -- plot the held-out sequence types -- %
figure
imagesc(gene_type); 
colormap(parula(3)); 
% cbh3 = colorbar; set(cbh3,'YTick', 1:3)
title('type'); set(gca, 'FontSize', 12)
xticks(''); xticklabels({''}); xlabel('held out')
ylabel('n=3175 sequences');
% print('-r400', 'splice_type', '-dpng')



% % -- look at the the Q_pomean -- %
% [(1:p)', Q_mat_pomean(:,traits_select)]
% 
% p_set1 = 1:27;
% p_set2 = 28:37;
% p_set3 = 38:p;
% 
% figure
% subplot(1,3,1)
% imagesc(Q_mat_pomean(p_set1,traits_select));
% subplot(1,3,2)
% imagesc(Q_mat_pomean(p_set2,traits_select));
% subplot(1,3,3)
% imagesc(Q_mat_pomean(p_set3,traits_select));



% -- consider the results given by the rule lists approach -- %
% -- consider the results given by the rule lists approach -- %
A_hat_binary = A_mat_pomean_bin(:,traits_select);
G_hat_binary = Q_mat_pomean_binary(:,traits_select);
z_tau_binmat

t_hat_EI = A_hat_binary(:,1) .* A_hat_binary(:,5);
t_hat_IE = z_tau_binmat(:,2);
t_hat_N = (1-z_tau_binmat(:,2)) .* (1-A_hat_binary(:,5));

t_hat_mat = [t_hat_EI, t_hat_IE, t_hat_N];

t_hat_mat1 = [t_hat_EI, t_hat_IE .* (1-t_hat_EI), t_hat_N];
t_hat_mat2 = [t_hat_EI .* (1-t_hat_IE), t_hat_IE, t_hat_N];

tabulate(sum(t_hat_mat, 2))
tabulate(sum(t_hat_mat1, 2))
tabulate(sum(t_hat_mat2, 2))

% gene_type
conf_mat_hat1 = zeros(3,3);
conf_mat_hat2 = zeros(3,3);

% confusion mat

for bb=1:B
    for cc=1:B
        conf_mat_hat1(bb,cc) = sum(gene_type==bb & t_hat_mat1(:,cc)) / sum(gene_type==bb);
        conf_mat_hat2(bb,cc) = sum(gene_type==bb & t_hat_mat2(:,cc)) / sum(gene_type==bb);
    end
end

       
% Generate Labels
fun_4f = @(x) sprintf('%0.4f', x);

t1 = num2cell(conf_mat_hat1); % extact values into cells
t1 = cellfun(fun_4f, t1, 'UniformOutput', 0);
t1 = cellfun(@num2str, t1, 'UniformOutput', false); % convert to string

% Generate Labels
t2 = num2cell(conf_mat_hat2); % extact values into cells
t2 = cellfun(fun_4f, t2, 'UniformOutput', 0);
t2 = cellfun(@num2str, t2, 'UniformOutput', 0); % convert to string


% Draw Image and Label Pixels
figure; 
subplot(121); imagesc(conf_mat_hat1); pbaspect([1 1 1])
colormap('gray'); colorbar; caxis([0 1])
text(x(:), y(:), t1, 'HorizontalAlignment', 'Center', 'FontSize', 12, 'Color', [0.5 0.5 0.5])
set(gca, 'FontSize', 14)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','N'}); yticklabels({'EI','IE','N'})
title('confusion matrix $t^{\star}$', 'Interpreter', 'latex')
xlabel('label given by $t^{\star}$', 'Interpreter', 'latex')
ylabel('held-out sequence types')
%
subplot(122); imagesc(conf_mat_hat2); pbaspect([1 1 1])
colormap('gray'); colorbar; caxis([0 1])
text(x(:), y(:), t2, 'HorizontalAlignment', 'Center', 'FontSize', 12, 'Color', [0.5 0.5 0.5])
set(gca, 'FontSize', 14)
xticks(1:B); yticks(1:B);
xticklabels({'EI','IE','N'}); yticklabels({'EI','IE','N'})
title('confusion matrix $t^{\dagger}$', 'Interpreter', 'latex')
xlabel('label given by $t^{\dagger}$', 'Interpreter', 'latex')
ylabel('held-out sequence types')

% print('-r500', 'conf_mat2', '-dpng')