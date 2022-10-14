% -- consider the results given by the rule-lists classifier -- %
% The rule-lists classification is carried out in Python
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
