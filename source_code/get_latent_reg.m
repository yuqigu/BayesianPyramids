function [conf_mat, prob_3cat_mat, model_coef] = get_latent_reg(covariates, gene_type)

n = length(gene_type); B = length(unique(gene_type));

[model_coef, ~, ~] = mnrfit(covariates, categorical(gene_type) );

prob_1_nume = exp(covariates * model_coef(2:end,1) + model_coef(1,1));
prob_2_nume = exp(covariates * model_coef(2:end,2) + model_coef(1,2));

denom = sum(1 + prob_1_nume + prob_2_nume, 2);

prob_3cat_mat = [prob_1_nume, prob_2_nume, ones(n,1)] ./ denom;


[~, ind_3cat_mat] = max(prob_3cat_mat, [], 2);


% clustering accuracy
conf_mat = zeros(B,B);
for bb=1:B
    for cc=1:B
        conf_mat(bb,cc) = sum(gene_type==bb & ind_3cat_mat==cc);
    end
end


end