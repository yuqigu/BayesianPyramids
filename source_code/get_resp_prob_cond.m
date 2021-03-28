function [rp_cond] = get_resp_prob_cond(theta_mat, I_full)
% produce R-matrix for all rows denoted by I, storing proportions of 
% positive response for each response pattern in I

[J, C] = size(theta_mat);
theta_mat0 = 1 - theta_mat;

rp_cond = zeros(2^J, C);
for i = 1:size(I_full, 1)
    rp_cond(i, :) = prod(theta_mat0(I_full(i,:)==0,:),1) .* ...
        prod(theta_mat(I_full(i,:)==1,:),1);
end


end