clear; clc; close all
rng(0513)
% -- global parameters -- %
% This function generate 2-latent-layer data under weaker signal strengths

n = 1750;  
rep = 50;
K0 = 4;
d = 4;

Q_true = [eye(K0); eye(K0); eye(K0); 
    1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1; 1 0 1 0; 0 1 0 1; ...
    1 1 1 0; 0 1 1 1];
p = size(Q_true, 1);

beta0_true = repmat([-3 -2 -1], [p 1]);

beta_max = [3 * ones(12,3); ...
            4 * ones(6, 3); ...
            6 * ones(2, 3)];

beta_mat_true_dense = zeros(p, K0, d-1);
for j=1:p
    for c=1:d-1
        beta_mat_true_dense(j,:,c) = (beta_max(j,c))/sum(Q_true(j,:));
    end
end
beta_mat_true = beta_mat_true_dense .* Q_true;


% generate data Y
Yn_big = zeros(n,p,rep); % data matrix
Yn_big_arr = zeros(n,d,p,rep); % binarized data array with 1 and 0

% Lambda, lower layer conditional prob. tables
% A_all = binary(0:(2^K0-1), K0);
A_all = get_I(K0, K0);
Lambda_arr_true = zeros(2^K0, p, d);
for j=1:p

    nume_temp = exp(beta0_true(j,:) + A_all * squeeze(beta_mat_true(j,:,:)));
    denom_temp = 1 + sum(nume_temp, 2);
    Lambda_arr_true(:,j,1:end-1) = nume_temp ./ denom_temp;
    
end
Lambda_arr_true(:,:,end) = 1 - sum(Lambda_arr_true(:,:,1:end-1),3);


%% proportions
% generate latent attribute patterns
ind_len = size(A_all, 1);
% nu_true = generate_prop(K0);
B = 2;
[nu_true, Bern_K_true, tau_true] = generate_latent_prop(B, K0);
z_mat_true = mnrnd(1, nu_true, n);
z_vec_true = z_mat_true * (1:2^K0)';


for g = 1:rep
    for j = 1:p
        % generate multinomial random variables
        for h=1:2^K0
            % size d*1
            mult_vec = reshape(Lambda_arr_true(h,j,:),1,[]);
            Yn_big_arr(z_mat_true(:,h)==1,:,j,g) = mnrnd(1, mult_vec, sum(z_mat_true(:,h)));
            Yn_big(z_mat_true(:,h)==1,j,g) = Yn_big_arr(z_mat_true(:,h)==1,:,j,g) * (1:d)';
        end
    end
end


filename = strcat('datasets_weak_rep', num2str(rep), '_K', num2str(K0), 'd', num2str(d), ...
    '_n', num2str(n), 'p', num2str(p),  '.mat');


save(filename);



