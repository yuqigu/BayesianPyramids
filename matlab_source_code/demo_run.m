% clear workspace
clear;

% set seed
rng(1)

% -- global parameters -- %
n = 500;  % sample size
K0 = 4; % number of binary latent variables
d = 4; % number of categories for each latent variable

% specify the graphical matrix between bottom two layers
Q_true = [eye(K0); eye(K0); eye(K0); 
    1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1; 1 0 1 0; 0 1 0 1; ...
    1 1 1 0; 0 1 1 1];
p = size(Q_true, 1); % number of observed variables in the bottom layer

% specify the continous parameters between bottom two layers
beta0_true = repmat([-4 -3 -2], [p 1]);
beta_max = [repmat([2 3 4], [12,1]); ...
            repmat([4 5 6], [6,1]); ...
            repmat([5 6 7], [2,1])];
beta_mat_true_dense = zeros(p, K0, d-1);
for j = 1:p
    for c = 1:d-1
        beta_mat_true_dense(j,:,c) = (beta_max(j,c)-beta0_true(j,c))/sum(Q_true(j,:));
    end
end
beta_mat_true = beta_mat_true_dense .* Q_true;

% Use beta and Q to define the lower layer conditional prob. tables
A_all = get_I(K0, K0);
Lambda_arr_true = zeros(2^K0, p, d);
for j = 1:p

    nume_temp = exp(beta0_true(j,:) + A_all * squeeze(beta_mat_true(j,:,:)));
    denom_temp = 1 + sum(nume_temp, 2);
    Lambda_arr_true(:,j,1:end-1) = nume_temp ./ denom_temp;
    
end
Lambda_arr_true(:,:,end) = 1 - sum(Lambda_arr_true(:,:,1:end-1),3);

% simulate a dataset Y
Y = zeros(n,p); % data matrix
Y_arr = zeros(n,d,p); % binarized data array with 1 and 0

% Generate the shallower latent layer
ind_len = size(A_all, 1);
% nu_true = generate_prop(K0);
B = 2;
[nu_true, Bern_K_true, tau_true] = generate_latent_prop(B, K0);
z_mat_true = mnrnd(1, nu_true, n);
z_vec_true = z_mat_true * (1:2^K0)';

% Generate the bottom data layer
for j = 1:p
    % generate multinomial random variables
    for h = 1:2^K0
        % size d*1
        mult_vec = reshape(Lambda_arr_true(h,j,:),1,[]);
        Y_arr(z_mat_true(:,h)==1,:,j) = mnrnd(1, mult_vec, sum(z_mat_true(:,h)));
        Y(z_mat_true(:,h)==1,j) = Y_arr(z_mat_true(:,h)==1,:,j) * (1:d)';
    end
end

% Y_data has size n * p, which is the observed categorical matrix
writematrix(Y, 'Y_data.csv') 

% Apply the Bayesian Pyramid method to the dataset 'Y_data.csv'
saved_file = bayes_pyramid('Y_data.csv');

load(saved_file)


