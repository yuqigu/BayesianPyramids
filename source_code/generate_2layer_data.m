clear; clc; close all
rng(0513)
% -- global parameters -- %

n = 2000;  
rep = 50;
K0 = 4;
d = 4;

Q_true = [eye(K0); eye(K0); eye(K0); 
    1 1 0 0; 0 1 1 0; 0 0 1 1; 1 0 0 1; 1 0 1 0; 0 1 0 1; ...
    1 1 1 0; 0 1 1 1];
p = size(Q_true, 1);


% % new beta_max
beta0_true = repmat([-4 -3 -2], [p 1]);

beta_max = [repmat([2 3 4], [12,1]); ...
            repmat([4 5 6], [6,1]); ...
            repmat([5 6 7], [2,1])];

beta_mat_true_dense = zeros(p, K0, d-1);
for j=1:p
    for c=1:d-1
        beta_mat_true_dense(j,:,c) = (beta_max(j,c)-beta0_true(j,c))/sum(Q_true(j,:));
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


%% independence
crv_true=zeros(p,p); crv_truevec=zeros(1,p*(p-1)/2); 

NMI_true=zeros(p,p); NMI_truevec=zeros(1,p*(p-1)/2);

ct_loop = 0;
for j1 = 1:p-1
    Lamj1 = squeeze(Lambda_arr_true(:,j1,:));
    pj1 = sum(bsxfun(@times,Lamj1,nu_true))'; Ij1 = - sum(pj1.*log(pj1)); 
    for j2 = j1+1:p
        ct_loop = ct_loop + 1;
        Lamj2 = squeeze(Lambda_arr_true(:,j2,:));
        pj2 = sum(bsxfun(@times,Lamj2,nu_true))'; Ij2 = - sum(pj2.*log(pj2));
        pj1j2 = bsxfun(@times,Lamj1,sqrt(nu_true))'*bsxfun(@times,Lamj2,sqrt(nu_true));
        
        tmp_MI = sum(sum(pj1j2.*log(pj1j2./(pj1*pj2')))); tmp_NMI = tmp_MI/sqrt(Ij1*Ij2);
        NMI_true(j1,j2) = tmp_NMI; NMI_true(j2,j1) = NMI_true(j1,j2);
        NMI_truevec(ct_loop) = tmp_NMI;
        crv = ((pj1j2 - pj1*pj2').^2)./(pj1*pj2'); tmp_crv = sqrt(sum(sum(crv/(d-1))));   
        crv_true(j1,j2) = tmp_crv; crv_true(j2,j1) = crv_true(j1,j2);
        crv_truevec(ct_loop) = tmp_crv;
    end
end

figure
subplot(121); imagesc(NMI_true); colorbar; title('NMI true'); pbaspect([1 1 1])
subplot(122); imagesc(crv_true); colorbar; title('Cramers V true'); pbaspect([1 1 1])



filename = strcat('datasets_weak_rep', num2str(rep), '_K', num2str(K0), 'd', num2str(d), ...
    '_n', num2str(n), 'p', num2str(p),  '.mat');


save(filename);



