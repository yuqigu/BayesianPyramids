function [ir] = get_ideal_resp(Q, Z)

% @param Q:   Q-matrix of size J * K
% @param Z:   latent attribute profiles of individuals; size N * K

% @return ir: binary ideal response matrix of size N * J

% This function is similar to get_Gamma(Q) but more general

[J, K] = size(Q);
ir = prod(bsxfun(@power, reshape(Z, [1 size(Z,1) K]), reshape(Q, [J 1 K])), 3)';

ir = ir';


end