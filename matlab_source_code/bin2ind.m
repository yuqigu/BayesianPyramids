function [ind] = bin2ind(bin_vec)
% This function outputs a vector of positive integers x in binary form
%
% @param x     : column vector of base 10 integers
% @param k     : output vector length (optional)
%
% @return r    : matrix with the binary forms of x (by row)
%

ind = sum(bsxfun(@times, bin_vec, 2.^((size(bin_vec,2)-1):-1:0)),2);

end

