function [r] = binary(x, varargin)
% This function outputs a vector of positive integers x in binary form
%
% @param x     : column vector of base 10 integers
% @param k     : output vector length (optional)
%
% @return r    : matrix with the binary forms of x (by row)
%

optargin = size(varargin,2);
switch optargin
    case 0 % no extra input
        base = 2;
        k = max(floor(log2(x))+1,1);
    case 1 % input k
        base = 2;
        kmax = max(floor(log2(max(x)))+1,1);
        k = varargin{1};
        assert(k>=kmax);
end

x = reshape(x,[],1);
divs = floor(bsxfun(@rdivide, x,base.^((k-1):-1:0)));
r = divs - [zeros(length(x),1), base*divs(:,1:(k-1))];

end

