function[I] = get_I(D,J)
% question: if we set D=2, do we get all the two-way interations of
% reponses? Yes!
    
index0 = 2.^((J-1):-1:0)'; % length M
indices = index0;
for d = 2:D
    newindices = bsxfun(@plus, indices', ...
        index0(d:J)); % outer sum, size (M-d+1) * M
    indices = [indices ; newindices(:)];
    [~, ordertemp, ~] = unique(indices, 'first');
    indices = indices(sort(ordertemp));
end
I = binary(indices, J);

I = [zeros(1, J); I];

end
