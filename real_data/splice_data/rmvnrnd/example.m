m = [2 0];
S = [1 0.9; 0.9 1];
X = rmvnrnd(m,S,1000,[-eye(2);eye(2)], [Inf; Inf; Inf; Inf]);
figure; plot(X(:,1),X(:,2),'.')

