function [L,C] = kmeanspp_ini2(X,k)
% The k-means++ initialization.
% Input X : V*D (D data points in V dimension)
% Output L : D dimensional initial labels
% Output C : k*V inital k centers

C = X(:,1+round(rand*(size(X,2)-1)));
L = ones(1,size(X,2));
for i = 2:k
    D = X-C(:,L);
    D = cumsum(dot(D,D,1));
    if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
    C(:,i) = X(:,find(rand < D/D(end),1));
    % We are doing C'*X which has same dim as X (k*s), i.e. we assume we
    % have memory for full(X)
    [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
end
C = C';
end