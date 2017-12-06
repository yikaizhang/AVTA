function [alpha,err]=frank_wolfe_2(p,A,epsilon)

%%*********Frank wolfe algorithm minimize eucidean distance(Greedy step)***
%Input:
%   p: mx1 vector The point which needs to be approximated
%   A: mxn matrix. Columns of A are points in convex hull.
%   epsilon: precision parameter
%Output:
%   alpha: nx1 vector. The weight of vertices to approximate p.
%   err: The l1-distance of p to the convex hull. 
%%***************************************************************

%% Initialize first iterate with closest points
[m,n]=size(A);
alpha=zeros(n,1);
diff=A-p(:,ones(n,1));
dis=sum(diff.^2,1);
min_ind=min(find(dis==min(dis)));
alpha(min_ind)=1;

p_p=A(:,min_ind);

gd=A'*(p_p-p);

mu=-1*min(gd);
lambda=gd+mu;
k=0;
gap=9999;

%% Keep iterating untill either p is epsilon close to convex hull or duality gap is small.
while gap>epsilon
    k=k+1;
    if(sqrt(sum((p-p_p).^2))<epsilon)
       break; % stop when p' is close to p
    end

    s=min(find(gd==min(gd)));
    gamma=(p_p-A(:,s))'*(p_p-p)/((p_p-A(:,s))'*(p_p-A(:,s))+eps);

    alpha=(1-gamma).*alpha;
    alpha(s)=alpha(s)+gamma;
    p_p=(1-gamma)*p_p+gamma*A(:,s);
    gd=A'*(p_p-p); % gradient
    lbd=gd-min(gd);
    gap=alpha'*lbd; %duality gap
end

err=sum(abs(p-p_p));

    