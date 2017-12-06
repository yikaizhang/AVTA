function [inorout,p_prime,alpha_coe,dist]=Tri_Algo(dataset,p,epsilon,varargin)

%%*********** Triangle algorithm********
%Input: 
%   dataset: dim x K  matrix
%   p: Query point
%   epsilon : precision parameter
%   alpha_0: If warm up initialization is used input K x 1 vector otherwise leave blank. 
%Output:
%   inorout: 1 if p is in the convex hull, 0 other wise
%   p_prime: An approximation of p or witness of p
%   alpha_coe:K x 1 vector. Weight of query point using convex combination of colums of dataset
%   dist: The distance between p and p'
%
%% Initialization
matA=dataset;
[m,n]=size(matA);
diffmat=matA-repmat(p,1,n);
eudis=sum(diffmat.^2,1);
if ~isempty(varargin)
    if nargin>4
        return;
    end
    alpha_0= deal(varargin{:});
    [ral cal]=size(alpha_0);
    if ral>cal
        alpha_0=alpha_0';
    end
    alpha=alpha_0;
    p_p=matA*alpha';
else
    tmparr=find(eudis==min(eudis));
    min_index=tmparr(1);
    p_p=matA(:,min_index);
    alpha=zeros(1,n);
    alpha(min_index)=1;
end
%% Situation when input matrix has only one point
distance=min(eudis);
if(n<=1)
    if min(eudis)~=0
        inorout=0;
        p_prime=p_p;
        alpha_coe=alpha;
        dist=min(eudis);
        return;
    else
        inorout=1;
        p_prime=p;
        alpha_coe=alpha;
        dist=min(eudis);
        return;
        
    end
end

%% Iterative step
inorout=1;
iter=0;
while(sqrt((p-p_p)'*(p-p_p))>epsilon)
    iter=iter+1;
    distance=sqrt((p-p_p)'*(p-p_p));
    gd=matA'*(p-p_p);
    p_norm=p'*p;
    p_p_norm=p_p'*p_p;
    dist_diff=(p_norm-p_p_norm)- 2*gd;
    index_pivot=find(dist_diff<=0);
    if length(index_pivot)==0
        found=0;%% No pivot found
    else
        v_index=find(gd==max(gd));
        beta=(p_p-matA(:,v_index))'*(p_p-p)/((p_p-matA(:,v_index))'*(p_p-matA(:,v_index))); %% compute update step size
        alpha=(1-beta)*alpha;
        alpha(v_index)=alpha(v_index)+beta;
        p_p=(1-beta)*p_p+beta*matA(:,v_index);
        found=1;
    end
    if(found==0)
        inorout=0;
        break;
    end
end
p_prime=p_p;
alpha_coe=alpha;
dist=distance;