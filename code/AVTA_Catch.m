function [M_hat,cluster_id2] = AVTA_Catch(A,K,varargin)
%%************ AVTA + Catchword*************
% 
% [M_hat, dominantTopic] = AVTA_Catch(A, K)
%
% Inputs: 
% A:Word-Doc matrix 
% K : Number of topics
%
% Output:
% M_hat : Final estimated topic matrix
% dominantTopic : Cluster (dominant topic) assignment of documents
%
% 
%%****************************************

[w,d] = size(A);
fprintf('\nNo of documents read:%d\n',d);

if ~isempty(varargin)   % Algo parameters
    [w0,eps1,eps2,rho,eps3,reps] = deal(varargin{:});
else
        w0 = 1/K;
        eps1 = 1/6;
        eps2 = 1/3;
        rho =1.2;
        eps3 = 1.0;
        reps = 3;  % will be changed later is data is too big or too small
end


%% Threshold on the document
tic;
[B, thres] = threshold_sparse(A', w0, eps1, w, d);

retained_docs = find(sum(B,1)~=0);
A= A(:,retained_docs);
[w,d] = size(A);

%% Compute the KxD weight matrix
start_time = tic;
sum_a=sum(A,1);
mat_App = A*spdiags(1./sum_a',0,d,d); % normalize w-d matrix 
size(mat_App)
if d <5000
    reps = 10;
else
    reps = 3;
end
[k_ind,c_val]=AVTA_rep(mat_App',K,2*K,0.05,reps); %% Find vertices and compute weight matrix 
CC=c_val; 
st_k=sort(k_ind);

%% Thresholding and normalize
CC(c_val<0.2/K)=0;
sum_c=sum(CC,1);
retained_docs = find(sum_c~=0); %% remove all zero vectors
A= A(:,retained_docs);
[w,d] = size(A);
mat_A = A*spdiags(1./sum_a',0,d,d); %normalized A
sum_a=sum(A,1);


CC = CC(:,retained_docs)*diag(1./sum_c(retained_docs)); %normalized weight matrix 



%% Initialize partition by majoirty weight before lloyd algorithm

cluster_id=zeros(1,K);

for i=1:d
        cluster_id(i)=min(find(CC(:,i)==max(CC(:,i))));
end
cluster_id(st_k)=1:K;


P1 = zeros(K,length(k_ind));  % Finding centers in original space

for k=1:K
    cols = find(cluster_id==k);
    P1(k,:) = sum(CC(:,cols),2)'./length(cols);
end


%% Lloyds on weight matrix
fprintf('Performing Lloyds on weight matrix with centers major weight \n')

[~, cluster_id2, ~,q2,~] = kmeans_fast(CC,P1,2,0);
clear B cluster_id
toc;

A = mat_A;


%% For topics without catchword, output mean of the cluster

P2 = zeros(w,K);    % this will be topic matrix without using cathwords
for k=1:K
    cols = find(cluster_id2==k);%
    if  length(cols)>0.01*d
        [importance,ori]=sort(CC(k,cols));
        cols=cols(ori(1:floor(0.01*d)));
    end
    
    P2(:,k) = sum(A(:,cols),2)./length(cols);
%     length(cols)
    if length(cols)==0
        P2(:,k)=A(:,st_k(k));% 
    end
end







%% Find Catchwords
fprintf('Finding Catchwords\n');
fractiles = zeros(w,K); % This will store the values g(i,l)

for l=1:K
    if sum(cluster_id2==l)>=2
        T = mat_A(:,cluster_id2==l)'; %columns of T are words
        T = sort(T,1,'descend'); % sort cols in descending

        fractiles(:,l) = T(min(floor(eps2*w0*d/2),size(T,1)),:); % compute the top words
    end
end
clear T;

catchword = false(w,K);
for l =1:K
    for i=1:w
        catchword(i,l) = false;
        fractile_1 = fractiles(i,l);
        isanchor = false;
        for l2 = 1:K
            if (l2==l), continue; end
            fractile_2 = fractiles(i,l2);
            isanchor = (fractile_1 > rho*fractile_2);
            if ~isanchor
                break
            end
        end
        if isanchor
            catchword(i,l)  = true;
        end
    end
end

catchy_topics = find(sum(catchword,1)~=0);

catchless = setdiff(1:K,catchy_topics);
A1_rowsum=sum(mat_A,2);
for l=1:K
    % check that frequency over catchwords is not very small
    
     if (~ismember(l,catchless) && sum(A1_rowsum(catchword(:,l))) <=0.01*d/(3*K))
        catchless = horzcat(catchless,l);
    end
end

if ~isempty(catchless)
    fprintf('Catchless topics: ');
    fprintf('%d ',catchless); fprintf('\n');
end



M_hat = zeros(w,K);
for l=1:K
    if ismember(l,catchless)
        M_hat(:,l) =P2(:,l);
        continue;
    end
    n = max(floor(eps3*w0*d/2),1); 
    [~,inds1]=sort(sum(mat_A(catchword(:,l),:),1),'descend');
    alpha1 = inds1(1:n);
    M_hat(:,l) = sum(mat_A(:,alpha1),2)*1.0/n;
end

toc;

end_time = toc(start_time);

%fprintf('\nAll Done, algorithm took %.2f seconds\n', end_time);
end
