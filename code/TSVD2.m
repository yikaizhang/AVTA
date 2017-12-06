function [M_hat,cluster_id2] = TSVD2(A,out,K,varargin)
%************ Thresolded SVD based K-Means for Topic Recovery *************
% 
% [M_hat, dominantTopic] = TSVD2s(A, outpath, K)
%
% Inputs:
% A: Word-Doc matrix.
% outpath : Output folder (final topic matrix will be written here, is not used. Can be used if uncommented)
% K : Number of topics
%
% Output:
% M_hat : Final estimated topic matrix
% dominantTopic : Cluster (dominant topic) assignment of documents
%
% Optional Inputs (Algorithm Parameters, need to be passed in order):
% (For details refer to paper, steps here refer to TSVD algo section 3.3)
% w_0 : double, lower bound of pobability of dominant topic, see step 1.a
%       (default 1/k)
% \epsilon_1 : double, \epsilon parameter for second thresolding in 
%              step 1(a) (default 1/6)
% \epsilon_2 : double, \epsilon_0 for computing g(.,.) in step 5(a)
%              (default 1/6)
% \gamma : double, for finding the set J(.) in step 5.b (default 1.1)
% \epsilon_3 : double, \epsilon_0 for finding top documents, step 6 
%              (default 1.0)
% reps : int, number of repitions for the first k-means step 4.a, more 
%        repititions are better though computationally expensive 
%        (default depends on data size: 1,3,5)
% 
% 
% Example:
% outpath = 'out/';
% K = 50;
% [M, labels] = TSVD2s(A, outpath, K);
%

%
%-------------------------------------------------------------------------
% Original Author: Trapit Bansal.
% 
%-------------------------------------------------------------------------
% 


if ~isempty(varargin)   % Algo parameters
    [w0,eps1,eps2,rho,eps3,reps] = deal(varargin{:});
else
        w0 = 1.0/K;
        eps1 = 1/6;
        eps2 = 1/3;
        rho = 1.1;
        eps3 = 1.0;
        reps = 3;  % will be changed later is data is too big or too small
end

% A = Amatrix(infile);  %Reads the data and returns the A matrix
[w,d] = size(A);
fprintf('\nNo of documents read:%d\n',d);

if isempty(varargin)
    if d <= 5000, reps = 5;
    elseif d >= 50000, reps = 3;
    else reps = 5;
    end
end
reps = 10;
% 
% if w*d > 1e8
%     save(strcat(outpath, '/A.mat'),'A');  % Save A for later use, if its large
% end

start_time = tic;

% A = A*spdiags(1./sum(A,1)',0,d,d); %normalized A

fprintf('Thresholding on A\n')

A = A'; % Row slicing is inefficient for sparse matrix
[B, thres] = threshold_sparse(A, w0, eps1, w, d);
% [B, thres] = threshold_sparse(A, w0, eps1, w, d);


% if w*d > 1e8, clear A; % clear to save space
% else A = A'; end
A = A'; 

% Uncomment the follwing line to write the thresholds (\zeta)
% dlmwrite(strcat(outpath,'/Thresholds'),full(thres),'delimiter','\n');
clear thres

retained_docs = find(sum(B,1)~=0);   %Columns which are not completely zero
B = B(:,retained_docs);

%% Computing the K-rank approximation for the matrix B
fprintf('Computing SVD Projection\n')

if w*d <= 5e7 && w > d
    [~,S,V] = svds(B,K);
    B_k = sparse(S)*V';
    clear S V;
else
    % computing BB^T then finding top eigenvectors!
    BBt = B*B';
%     BBt = spSymProd(B);
    [U,~] = eigs(BBt,[],K,'lm',struct('issym',1));
    clear BBt;
    B_k = U'*B;
    clear U;
end

% 
%% K-means on projected matrix
fprintf('Performing k-means on columns of B_k\n')

q_best = Inf;
cluster_id = [];
for r = 1:reps
    [~,ini_center] = kmeanspp_ini(B_k,K);
    [~, c_id, ~,q2,~] = kmeans_fast(B_k,ini_center,2,reps > 1);
    if q2 < q_best
        cluster_id = c_id; 
        q_best = q2;
    end
end
% fprintf('Quality:%.4f\n',q_best);

clear B_k ini_center c_id;

P1 = zeros(K,w);  % Finding centers in original space
for k=1:K
    cols = find(cluster_id==k);
    P1(k,:) = sum(B(:,cols),2)./length(cols);
end

% Uncomment following line to write clusting info
% doc_ids = zeros(d,1) - 1;
% doc_ids(retained_docs) = cluster_id; %removed docs are -1
% dlmwrite(strcat(outpath,'/P1'),full(P1),'delimiter',' ','precision','%.6f');
% dlmwrite(strcat(outpath,'/doc_cluster_id'),full(doc_ids),'\n');
% dlmwrite(strcat(outpath,'/clusterID'),doc_ids,'\n');
% clear doc_ids

%% Lloyds on B with start from B_k
fprintf('Performing Lloyds on B with centers from B_k clustering\n')

[~, cluster_id2, ~,q2,~] = kmeans_fast(B,P1,2,0);
clear B cluster_id
% fprintf('Quality:%.4f\n',q2);

% 
% if w*d > 1e8, load(strcat(outpath, '/A.mat')); end
A = A*spdiags(1./sum(A,1)',0,d,d); %normalized A
A1_rowsum = full(sum(A,2));

P2 = zeros(w,K);    % this will be topic matrix without using cathwords
for k=1:K
    cols = retained_docs(cluster_id2==k);
    P2(:,k) = sum(A(:,cols),2)./length(cols);
end

% Uncomment following two lines to write clustering info
% dlmwrite(strcat(outpath,'/P2'),full(P2'),'delimiter',' ','precision','%.6f');
% dlmwrite(strcat(outpath,'/clusterID2'),cluster_id2,'\n');

%% Find Catchwords
fprintf('Finding Catchwords\n');

fractiles = zeros(w,K); % This will store the values g(i,l)

for l=1:K
    T = A(:,retained_docs(cluster_id2==l))'; %columns of T are words
    % sorting on columns is faster for sparse matrix
    T = sort(T,1,'descend'); % sort cols in descending
    fractiles(:,l) = T(min(floor(eps2*w0*d/2),size(T,1)),:);
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
for l=1:K
    % check that frequency over catchwords is not very small
    if (~ismember(l,catchless) && sum(A1_rowsum(catchword(:,l))) <= 0.01*d/(2*K))
        catchless = horzcat(catchless,l);
    end
end

if ~isempty(catchless)
    fprintf('Catchless topics: ');
    fprintf('%d ',catchless); fprintf('\n');
end

% Uncomment the following line to write the catchords indicator matrix
% dlmwrite(strcat(outpath,'/catchwords'),full(catchword),'delimiter',' ');

M_hat = zeros(w,K);
for l=1:K
    if ismember(l,catchless)
        M_hat(:,l) = P2(:,l);
        continue;
    end
    n = max(floor(eps3*w0*d/2),1); % new - 08/14
    [~,inds1]=sort(sum(A(catchword(:,l),:),1),'descend');
    alpha1 = inds1(1:n);
    M_hat(:,l) = sum(A(:,alpha1),2)*1.0/n;
end



end_time = toc(start_time);
% 
% fprintf('Writing topic matrix\n');
dlmwrite(strcat(out,'/M_hat'),full(M_hat'),'delimiter',' ','precision','%.6f');

fprintf('\nAll Done, algorithm took %.2f seconds\n', end_time);
end
