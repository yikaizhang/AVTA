function [k_index,all_alpha]=AVTA_rep(Q_bar,K,proj_dim,gamma,rep_t)
%%*****************A practical implementation of AVTA*****************
% Repeat AVTA for multiple times and pick the best one for further
% application. 
% Input: 
%       Q_bar? # of pts x dim matrix. Columns as points.
%       K? # of vertices for output
%       proj_dim: # of dimension for random projection.
%       gamma: precision parameter.
%       rep_t: Number of repetance.
%Output:
%       k_index: Index of vertices
%       all_alpha: Weights of representing points using vertices.
%%**********************************

disp('start computing vertices')
err_vec=zeros(1,rep_t);
k_ind_all=zeros(rep_t,K);
rand_mat=zeros(proj_dim,size(Q_bar,1),rep_t);
for i=1:rep_t
    [k_index,all_alpha_i,err,rnd_mat]=AVTA_prune(Q_bar,K,proj_dim,gamma);
    err_vec(i)=err;
    k_ind_all(i,:)=k_index;
    rand_mat(:,:,i)=rnd_mat';
    %freq(k_index)=freq(k_index)+1;
end

disp('finish compute vertices')
min_err_ind=find(err_vec==max(err_vec));

k_index=k_ind_all(min_err_ind,:);


disp('start computing alpha');
data_Q=rand_mat(:,:,min_err_ind);
all_alpha=zeros(K,size(Q_bar,1));
vtx_data=data_Q(:,k_index);
for i=1:size(Q_bar,1)
    p=data_Q(:,i);
    [alp dis]=frank_wolfe_2(p, vtx_data, 0.01);
    all_alpha(:,i)=alp;
    %i
end
disp('finish computing alpha');
    