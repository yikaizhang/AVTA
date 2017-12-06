function [index,ta_freq]=AVTA_MultipleRp(data_mat,K,proj_num,proj_dim)
%%************Using multiple random projection to compute vertices*********
% Inputs: 
% A: sample num x dim matrix. Rows as sample points.
% K : Number of topics
% proj_num: The number of projections.
% proj_dim: The dimension of random projection.
% Output:
% ta_freq : The frequency vector. i-th entry is the frequency of i-th point
% detected as a vertex.
% 
%
%%****************************




matA_p=data_mat';
[m,n]=size(matA_p);
tmp_ta_freq=zeros(proj_num,n);


parfor proj=1:proj_num
    disp(['Now computing # '  num2str(proj) ' projection.']);
    %% Random projection
    dim_project=proj_dim;
    rand_matrix=random('Normal',0,1,m,dim_project);
    projected_Q=matA_p'*rand_matrix;
    projected_matA=projected_Q';

    recover_index_ta=AVTA_K(projected_matA',K);
    tmp_logic=zeros(1,n);
    tmp_logic(recover_index_ta)=1;
    tmp_ta_freq(proj,:)=tmp_logic;
end
ta_freq=sum(tmp_ta_freq,1);
[or_val_ta,or_ind_ta]=sort(ta_freq);
index=sort(or_ind_ta((end-K+1):end));
% 
% 
% 
% %% Cleaning step
% clean_gamma=now_gamma/10;
% clean_ed=0;
% while clean_ed<2*K
%     clean_ed=clean_ed+1;
%     [or_val_ta or_ind_ta]=sort(ta_freq,'descend');
%     ta_freq=or_ind_ta(1:2*K);
%     working_index=or_ind_ta(1:3*K);
%     working_mat=matA_p(:,working_index);
%     now_index=ta_freq(clean_ed);
%     working_vec=matA_p(:,now_index);
%     dist_mat=working_mat-working_vec(:,ones(1,3*K));
%     dist_val=sqrt(sum(dist_mat.*dist_mat,1));
%     index_2cl=find(dist_val<clean_gamma);
%     index_2cl_set=working_index(index_  2cl);
%     tmp_sum=sum(ta_freq(index_2cl_set));
%     ta_freq(index_2cl_set)=0;
%     ta_freq(now_index)=tmp_sum;
% end
% 
% [or_val_ta,or_ind_ta]=sort(ta_freq);
% ta_freq=sort(or_ind_ta((end-K+1):end));
% 
% 



