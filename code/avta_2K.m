function [index]=avta_K(mat_A,proj_dim,K)

%%*********Compute K vertices using avta_K
%Input: 
%   mat_A: sample number x dimension matrix. Rows are points.
%   proj_dim: The dimension of random projection.
%   K: Number of vertices.
%Output:
%   index: The index of vertices.

[NN MM]=size(mat_A);

%% Random projection
dim_project=proj_dim;
rand_matrix=random('Normal',0,1,MM,dim_project);

projected_Q=mat_A*rand_matrix;
mat_A=projected_Q;
Qmat=mat_A;

%% Compute vertices
index_set=[];
epsilon=0.5;
eudis=sum(mat_A.*mat_A,2);
tmparr=find(eudis==max(eudis));
max_index=tmparr(1);
index_set=logical(zeros(1,NN));
vtx_ind_1=max_index;
index_set(vtx_ind_1)=1;
sum(index_set)

S_set=logical(ones(1,NN));
candidate_set=S_set;
candidate_set(index_set)=0;
index_series=1:NN;
alpha_zero=ones(sum(index_set),NN);

pre_index=index_set;

%% Stop if K vertices found
while sum(index_set)<K
    S_set=logical(ones(1,NN));
    candidate_set=S_set;
    candidate_set(index_set)=0;
    
    
    while sum(candidate_set)>0
        [kk NN]=size(alpha_zero);
        if length(index_series(index_set))~=kk && sum(pre_index)~=0
            pre_series=index_series(pre_index);
            now_index_set=index_series(index_set);
            [tmp_srt,tmp_idx]=sort([pre_series,setdiff(now_index_set,pre_series)]);
            
            tmp_al=alpha_zero;
            alpha_zero=zeros(length( now_index_set),NN);
            num_new=length(setdiff(now_index_set,pre_series));
            alpha_zero(1:(end-num_new),:)=tmp_al;
            alpha_zero=alpha_zero(tmp_idx,:);
            
            
            pre_index=index_set;
            
        end
        
        
      
        
        
        current_index=index_series(candidate_set);
        current_size=length(current_index);
        rnd_ind=unidrnd(current_size,1,1); 
        this_index=current_index(rnd_ind);
        this_data=mat_A(index_set,:)';
        p=mat_A(this_index,:)';
        [inorout,p_prime,alpha_coe,dist]=anti_ta_warm(this_data,p,epsilon,alpha_zero(:,this_index));
        if inorout==1
            candidate_set(this_index)=0;
            alpha_zero(:,this_index)=alpha_coe';
        else



            direction=p_prime-p;
            S_index=index_series(candidate_set);
            S_data=Qmat(S_index,:);

            projected_val=S_data*direction;
            index_2=min(find(projected_val==min(projected_val)));
            index_candidate_2=S_index(index_2);

            
            if index_set(index_candidate_2)==0
                index_set(index_candidate_2)=1;
                if sum(index_set)>=K
                    break;
                    
                end
                candidate_set(index_candidate_2)=0;
            end

        end
    end
    epsilon=epsilon*0.7;
end

index=index_series(index_set);

%

