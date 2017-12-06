%%%%input dataset%%
%%%output vertices found%%
function [k_index,all_alpha,dist,rand_Q]=AVTA_prune(Q_bar,K,proj_dim,gamma)
%%*************Compute K vertices of a given set of points*****************
% Input: 
%       Q_bar: # of pts x dim matrix. Columns as points.
%       K: # of vertices for output
%       proj_dim: # of dimension for random projection.
%       gamma: precision parameter.
%Output:
%       k_index: Index of vertices
%       all_alpha: Weights of representing points using vertices.
%       dist: K x1 vector. Distance of every vertex to convex hull of other
%       vertices.
%       rand_Q: The projected data matrix.
%%*********************************************************



[NN MM]=size(Q_bar);


%% Random projection
dim_project=proj_dim;
rand_matrix=random('Normal',0,1,MM,dim_project);
projected_Q=Q_bar*rand_matrix;
Q_bar=projected_Q;
rand_Q=projected_Q;
Q_bar=Q_bar;

index_set=[];


epsilon=2;

eudis=sum(Q_bar.*Q_bar,2);
tmparr=find(eudis==max(eudis));
max_index=tmparr(1);


index_set=logical(zeros(1,NN));
vtx_ind_1=max_index;

index_set(vtx_ind_1)=1;

S_set=logical(ones(1,NN));
candidate_set=S_set;
candidate_set(index_set)=0;
index_series=1:NN;
alpha_zero=ones(sum(index_set),NN);

pre_index=index_set;

%% compute a super set of vertices
while epsilon>gamma || sum(index_set)<K
    S_set=logical(ones(1,NN));
    candidate_set=S_set;
    candidate_set(index_set)=0;
    
    if sum(index_set)>=K
        break;
    end
%     disp(['number of vtx of convex hull is ' num2str(sum(index_set))]);
    %% Check if every point is in distance of 'epsilon' to current convex hull
    while sum(candidate_set)>0
%         disp(['number of pts in the convex hull is ' num2str(sum(candidate_set))]);
        current_index=index_series(candidate_set);
        current_size=length(current_index);
        rnd_ind=unidrnd(current_size,1,1); 
        this_index=current_index(rnd_ind);
        this_data=Q_bar(index_set,:)';
        p=Q_bar(this_index,:)';
        [inorout,p_prime,alpha_coe,dist]=anti_ta_warm(this_data,p,epsilon,alpha_zero(:,this_index));
        if inorout==1 
            candidate_set(this_index)=0;
            alpha_zero(:,this_index)=alpha_coe';
        else
            %% Find the farthest point rpt the witness direction.
            direction=p_prime-p;
            S_index=index_series(candidate_set);
            S_data=Q_bar(S_index,:);

            projected_val=S_data*direction;
            index_2=min(find(projected_val==min(projected_val)));
            index_candidate_2=S_index(index_2);
            candidate_set(index_candidate_2)=0;
            
            
            %% Adding new vertices.
            if index_set(index_candidate_2)==0
                index_set(index_candidate_2)=1;
                
                
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
                
            end

        end
    end
    %epsilon
    epsilon=epsilon*0.9;
end


super_index=index_series(index_set);
%disp(['number of vertices before prune is ' num2str(length(super_index))]);
superset_data=Q_bar(super_index,:)';


%% pruning 

if length(super_index)>K
    index_set=vtx_prune(superset_data,K,'forward');
    k_index=index_series(index_set);
    %disp('finish prune');
else
    k_index=super_index;
end



%% compute vertex distance 
all_alpha=zeros(K,size(Q_bar,1));
dis_vec=zeros(1,K);

for i=1:K
    p=Q_bar( k_index(i),:)';
    tmp_ind=setdiff(k_index,k_index(i));
    current_cols=Q_bar(tmp_ind,:)';
    [alp dis]=frank_wolfe_2(p, current_cols, 0.01);
    dis_vec(i)=dis;
end
dist=mean(dis_vec);

    




