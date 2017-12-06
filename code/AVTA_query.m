function [index query_result p_p alpha]=AVTA_query(data_mat,query,eps_err)
%%******This program returns answer of convex hull query*********
%Input:
%   data_mat: data_num x dimension matrix. Rows as points.
%   query: dimension x 1 vector. The query point.
%   eps_err: precision parameter
%Output:
%   index: index of the vertices detected by the algorithm?
%   query_result: return 1 if query point is in the convex hull. 0
%   otherwise?
%   p_p: If query point is in the convex hull, give an epsilon approximate
%   otherwise p_p is a witness.
%   alpha: The weightes of vertices to represent p_p.
%%******************************************

%% Initialize the first vertex by the point with max l2-norm
Qmat=data_mat;
[NN MM]=size(Qmat);
index_set=[];
now_eps=max(2,eps_err+1);
eudis=sum(data_mat.*data_mat,2);
tmparr=find(eudis==max(eudis));
max_index=tmparr(1);


index_set=logical(zeros(1,NN));
vtx_ind_1=max_index;

index_set(vtx_ind_1)=1;

candidate_set=logical(ones(1,NN)); %% Indicator vector: 1 if i-th point is not in the current gamma distance of convex hull.
candidate_set(index_set)=0;
index_series=1:NN;
alpha_0=1;
pre_index=0;
p=query;
while now_eps>eps_err
    
    now_index_set=index_series(index_set);%% current vertices set
    this_data=data_mat(now_index_set,:)'; 
    
    %sum(index_set)
    %% when new vertex us added into current vertices set, update the size of weight matrix(representaion of every points using current vertices)
    if length(index_series(index_set))~=length(alpha_0) && sum(pre_index)~=0
        pre_series=index_series(pre_index);
        [tmp_srt,tmp_idx]=sort([pre_series,setdiff(now_index_set,pre_series)]);
        tmp_al=alpha_0;
        alpha_0=zeros(length( now_index_set),1);
        num_new=length(setdiff(now_index_set,pre_series));
        alpha_0(1:(end-num_new))=tmp_al;
        alpha_0=alpha_0(tmp_idx);
    end
    %% Triangle algorithm using previous iterate as warm up initialization.
    [inorout,p_prime,alpha_coe,dist]=TA_anti_warm(this_data,p,now_eps,alpha_0);
    %[inorout,p_prime,alpha_coe,dist]=anti_ta_warm(this_data,p,now_eps,alpha_0);
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%
    %disp('finish check if p is inside')
    
    %%%%%%%%%%%%%%%%%%%%%%
    
    %% compute new vertices when a separartion hyperplane is found
    alpha_0=alpha_coe';
    %alpha_coe
    if inorout==1
        now_eps=now_eps*0.9;
        alpha=alpha_coe;
        query_result=1;
    else
        direction=p_prime-p;
        S_index=index_series(candidate_set);
        S_data=Qmat(S_index,:);

        projected_val=S_data*direction;
        if min(projected_val)>p'*direction
            query_result=0;
            alpha=0;
            break;
        end
        %index_1=max(find(projected_val==max(projected_val)));%% The point with maxed and min
        index_2=min(find(projected_val==min(projected_val)));%% projected value must be vertices
        %index_candidate_1=S_index(index_1);
        index_candidate_2=S_index(index_2);

        pre_index=index_set;
        if index_set(index_candidate_2)==0
            %%%%%%%%%%%%%%%%%%%%%%
            %disp('found new vertex')

            %%%%%%%%%%%%%%%%%%%%%%
            
            
            index_set(index_candidate_2)=1;
            
            candidate_set(index_candidate_2)=0;
        end

    end

end

index=index_series(index_set);
p_p=p_prime;
