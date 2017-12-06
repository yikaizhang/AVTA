%%%%input dataset%%
%%%output vertices found%%
function [index]=AVTA_eps(mat_p,GAMMA,varargin)

%%*********Compute K vertices using avta_K
%Input: 
%   mat_p: sample number x dimension matrix. Rows are points.
%   GAMMA: The precision parameter
%   varargin: ('Rp', dim of projection),  if using random projection, 
%Output:
%   index: index of vertices
%%***************************************
nargin
if ~isempty(varargin)
    if nargin>=4
        [Type Para]=deal(varargin{:});
    else
        Type=deal(varargin{:});
        Para=size(mat_p,2)*2;
    end
end

if nargin>=4 && strcmp(Type,'Rp')==1
    [NN MM]=size(mat_p);
    rand_matrix=random('Normal',0,1,MM,Para);
    mat_p=mat_p*rand_matrix;
end



[NN,MM]=size(mat_p);

index_set=[];


epsilon=2;

eudis=sum(mat_p.*mat_p,2);
tmparr=find(eudis==max(eudis));
max_index=tmparr(1);


index_set=logical(zeros(1,NN));
vtx_ind_1=max_index;

index_set(vtx_ind_1)=1;
%sum(index_set)

S_set=logical(ones(1,NN));
%S_set(index_set)=1;
candidate_set=S_set;
candidate_set(index_set)=0;
index_series=1:NN;
alpha_zero=ones(sum(index_set),NN);

pre_index=index_set;

while epsilon>GAMMA
    S_set=logical(ones(1,NN));
    %S_set(index_set)=1;
    candidate_set=S_set;
    candidate_set(index_set)=0;
    
    
    
    while sum(candidate_set)>0
        [kk NN]=size(alpha_zero);
        if length(index_series(index_set))~=kk && sum(pre_index)~=0
            %pre_index
            %sum(index_set)
            pre_series=index_series(pre_index);
            now_index_set=index_series(index_set);
            [tmp_srt,tmp_idx]=sort([pre_series,setdiff(now_index_set,pre_series)]);
            
            tmp_al=alpha_zero;
            alpha_zero=zeros(length( now_index_set),NN);
            num_new=length(setdiff(now_index_set,pre_series));
            
%             size(tmp_al)
%             size(alpha_zero)
%             length(now_index_set)
%             num_new
            alpha_zero(1:(end-num_new),:)=tmp_al;
            alpha_zero=alpha_zero(tmp_idx,:);
            
            
            pre_index=index_set;
            
        end
        
        
      
        
        
%          disp('# of pts to be tested :')
%          sum(candidate_set)%% Uncomment if output # of vertices && # of
%          disp('# of vertices is :')%%pts remains
%          length_ind=sum(index_set)%%
        current_index=index_series(candidate_set);
        current_size=length(current_index);
        rnd_ind=unidrnd(current_size,1,1); 
        this_index=current_index(rnd_ind);
        this_data=mat_p(index_set,:)';
        p=mat_p(this_index,:)';
        [inorout,p_prime,alpha_coe,dist]=TA_anti_warm(this_data,p,epsilon,alpha_zero(:,this_index));
        %[inorout,p_prime,alpha_coe,dist]=anti_ta_warm(this_data,p,epsilon,alpha_zero(:,this_index));
        if inorout==1
            %S_set(this_index)=0;
            candidate_set(this_index)=0;
%             22
%             size(alpha_zero(:,this_index))
%             
%             size(this_data)
%             size(alpha_coe')
            alpha_zero(:,this_index)=alpha_coe';
        else



            direction=p_prime-p;
            S_index=index_series(candidate_set);
            S_data=mat_p(S_index,:);

            projected_val=S_data*direction;
            index_2=min(find(projected_val==min(projected_val)));
            index_candidate_2=S_index(index_2);

            
            if index_set(index_candidate_2)==0
                index_set(index_candidate_2)=1;

                candidate_set(index_candidate_2)=0;
            end

        end
    end
    %%disp('current epsilon is : ')
    %%epsilon
    epsilon=epsilon/2;
end
index=index_series(index_set);

