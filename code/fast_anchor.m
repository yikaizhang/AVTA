function [index]=fast_anchor(Data_m,K,randpro)
%%*****Using fast anchor to compute vertices*************
%Input:
%   Data_m: sample number x dimension matrix. Rows of 'Data_m' as a point.
%   K: Number of vertices.
%   randpro: Using randomg projection if set to 1, 0 otherwise,.
%Output:
%   index: The index of vertices.
%%******************************************
[samplen dimen ]=size(Data_m);
if randpro==1
    dim_project=500;
    rand_matrix=random('Normal',0,1,dimen,dim_project);
    projected_Q=Data_m*rand_matrix;
    data_m=projected_Q;
else
    data_m=Data_m;
end



[mm,nn]=size(data_m);
pqnorm = sum(data_m.^2,2);
S_index=logical(zeros(1,mm));
Candi_index=logical(ones(1,mm));
max_norm=find(pqnorm==max(pqnorm));
S_index(max_norm)=1;
ind_series=1:mm;



for i=2:K
    S_mat= data_m(S_index,:)';
    [Qv,Rt]=qr(S_mat);
    [rtn,rtm]=size(Rt);
    if rtm>size(Qv,2)
        break;
    end
    ppq=data_m*Qv(:,1:rtm);
    ppqnorm=sum(ppq.^2,2);
    dis=pqnorm-ppqnorm;
    max_ind=min(find(dis==max(dis)));
    S_index(max_ind)=1;

end
for i=1:(K-1)
    temp_index=S_index;
    ind_set_tmp=ind_series(S_index);
    if i>length(ind_set_tmp)
        break;
    end
    temp_index(ind_set_tmp(i))=0;
    S_mat= data_m(temp_index,:)';
    [Qv,Rt]=qr(S_mat);
    [rtn rtm]=size(Rt);
    ppq=data_m*Qv(:,1:rtm);
    ppqnorm=sum(ppq.^2,2);
    dis=pqnorm-ppqnorm;
    max_ind=min(find(dis==max(dis)));
    temp_index(max_ind)=1;
    S_index=temp_index;

end
index=ind_series(S_index);