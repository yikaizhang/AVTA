
function [at_A]=recover_l2_avta(occur,K)

%%*********recoverl2 using avta to find vertices*********
%Input:
%   occur: w-d word document occurance matrix
%   K: number of topics
%Output:
%   at_A: word topic distribution
%   
%%***************************************************************

occur=full(occur);
[num_word,num_doc]=size(occur);

nd=sum(occur,1); %number of documents



%% compute the normalized word-word co-occurance matrix

hh=zeros(num_word,num_doc);
avg_hd=zeros(num_word,1);
for i=1:num_doc
%      i
    hd=occur(:,i);
    hd_t=hd./sqrt(nd(i)*(nd(i)-1));
    hd_hat=(hd)/(nd(i)*(nd(i)-1));
    hh(:,i)=hd_t;
    
    avg_hd=avg_hd+hd_hat;
end
QQ=hh*hh'-diag(avg_hd);

%% Normalize word-word co-occurance matrix

Q_bar=zeros(num_word,num_word);
row_index=logical(ones(1,num_word));
sum_row=sum(QQ,2);

for i=1:num_word
%     i
    sum_row(i)=sum(QQ(i,:));
    Q_bar(i,:)=QQ(i,:)./(sum_row(i)+eps);
%     if(sum(isnan(Q_bar(i,:)))>0)
%         row_index(i)=0;
%         pause();%% if a word never appears, stop and redo the data preprocessing step
%     end
    
end 

nips_Q_bar=Q_bar;


%% compute vertices using avta
at_index=AVTA_K(nips_Q_bar,K,'rp',K);


%% compute weight by minimizing l2-distance
data_at=nips_Q_bar(at_index,:);
coe_mat_at=zeros(size(data_at,2),K);
for i=1:size(data_at,2)
    [alpha_at,err_at]=frank_wolfe_2(nips_Q_bar(i,:)',data_at',0.0001);
    coe_mat_at(i,:)=alpha_at';
end


%% Recover word-topic matrix
at_Ap=diag(sum_row)*coe_mat_at;
at_sum_cA=sum(at_Ap,1);
at_A=at_Ap./(at_sum_cA(ones(size(data_at,2),1),:));






