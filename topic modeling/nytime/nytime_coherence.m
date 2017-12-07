%%*******************************************************
%This program runs the nytime corpus experiment. It computes
%coherence of topics of algorithms learned from 'NYtime' dataset.
%This dataset has pruned vocalburary.
%'nytime_subset.mat' is the word- doc occurance matrix with pruned vocalb.
%It has 30000 documents which is uniformly sampled from whole dataset.
%'nytime_subset.mat' should be included in current path.
%Experiment setting:
%   training set size:20000.
%   testing set size: 10000
%%******************************************************

%% Experiment setting
load('nytime_subset.mat');
[word_num,doc_num]=size(occur);
tr_n=20000; %% training set size
Rep=10;


%% compute coherence
K=50;%number of topic
te_coe_tvt=zeros(1,K);
te_coe_tsvd=zeros(1,K);
te_coe_aa=zeros(1,K);
te_coe_at=zeros(1,K);
time_aa=0;
time_at=0;
time_tvt=0;
time_tsvd=0;
eps_val=0.01;
TAVTK=0;
%% 
for re=1:Rep
    tr_index=randperm(doc_num,tr_n);
    te_index=setdiff(1:doc_num,tr_index);
    tr_data=occur(:,tr_index);
    te_data=occur(:,te_index);
    this_occur=occur(:,te_index);

    [m,doc_n]=size(this_occur);
    occur_mat=zeros(m,doc_n);
    occur_mat(this_occur>0)=1;
    sum_occ_word=sum(occur_mat,2);
      
     
    tmp_start=tic;
    [tsvd cluster_id2]=TSVD2(tr_data,pwd,K);
    time_tsvd=time_tsvd+toc(tmp_start);
    
    tmp_start=tic;
    [aa_A]=recover_l2(tr_data,K);
    time_aa=time_aa+toc(tmp_start);
    
    tmp_start=tic;
    [at_A]=recover_l2_avta(tr_data,K);
    time_at=time_at+toc(tmp_start);
    
    tmp_start=tic;
    [tvt cluster_id2]=AVTA_Catch(tr_data,K);
    time_tvt=time_tvt+toc(tmp_start);


    for ii=1:K
        [stdv,odv]=sort(tsvd(:,ii),'descend');
        index_word=odv(1:5);

        for jj=1:length(index_word)
            for kk=1:length(index_word)
                te_coe_tsvd(ii)=te_coe_tsvd(ii)+log((sum(occur_mat(index_word(jj),:).*occur_mat(index_word(kk),:))+eps_val)/(sum_occ_word(index_word(kk))));
            end
        end
    end

    for ii=1:K
        [stdv,odv]=sort(aa_A(:,ii),'descend');
        index_word_aa=odv(1:5);

        for jj=1:length(index_word_aa)
            for kk=1:length(index_word_aa)
                te_coe_aa(ii)=te_coe_aa(ii)+log((sum(occur_mat(index_word_aa(jj),:).*occur_mat(index_word_aa(kk),:))+eps_val)/(sum_occ_word(index_word_aa(kk))));
            end
        end
    end

    for ii=1:K
        [stdv,odv]=sort(at_A(:,ii),'descend');
        index_word_at=odv(1:5);

        for jj=1:length(index_word_at)
            for kk=1:length(index_word_at)
                te_coe_at(ii)=te_coe_at(ii)+log((sum(occur_mat(index_word_at(jj),:).*occur_mat(index_word_at(kk),:))+eps_val)/(sum_occ_word(index_word_at(kk))));
            end
        end
    end

    for ii=1:K
        [stdv,odv]=sort(tvt(:,ii),'descend');
        index_word=odv(1:5);

        for jj=1:length(index_word)
            for kk=1:length(index_word)
                te_coe_tvt(ii)=te_coe_tvt(ii)+log((sum(occur_mat(index_word(jj),:).*occur_mat(index_word(kk),:))+eps_val)/(sum_occ_word(index_word(kk))));
            end
        end
    end
    disp(['Fast anchor+recoverl2 coherence:' num2str(mean(te_coe_aa./re))])
    disp(['AVTA+recoverl2 coherence:' num2str(mean(te_coe_at./re))])
    disp(['Tsvd coherence:' num2str(mean(te_coe_tsvd./re))])
    disp(['AVTA+Catchword coherence :' num2str(mean(te_coe_tvt./re))])
    
    disp(['Fast anchor+recoverl2 coherence:' num2str(mean(te_coe_aa./re))])
    disp(['AVTA+recoverl2 coherence:' num2str(mean(te_coe_at./re))])
    disp(['Tsvd coherence:' num2str(mean(te_coe_tsvd./re))])
    disp(['AVTA+Catchword coherence :' num2str(mean(te_coe_tvt./re))])

    disp(['Fast anchor+recoverl2 time:' num2str(time_aa/re)])
    disp(['AVTA+recoverl2 time:' num2str(time_at/re)])
    disp(['Tsvd time:' num2str(time_tsvd/re)])
    disp(['AVTA+Catchword time :' num2str(time_tvt/re)])
end


