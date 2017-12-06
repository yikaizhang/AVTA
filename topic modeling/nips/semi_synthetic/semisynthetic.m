%%*******************************************************
%   This program runs the semi synthetic data experiment.
%%*******************************************************
%parpool(4); If parallel is used, uncomment second loop->parfor re=1:rep.
load('nips_2000_wd.mat');% load the word topic distribution
rep=8;% repeat t times
number_doc_list=(1:10)*5000; % Experiment runs with different number of documents from 5000-50000

%% Initialization
tsvd_err=zeros(rep,length(number_doc_list)); %error of tsvd
tvt_err=zeros(rep,length(number_doc_list));  %error of AVTA+Catch
aa_err=zeros(rep,length(number_doc_list));   %error of l2+Fast anchor
at_err=zeros(rep,length(number_doc_list));   %error of l2+AVTA

time_tsvd=zeros(1,length(number_doc_list));  %running time of tsvd
time_tvt=zeros(1,length(number_doc_list));   %running time of AVTA+Catch
time_aa=zeros(1,length(number_doc_list));    %running time of l2+Fast anchor
time_at=zeros(1,length(number_doc_list));    %running time of l2+AVTA
currentFolder = pwd;  %use current file to wirte results
savefile = [pwd datestr(now,30) '.mat'];
writefile=pwd;

for iter=1:length(number_doc_list)
	parfor re=1:rep
    %for re=1:rep
        disp(['iter=',num2str(iter)]);
        disp(['re=',num2str(re)]);
        K=50;
        sample_num=number_doc_list(iter);
        word_num=1000;
        disp('generate data');
        dataset=gibb_art(scal,50,word_num,sample_num);
        disp('end generate data');
        [art_num_word,art_num_doc]=size(dataset);
        tic
        [tvt cluster_id2]=AVTA_Catch(dataset,50);
        tvtt=toc;
        time_tvt(re,iter)=tvtt;
        tic
        [tsvd cluster_id2]=TSVD2(dataset,writefile,K);
        tsvdt=toc;
        time_tsvd(re,iter)=tsvdt;
        tic
        [aa_A]=recover_l2(dataset,K);
        aat=toc;
        time_aa(re,iter)=aat;
        tic
        [at_A]=recover_l2_avta(dataset,K);
        att=toc;
        time_at(re,iter)=att;
        error_tvt=zeros(1,K);
        for i=1:K 
            %i
            orin_vtx=scal(:,i);
            diffmat=tvt-orin_vtx(:,ones(K,1));
            distance_v=sum(abs(diffmat),1);
            error_tvt(i)=min(distance_v);
        end
        tvt_err(re,iter)=mean(error_tvt);
        error_tvd=zeros(1,K);
        for i=1:K 
            %i
            orin_vtx=scal(:,i);
            diffmat=tsvd-orin_vtx(:,ones(K,1));
            distance_v=sum(abs(diffmat),1);
            error_tvd(i)=min(distance_v);
        end         
       tsvd_err(re,iter)=mean(error_tvd);
       error_aa=zeros(1,K);
        for i=1:K 
            orin_vtx=scal(:,i);
            diffmat=aa_A-orin_vtx(:,ones(K,1));
            distance_v=sum(abs(diffmat),1);
            error_aa(i)=min(distance_v);
        end         
       aa_err(re,iter)=mean(error_aa);
       error_at=zeros(1,K);
        for i=1:K 
            orin_vtx=scal(:,i);
            diffmat=at_A-orin_vtx(:,ones(K,1));
            distance_v=sum(abs(diffmat),1);
            error_at(i)=min(distance_v);
        end         
       at_err(re,iter)=mean(error_at);

    end
end
save(savefile,'number_doc_list','tsvd_err','tvt_err','at_err','aa_err','time_tvt','time_tsvd','time_aa','time_at');



%% Plot the error
hold on;
xlabel('Number of documents') % x-axis label
ylabel('l1 error') % y-axis label
plot(number_doc_list,sum(tsvd_err,1)./rep,'--','Color',[1 0 1],'DisplayName','TSVD','LineWidth',1.5);
plot(number_doc_list,sum(tvt_err,1)./rep,'^-','Color',[0 0 1],'DisplayName','AVTA+Catchword','LineWidth',1.5);
plot(number_doc_list,sum(aa_err,1)./rep,'*-','Color',[0.5 0 1],'DisplayName','Fast Anchor+Recoverl2','LineWidth',1.5);
plot(number_doc_list,sum(at_err,1)./rep,'--','Color',[0.5 0.5 0.5],'DisplayName','AVTA+Recoverl2','LineWidth',1.5);
hold off;

legend('show','Location','east')%,'Orientation','horizontal')



%% Plot the running time
hold on;
xlabel('Number of documents') % x-axis label
ylabel('Running time') % y-axis label
plot(number_doc_list,sum(time_aa,1),':*','Color',[0.5 0 1],'DisplayName','Fast anchor+Recoverl2','LineWidth',1.5);
plot(number_doc_list,sum(time_tvt,1),':^','Color',[0 0 1],'DisplayName','AVTA+Catch Word','LineWidth',1.5);
plot(number_doc_list,sum(time_tsvd,1),'-.','Color',[1 0 1],'DisplayName','TSVD','LineWidth',1.5);
plot(number_doc_list,sum(time_at,1),'--','Color',[0.5 0.5 0.5],'DisplayName','AVTA+Recoverl2','LineWidth',1.5);
hold off;

legend('show','Location','east')%,'Orientation','horizontal')



%% Plot the running time
hold on;
xlabel('Number of documents') % x-axis label
ylabel('Running time(log)') % y-axis label
plot(number_doc_list,log(sum(time_aa,1)),':*','Color',[0.5 0 1],'DisplayName','Fast anchor+Recoverl2','LineWidth',1.5);
plot(number_doc_list,log(sum(time_tvt,1)),':^','Color',[0 0 1],'DisplayName','AVTA+Catch Word','LineWidth',1.5);
plot(number_doc_list,log(sum(time_tsvd,1)),'-.','Color',[1 0 1],'DisplayName','TSVD','LineWidth',1.5);
plot(number_doc_list,log(sum(time_at,1)),'--','Color',[0.5 0.5 0.5],'DisplayName','AVTA+Recoverl2','LineWidth',1.5);
hold off;

legend('show','Location','east')%,'Orientation','horizontal')
