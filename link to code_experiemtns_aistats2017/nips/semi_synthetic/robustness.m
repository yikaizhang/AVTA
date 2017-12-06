%%*******************************************************
%   This program runs the robustness experiment.
%%*******************************************************
%parpool(4); If parallel is used, uncomment second loop->parfor re=1:rep.
load('nips_2000_wd.mat');% load the word topic distribution
rep=10;% repeat t times. 
sample_num=10000;
var_iter=(1:10)*0.005;
tsvd_err=zeros(rep,length(var_iter)); %error of tsvd
tvt_err=zeros(rep,length(var_iter));  %error of AVTA+Catch
aa_err=zeros(rep,length(var_iter));   %error of l2+Fast anchor
at_err=zeros(rep,length(var_iter));   %error of l2+AVTA
time_tsvd=zeros(1,length(var_iter));  %running time of tsvd
time_tvt=zeros(1,length(var_iter));   %running time of AVTA+Catch
time_aa=zeros(1,length(var_iter));    %running time of l2+Fast anchor
time_at=zeros(1,length(var_iter));    %running time of l2+AVTA
currentFolder = pwd;  %use current file to wirte results
savefile = [pwd datestr(now,30) 'rbst.mat'];
writefile=pwd;

for iter=1:length(var_iter)
    %parfor re=1:rep
    for re=1:rep
        disp(['iter=',num2str(iter)]);
        disp(['re=',num2str(re)]);
        K=50;
        scale=var_iter(iter);
        word_num=1000;
        disp('generate data');
        dataset=gibb_art(scal,50,word_num,sample_num);
        dataset=dataset+rand(length(dataset(:,1)),sample_num)*scale;
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
            orin_vtx=scal(:,i);
            diffmat=tvt-orin_vtx(:,ones(K,1));
            distance_v=sum(abs(diffmat),1);
            error_tvt(i)=min(distance_v);
        end
        tvt_err(re,iter)=mean(error_tvt);
        error_tvd=zeros(1,K);
        for i=1:K 
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
save(savefile,'number_doc_list','tsvd_err','tvt_err','tvt_K','aa_err','time_tvt','time_tsvd','time_aa','time_at');


%% Plot the error

hold on;
xlabel('Number of documents') % x-axis label
ylabel('Running time') % y-axis label
plot(var_iter,sum(time_aa,1),':*','Color',[1 0 1],'DisplayName','Fast Anchor+Recoverl2','LineWidth',1.5);
plot(var_iter,sum(time_tvt,1),':^','Color',[0 0 1],'DisplayName','AVTA+Catch Word','LineWidth',1.5);
plot(var_iter,sum(time_tsvd,1),'-.','Color',[0.5 0.5 0.5],'DisplayName','TSVD','LineWidth',1.5);
plot(var_iter,sum(time_at,1),'--','Color',[0.5 0.5 0.5],'DisplayName','AVTA+Recoverl2');

hold off;

legend('show','Location','east')%,'Orientation','horizontal')

%% Plot the running time
hold on;
xlabel('Scale of Noise') % x-axis label
ylabel('L1 err') % y-axis label
plot(var_iter,sum(tsvd_err,1)./rep,'--','Color',[1 0 1],'DisplayName','TSVD');
plot(var_iter,sum(tvt_err,1)./rep,'^-','Color',[0 0 1],'DisplayName','AVTA+Catchword');
plot(var_iter,sum(aa_err,1)./rep,'*-','Color',[0.5 0 1],'DisplayName','Fast Anchor+Recoverl2');
plot(var_iter,sum(at_err,1)./rep,'--','Color',[0.5 0.5 0.5],'DisplayName','AVTA+Recoverl2');
hold off;
legend('show','Location','east')%,'Orientation','horizontal')

%% Plot the range of l1-error
x = var_iter;
y_tvt = sum(tvt_err)./rep;
y_tvt_pos = max(tvt_err)-y_tvt;
y_tvt_neg=y_tvt-min(tvt_err); 
y_tsvd = sum(tsvd_err)./rep;
y_tsvd_pos = max(tsvd_err)-y_tsvd;
y_tsvd_neg=y_tsvd-min(tsvd_err); 
y_aa = sum(aa_err)./rep;
y_aa_pos = max(aa_err)-y_tsvd;
y_aa_neg=y_tsvd-min(tsvd_err); 
y_at = sum(at_err)./rep;
y_at_pos = max(at_err)-y_tsvd;
y_at_neg=y_tsvd-min(tsvd_err); 
hold on
xlabel('scale of noise') % x-axis label
ylabel('err range in repeat experiments ') % y-axis label
myC= [0 0 1
  1 0 0
  1 0.4 0
  0 0.8 1
  0.6 0 1
  0 1 0 ];
H=bar(x,[y_tvt_pos+y_tvt_neg;y_tsvd_pos+y_tsvd_neg;y_aa_pos+y_aa_neg;y_at_pos+y_at_neg]');
for k=1:length(H)
  set(H(k),'facecolor',myC(k,:))
end
AX=legend(H, {'AVTA+Catch Word','TSVD','Fast Anchor+Recoverl2','AVTA Recoverl2'}, 'Location','Best','FontSize',8);
hold off;





