%%***************This program compares AVTA quick hull on vertices computing problem********************
% Compare the efficiency of AVTA and the quick hull algorithm
% on vertices computing problem. The comparison has dimension varing from
% 2~12.
%%***********************************

%% Parameter
Num_of_vertex=100;
Num_of_points=500;
Dim=50;

NUM_ITER=11;
dim_vec=(1:NUM_ITER)+1;

%% Experiment setting
K=Num_of_vertex;
mean_err_ta=zeros(1,NUM_ITER);
mean_err_fast=zeros(1,NUM_ITER);
time_ta=zeros(1,NUM_ITER);
time_qhull=zeros(1,NUM_ITER);
vtx_num_ta=zeros(1,NUM_ITER);
vtx_num_qhull=zeros(1,NUM_ITER);
for iter=1:NUM_ITER
    for jj=1:20
        iter
        dim_i=dim_vec(iter);


        normal_val=Random_pts(Dim,Num_of_vertex,'normal')*10;


        %% Generate data in convex hull
        inhulldata=Random_nn(normal_val,Num_of_points,'unif');
        matA=[inhulldata,normal_val];
        [m,n]=size(matA);
        normal_val_index=(Num_of_points+1):n;

        %% Mix data
        rndindx=randperm(n,n);
        [or_val or_ind]=sort(rndindx);
        vertices_index=or_ind(normal_val_index);
        matA=matA(:,rndindx);
    
        %%  Initialization
        epsilon=0.001;
        alpha=0;
        tof=0;
        matA_p=matA;
        
        %% AVTA
        t1=clock;
        index_this=AVTA_eps(matA_p',0.05);
        t2=clock;

        time_ta(iter)=time_ta(iter)+etime(t2,t1);
        vtx_num_ta=length(index_this);
        
        %% Quickhull
        t1=clock;
          QH = convhulln(matA_p');
        t2=clock;
        time_qhull(iter)=time_qhull(iter)+etime(t2,t1);
        vtx_num_qhull=length(unique(QH));
    end
    

end
% hold on;
% plot(dim_vec,mean_err_ta,'Color',[0 0 1],'DisplayName','Triangle Algorithm');
% plot(dim_vec,mean_err_fast,'--','Color',[1 0 1],'DisplayName','FastAnchorWord');
% title('perturbation variance vs recovery error(1000 reperturbation)')
% 
% hold off;
% 
% legend('show','Location','northwest')%,'Orientation','horizontal')


hold on;
plot(dim_vec,time_ta,'Color',[0 0 1],'DisplayName','AVTA');
plot(dim_vec,time_qhull,'--','Color',[1 0 1],'DisplayName','Quickhull');
title('dim vs running time')

hold off;

legend('show','Location','west')%,'Orientation','horizontal')

savefile = ['perturbation_tests'  datestr(now,30) '.mat']
save(savefile, 'mean_err_fast', 'mean_err_ta','runningtime');



