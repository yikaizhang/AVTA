%%*******************
% Compare the efficiency of AVTA and Fast Anchor on computing 
% vertices problem.
% 
% Experimental setting:
%	Simplex : 	Num_of_vertex=100;
%               Num_of_points=500;
%               Dim=100;
%
%
%%*******************

%% vtx generate
Num_of_vertex=50;
Num_of_points=1000;
Dim=100;
K=Num_of_vertex;
normal_val=Random_pts(Dim,Num_of_vertex,'normal')*10;


%% Generate interior pts
inhulldata=Random_nn(normal_val,Num_of_points,'unif');
matA=[inhulldata,normal_val];
[m,n]=size(matA);
normal_val_index=(Num_of_points+1):n;

rndindx=randperm(n,n);
[or_val or_ind]=sort(rndindx);
vertices_index=or_ind(normal_val_index);
matA=matA(:,rndindx);


%% experiment setting
NUM_ITER=5;

mean_err_ta=zeros(1,NUM_ITER);
mean_err_fast=zeros(1,NUM_ITER);
time_ta=zeros(1,NUM_ITER);
time_aa=zeros(1,NUM_ITER);

for iter=1:NUM_ITER
    iter
%%AVTA
    t1=clock;
    index_this=AVTA_K(matA_p',K);
    t2=clock;
    
    time_ta(iter)=etime(t2,t1);
%%Fast anchor
    t1=clock;
    recover_error_anchor=zeros(1,K);
    recover_index_anchor=fast_anchor(matA_p',K,0);
    t2=clock;
    time_aa(iter)=etime(t2,t1);
end

hold on;
xlabel('different datasets') % x-axis label
ylabel('running time') % y-axis label
plot(1:NUM_ITER,time_ta,'Color',[0 0 1],'DisplayName','Triangle Algorithm');
plot(1:NUM_ITER,time_aa,'--','Color',[1 0 1],'DisplayName','FastAnchorWord');
title('different datasets vs running time')
hold off;
legend('show','Location','west')%,'Orientation','horizontal')

savefile = ['perturbation_tests'  datestr(now,30) '.mat']
save(savefile, 'mean_err_fast', 'mean_err_ta','time_ta','time_aa');



