%%*******************
% Compare the robustnes of AVTA and Fast Anchor on computing 
% vertices problem.
% The comparison has # of vertices varying from 10-100
%
%
%%*******************




%% experiment setting
NUM_ITER=10;
num_vtx_vec=(1:NUM_ITER)*10;
mean_err_ta=zeros(1,NUM_ITER);
mean_err_fast=zeros(1,NUM_ITER);
time_ta=zeros(1,NUM_ITER);
time_aa=zeros(1,NUM_ITER);

for iter=1:NUM_ITER
    iter
    %% Generate data
    Num_of_vertex=num_vtx_vec(iter);
	%% vtx generate
    Num_of_points=500;
    Dim=50;
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
    %% Perturbation
    variance_n=0;
    noise_val=random('Normal',0,variance_n,m,n);
    matA_p=matA+noise_val;
    %matA_p=matA;
    
    


    %% AVTA
    t1=clock;
    index_this=AVTA_K(matA_p',K);
    t2=clock;
    
    time_ta(iter)=etime(t2,t1);
    
    recover_error_ta=zeros(1,K);
    recover_index_ta=index_this;
    recover_mat=matA_p(:,recover_index_ta);
    for i=1:K
        orin_vtx=matA(:,vertices_index(i));
        [alpha err]=frank_wolfe(orin_vtx,recover_mat,0.01);
        recover_error_ta(i)=err;
    end
    mean_err_ta(iter)=mean(recover_error_ta);
    %% Fast anchor
    t1=clock;
    recover_error_anchor=zeros(1,K);
    recover_index_anchor=fast_anchor(matA_p',K,0);
    t2=clock;
    time_aa(iter)=etime(t2,t1);
    recover_mat=matA_p(:,recover_index_anchor);
    for i=1:K
        orin_vtx=matA(:,vertices_index(i));
        [alpha err]=frank_wolfe(orin_vtx,recover_mat,0.01);
        recover_error_anchor(i)=err;
    end
    mean_err_fast(iter)=mean(recover_error_anchor);
    

end
hold on;
xlabel('# of vertices') % x-axis label
ylabel('recovery error') % y-axis label
plot(num_vtx_vec,mean_err_ta,'^-','Color',[0 0 1],'DisplayName','AVTA');
plot(num_vtx_vec,mean_err_fast,'*-','Color',[1 0 1],'DisplayName','Fast Anchor Word');
title('# of vertices vs recovery error')

hold off;

legend('show','Location','northwest')%,'Orientation','horizontal')




