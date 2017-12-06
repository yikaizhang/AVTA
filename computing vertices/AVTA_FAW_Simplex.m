%%*******************
% Compare the robustnes of AVTA and Fast Anchor on computing perturbed
% vertices problem.
% The comparison has variance scale varying from 0.3-3. 
% Experimental setting:
%	Simplex : 	Num_of_vertex=100;
%			Num_of_points=500;
%			Dim=100;
%
%
%%*******************


Num_of_vertex=100;
Num_of_points=500;
Dim=100;
K=Num_of_vertex;



%% experiment setting
NUM_ITER=10;
Rep=10;
variance_vec=(1:NUM_ITER)*0.3;
mean_err_ta=zeros(1,NUM_ITER);
mean_err_fast=zeros(1,NUM_ITER);
time_ta=zeros(1,NUM_ITER);
time_aa=zeros(1,NUM_ITER);

for iter=1:NUM_ITER
    for rep=1:Rep
        iter
        %% vtx generate

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
        variance_n=variance_vec(iter);
        noise_val=random('Normal',0,variance_n,m,n);
        matA_p=matA+noise_val;
        %matA_p=matA;

        t1=clock;

        NUM_OF_GAMMA=20;
        number_of_vertices=zeros(1,NUM_OF_GAMMA);

        gamma_val=zeros(1,NUM_OF_GAMMA);
        now_gamma=5;
        stop_rule=K;
        length_use=0;




    %%AVTA
        index_this=AVTA_MultipleRp(matA_p',K,20,K);
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
end


hold on;
xlabel('perturbation variance') % x-axis label
ylabel('recovery error') % y-axis label
plot(variance_vec,mean_err_ta,'^-','Color',[0 0 1],'DisplayName','AVTA');
plot(variance_vec,mean_err_fast,'*-','Color',[1 0 1],'DisplayName','Fast Anchor Word');
title('perturbation variance vs recovery error')

hold off;

legend('show','Location','northwest')%,'Orientation','horizontal')

% 
% hold on;
% plot(variance_vec,time_ta,'Color',[0 0 1],'DisplayName','Triangle Algorithm');
% plot(variance_vec,time_aa,'--','Color',[1 0 1],'DisplayName','FastAnchorWord');
% title('perturbation variance vs running time')
% 
% hold off;
% 
% legend('show','Location','west')%,'Orientation','horizontal')
% 
% savefile = ['perturbation_tests'  datestr(now,30) '.mat']
% save(savefile, 'mean_err_fast', 'mean_err_ta','time_ta','time_aa');
% 
% 
% 
