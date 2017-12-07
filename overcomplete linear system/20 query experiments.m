%%***************This program compares AVTA with simplex  method on the linear system feasibility problem********************
% Compare the efficiency of AVTA+ simplex and the simplex method on non negative
% linear system feasibility problem. In the experiment, the total number of
% query is 20.
%%***********************************
%% ---Parameters

NUM_query=(1:20);


%% Initialization 
time_lp=zeros(1,length(NUM_query));
time_avta_lp=zeros(1,length(NUM_query));
%% Generate data

Num_of_vertex=100;
Dim=100;

error_eps=0.01;
Num_pts=50000;
A=rand(Dim,Num_of_vertex);
mat_1=Random_cvx(A,Num_pts,'nn over cmplt')*2;
mat_A=[mat_1,A];


[M_con,N_var]=size(mat_A);
mat_A=mat_A(:,randperm(N_var,N_var));



%% Dimension reduction using AVTA
tic
rd_vec=random('normal',0,1,M_con,1);
rd_vec=rd_vec./(sqrt(sum(rd_vec.^2)));
inner_val=rd_vec'*mat_A;
b=1;
mat_AA=zeros(M_con,N_var);
for i=1:N_var
    tmp_vec=mat_A(:,i)./(sqrt(sum(mat_A(:,i).^2)));
    scale_val=b/(rd_vec'*tmp_vec);
    mat_AA(:,i)=scale_val*tmp_vec;
end

min_d=sqrt(min(sum(mat_AA.^2)));
[index]=AVTA_eps(mat_AA',min_d*0.1)
t_avta_lp=toc;
time_avta_lp=time_avta_lp+t_avta_lp;
mat_a=mat_A(:,index);

%% Compute feasibility for multiple queries
for iter=1:length(NUM_query)
    
    %% Generate feasible query points 
    if iter>=length(NUM_query)*0.5
        X=rand(N_var,1)*5;
        b=mat_A*X;
    else
        b=rand(Dim,1);
    end
    
    %% Simplex method
    [n_1,n_2]=size(mat_A);
    x_y_size=n_1+n_2;

    c_coe=[zeros(n_2,1);ones(n_1,1)];
    Aeq=[mat_A eye(n_1)];

    beq=b;
    lb=zeros(x_y_size,1);
    tic
    [bb,fval] = linprog(c_coe,[],[],Aeq,beq,lb,[]);
    t_lp=toc;
    time_lp(iter:end)=time_lp(iter:end)+t_lp;
    
    %% Simplex method on reduced problem
    [n_1,n_2]=size(mat_a);
    x_y_size=n_1+n_2;
    c_coe=[zeros(n_2,1);ones(n_1,1)];
    % Aeq=;
    Aeq=[mat_a eye(n_1)];

    beq=b;
    lb=zeros(x_y_size,1);
    tic
    [bb,fval] = linprog(c_coe,[],[],Aeq,beq,lb,[]);
    t_av_lp=toc;
    time_avta_lp(iter:end)=time_avta_lp(iter:end)+t_av_lp;
    
end

hold on;

plot(NUM_query,time_lp,'Color',[0 0 1],'DisplayName','LP');
plot(NUM_query,time_avta_lp,'--','Color',[1 0 1],'DisplayName','AVTA+LP');
title('Running time vs Size of the problem(Feasible)')

hold off;


legend('show','Location','northwest')%,'Orientation','horizontal')
