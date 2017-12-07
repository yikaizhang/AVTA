%%***************This program compares AVTA with frank_wolfe , Triangle algorithm and simplex  method********************
% Compare the efficiency of AVTA, frank_wolfe, Traignle algorithm and the simplex method on 
% solving convexhull query problem. The comparison has number of redudant
% points(points in the convex hull) varying from 5000-500000.
%%***********************************
%% ---Parameters
Num_iter=10; %Number of different parameters
Num_of_vertex=100;
Num_rep=8;
Num_of_points=500000;
Dim=50;
Num_pts=((1:Num_iter).^2)*5000;


%% Initialization
lp_time=zeros(Num_rep,Num_iter);
frank_time=zeros(Num_rep,Num_iter);
ta_time=zeros(Num_rep,Num_iter);
avta_time=zeros(Num_rep,Num_iter);

for iter=1:Num_iter
    iter
    
    Num_of_points=Num_pts(iter);
	parfor rep=1:Num_rep
    %for rep=1:Num_rep
        %% Generate vertices using high dimensional gaussian
        normal_val=Random_pts(Dim,Num_of_vertex,'normal');



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

        query_point=Random_nn(normal_val,1,'unif');%% The query point is in the convex hull

        norm_matA=sum(matA.*matA);
        max_norm=sqrt(max(norm_matA));	
        error_eps=0.0002*max_norm;
	
        %% Frank wolfe
        tic
        [alpha err]=frank_wolfe(query_point,matA,error_eps);
        frank_t=toc;
        frank_time(rep,iter)=frank_t;
        
        %% AVTA
        tic
        [index query_result p_p alpha]=AVTA_query(matA',query_point,error_eps);

        avtat=toc;
        avta_time(rep,iter)=avtat;

        %% Triangle algorithm
        tic
        [inorout,p_prime,alpha_coe,dist]=Tri_Algo(matA,query_point,error_eps);
        tat=toc;
        
        ta_time(rep,iter)=tat;

        %% Linear programming
        n_var=n;
        n_con=m;
        x_y_size=n_con+n_var;
        cvx_aeq=ones(1,n_var);
        Aeq=[matA,eye(n_con);ones(1,n_var),zeros(1,n_con)];
        beq=[query_point;1];
        lb=zeros(x_y_size,1);
        tic
        [bb,fval] = linprog([zeros(n_var,1);ones(n_con,1)],[],[],Aeq,beq,lb,[]);
        lp_t=toc;
        lp_time(rep,iter)=lp_t;
     end
end



%% Plot the running time
hold on;

plot(Num_pts,sum(frank_time,1),'+','Color',[0 1 0],'DisplayName','FW');
plot(Num_pts,sum(ta_time,1),'Color',[0 0 1],'DisplayName','TA');
plot(Num_pts,sum(avta_time,1),'--','Color',[1 0 1],'DisplayName','AVTA');
plot(Num_pts,sum(lp_time,1),'-.','Color',[0.5 0.5 0.5],'DisplayName','lp');
title('Running time vs #of points( p in convex hull)')

hold off;


legend('show','Location','northwest')%,'Orientation','horizontal')


% 
% 
% %% Plot the running time
% hold on;
% 
% plot(Num_pts,log(frank_time),'+','Color',[0 1 0],'DisplayName','FW');
% plot(Num_pts,log(ta_time),'Color',[0 0 1],'DisplayName','TA');
% plot(Num_pts,log(avta_time),'--','Color',[1 0 1],'DisplayName','AVTA');
% plot(Num_pts,log(lp_time),'-.','Color',[0.5 0.5 0.5],'DisplayName','lp');
% title('log(Running time) vs #of points( p in convex hull)')
% 
% hold off;
% 
% 
% legend('show','Location','northwest')%,'Orientation','horizontal')
% 
