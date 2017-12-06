function [mat_A]=Random_nn(vtx_A,n,varargin)
%******* Randomly generate data points in convex hull using vertices******
% 
% [mat_A]=Random_cvx(vtx_A,n,vargin)
%
% Inputs:
% vtx_A: dim x # of verticex  matrix. Columns as vertices
% n: number of points in the convex hull
% vargin:
% 'dir',para : Use the standard dirichlet distribution with para as
% parameter
% 
% 'unif' : Use the scaled standard uniform distribution
% 
% 'nn over cmplt': Use the unscale standard uniform distribution
%
% Output:
% mat_A: mxn matrix where columns as points.
%------------------------------------------------------------------

num_vtx=size(vtx_A,2);
if ~isempty(varargin)
    if nargin>=4
        [Type Para]=deal(varargin{:});
    else
        Type=deal(varargin{:});
        if strcmp(Type, 'dir')==1
            Para=0.1*ones(num_vtx,1);
        end
    end
else
    Type='unif';
end

if strcmp(Type, 'unif')==1
    coe_A_raw=random('unif',0,1,num_vtx,n); %generate iid gaussian 
    norm_mat_A_raw=sum(coe_A_raw,1);
    coe_A= coe_A_raw./(norm_mat_A_raw(ones(size(vtx_A,2),1),:));%scale points to uniball
elseif strcmp(Type, 'dir')==1
    coe_A=dirichletrnd(Para, n); %generate weights of vertices using dirichlet distribution
elseif strcmp(Type, 'nn over cmplt')==1
    coe_A=random('unif',0,1,num_vtx,n); %generate iid gaussian 
end

mat_A=vtx_A*coe_A;


