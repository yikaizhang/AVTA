function [index_set]=vtx_prune(superset_data,K,varargin)
%%********Given a super set of vertices, prune and output K vertices**********
%Input:
%       superset_data : dim x # of vertices. Columns as vertices. 
%       K: # of vertices to output.
%       varargin: 'forward': incremental, 'backward': decremental.
%Output: 
%       index_set: Logical variable. If i-th entry is '1' then i-th point
%       is output as vertex ('0' otherwise).
%%******************
[NN,MM]=size(superset_data);

if ~isempty(varargin)
    [Type]=deal(varargin{:});
else
    Type='backward';
end

if MM>K
    if strcmp(Type, 'forward')==1
        index_set=logical(zeros(1,MM));
        norm_vtx=sqrt(sum(superset_data.^2,1));
        max_ind=min(find(norm_vtx==max(norm_vtx)));
        index_set(max_ind)=1;
        index_series=1:MM;
        while sum(index_set)<K
            sum(index_set)
            diff_set=logical(abs( index_set-1));
            current_indexset=index_series(diff_set);
            dist_vtx=zeros(1,length(current_indexset));
            for i=1:length(current_indexset)
                p=superset_data(:,current_indexset(i));
                current_cols=superset_data(:,index_set);
                [alp dis]=frank_wolfe_2(p, current_cols, 0.01);
                dist_vtx(i)=dis;
            end
             max_dist_index=min(find(dist_vtx==max(dist_vtx)));
            index_set(current_indexset(max_dist_index))=1;
        end
       
    else
        %disp('start prune');
        [NN MM]=size(superset_data);
        index_set=logical(ones(1,MM));
        index_series=1:MM;
        while sum(index_set)>K
            %disp(['current number of vertices is ' num2str(sum(index_set))])
            current_indexset=index_series(index_set);
            dist_vtx=zeros(1,length(current_indexset));

            for i=1:length(current_indexset)
                p=superset_data(:,current_indexset(i));
                tmp_ind=setdiff(current_indexset,current_indexset(i));
                current_cols=superset_data(:,tmp_ind);
                [alp dis]=frank_wolfe_2(p, current_cols, 0.01);
                %dis
    %             sum(p)
                dist_vtx(i)=dis;
            end
            min_dist_index=min(find(dist_vtx==min(dist_vtx)));
            index_set(current_indexset(min_dist_index))=0;
        end
    end
else
    index_set=logical(ones(1,MM));
end
