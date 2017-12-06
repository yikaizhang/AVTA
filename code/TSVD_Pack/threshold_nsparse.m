function [B, thres] = threshold_nsparse(A, omega1, omega2, w, d)
% Assuming input will be (d*w), that is transposed, so A will not be copied

th1 = min(floor(omega1*d/2),d);

th2 = 3*omega2*omega1*d;

nnz_1 = floor(nnz(A)/10);
thres = zeros(w,1);
id_cols = zeros(nnz_1, 1);
id_rows = zeros(nnz_1, 1);
values = zeros(nnz_1, 1);
allocs = nnz_1;

nnz_B = 0;

dw = sum(A,2); % words per doc
m = floor(mean(dw));  % average words per doc
noscale = false;
if all(dw == m), noscale = true; end

for i =1:w
    t1 = A(:,i);
    if ~noscale, t1 = round((t1 ./ dw) * m); end  % scale if m is not constant
    t = sort(t1,'descend');
    index = th1;
    
    zeta = t(index);

    while zeta ~= 0
    	if (sum(t == zeta) < th2)
            break;
    	else 
    	  index = find(t < zeta);
%           find(t < zeta)
          zeta = t(index(1));
    	end
    end
    thres(i) = zeta;
    if thres(i) <= 1
        t2 = t(t~=0);
        if ~isempty(t2)
            thres(i) = min(t2);
        end
    end
    rdocs = find(t1 >= thres(i));
    ndocs = length(rdocs);
    id_cols(nnz_B+1:nnz_B+ndocs) = rdocs;
    id_rows(nnz_B+1:nnz_B+ndocs) = ones(1,ndocs)*i;
    values(nnz_B+1:nnz_B+ndocs) = (thres(i)^0.5)*ones(1,ndocs);
    nnz_B = nnz_B + ndocs;
    if nnz_B >= allocs/2
        id_cols(nnz_B + nnz_1) = 0;
        id_rows(nnz_B + nnz_1) = 0;
        values(nnz_B + nnz_1) = 0;
        allocs = nnz_B + nnz_1;
    end
        
end
B = sparse(id_rows(1:nnz_B),id_cols(1:nnz_B),values(1:nnz_B),w,d);
B=full(B);
clear id_rows id_cols values;
end