function centdist = alldist ( centers )

% output: matrix of all pairwise distances
% input: data points (centers)

  k = size(centers,1);
  centdist = zeros(k,k);
  for j = 1:k
    centdist(1:j-1,j) = calcdist2(centers(1:j-1,:),centers(j,:));
  end
  centdist = centdist+centdist';

  return
end

