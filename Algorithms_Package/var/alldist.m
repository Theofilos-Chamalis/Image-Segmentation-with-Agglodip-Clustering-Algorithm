%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% input: data points (centers)
% output: matrix of all pairwise distances
% source: fastkmeans package from Matlab File Exchange
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function centdist = alldist(centers)

    k = size(centers,1);
    centdist = zeros(k,k);
    for j = 1:k
        centdist(1:j-1,j) = calcdist(centers(1:j-1,:),centers(j,:));
    end
    centdist = centdist+centdist';
end
