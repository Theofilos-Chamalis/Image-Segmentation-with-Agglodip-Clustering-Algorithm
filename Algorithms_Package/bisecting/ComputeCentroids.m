%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Computes the centroids based on a partition table.
%------------
% Input
%   X:      the dataset (row vectors)
%   gIdx:  the matrix that contains the object-to-cluster assingment
%   k:      the number of cluster
%   mode: if mode=1 then normalization is performend wrt norm-2.
% Output
%   c:                the computed centroids
%   clmemberIDs: a cell with the object ids assigned to each cluster.
%   clmembers  :  the size of each cluster 
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c, clmemberIDs, clmembers] = ComputeCentroids(X, gIdx, k, mode)

if (iscell(gIdx))
    clmemberIDs = gIdx;
    for t=1:k,
        clmembers(t) = length(clmemberIDs{t});
    end
else
    clmemberIDs = cell(k,1);
    clmembers   = zeros(k,1);
    for t=1:k,
        clmemberIDs{t} = find(gIdx==t);
        clmembers(t) = length(clmemberIDs{t});
    end
end

c = zeros(k,size(X,2));
for t=1:k,
    c(t,:) = sum(X(clmemberIDs{t},:),1) / clmembers(t); 
end

% normalize initial prototypes for cosine simfunc
if (mode == 1),  c = bsxfun(@rdivide, c, sum(c.^2, 2).^(0.5)); end 
