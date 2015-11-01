%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Standardize a dataset of real vectors.
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function X = standardize(X, dim)

    n = size(X,dim);
    X = bsxfun(@minus, X, sum(X,dim) / n);
    X = bsxfun(@rdivide, X, std(X, 0, dim));
end
