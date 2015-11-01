%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Computes the RBF kernel (similarity matrix) for the row vectors of X and 
% with a provided sigma parameter.
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function K = RBFkernel (X, sigma)

    K = exp(-sqdist_rows(X).^2 ./ (2*sigma^2)); % rbf kernel
end