%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Finds all the index maps in A that contains a value included in T.
% A is the array to search into
% T is the array of target values (non-decreasing order)
% example:
%   find1DIndexMaps([1 2 2 1 1 4 5 1 3 4 2], [1 3])
%   [1 0 0 1 1 0 0 1 0 0 0
%    0 0 0 0 0 0 0 0 1 0 0]'
%------------
% Copyright (C) Argyris Kalogeratos, February 2012.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indexMap = find1DIndexMaps(A, T)

    indexMap = bsxfun(@eq, A(:), T);
end
