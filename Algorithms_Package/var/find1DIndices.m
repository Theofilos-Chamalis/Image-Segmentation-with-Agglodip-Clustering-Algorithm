%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Finds all the indices in A that contains a value included in T.
% A is the array to search into
% T is the array of target values (non-decreasing order)
% example:
%   find1DIndices([1 2 2 1 1 4 5 1 3 4 2], [1 3])
%------------
% Copyright (C) Argyris Kalogeratos, February 2012.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function index = find1DIndices(A, T)

    [r,c] = find( bsxfun(@eq, A(:), T) );
    index = mat2cell(r, histc(T(c), T), 1);
end