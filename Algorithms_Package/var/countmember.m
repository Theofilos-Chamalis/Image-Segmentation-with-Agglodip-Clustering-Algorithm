%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% COUNTMEMBER - count members
%
%   C = COUNTMEMBER(A,B) counts the number of times the elements of array A are
%   present in array B, so that C(k) equals the number of occurences of
%   A(k) in B. A may contain non-unique elements. C will have the same size as A.
%   A and B should be of the same type, and can be cell array of strings.
%   unique_vals is true when all elements of A are unique and in a
%   monotonically non-decreasing order.
%
%
%   Examples:
%     countmember([1 2 1 3],[1 2 2 2 2]) 
%        -> 1     4     1     0
%     countmember({'a','b','c'},{'a','x','a'}) 
%        -> 2     0     0
%
%   See also ISMEMBER, UNIQUE, HISTC
% source: Matlab File Exchange
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function C = countmember(A,B, unique_vals)

    error(nargchk(3,3,nargin));

    if ~isequal(class(A),class(B)),
        error('Both inputs should be the same class.');
    end
    if isempty(B),      C = zeros(size(A));
        return;
    elseif isempty(A),  C = [];
        return;
    end

    if unique_vals
        % assign each element in B a number corresponding to the element of A
        %[~, L] = ismember_sorted(B,A(:));   % same as: %[~, L] = ismember(B,A(:));  but faster from lightspeed tlbx
        [~,L] = ismember(B,A(:)); % same to: L = match(B, A(:));
        % count these numbers
        C = histc(L(:),1:length(A));
        if(size(C) ~= size(A))
            C = C';
        end
    else
        % which elements are unique in A,
        % also store the position to re-order later on
        [AU,~, j] = unique(A(:)); 
        % assign each element in B a number corresponding to the element of A
        [~, L] = ismember(B,AU); 
        % count these numbers
        N = histc(L(:),1:length(AU));
        % re-order according to A, and reshape
        C = reshape(N(j),size(A));
    end
end