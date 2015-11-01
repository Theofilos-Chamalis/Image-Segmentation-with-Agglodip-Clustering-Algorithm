%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Computes pairwise squared Euclidean distances between points
% considers data vectors as rows of a and b
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function d = sqdist_rows(a,b)

    if (nargin == 2) 
        aa = sum(a.*a, 2); bb = sum(b.*b, 2); ab = a * b';
        d = abs(repmat(aa,[1 size(bb,1)]) + repmat(bb',[size(aa,1) 1]) - 2*ab);
        d(1:size(d,1)+1:end) = 0;    
        d = sqrt(d);
    else
        aa = sum(a.*a, 2); ab = a * a';
        d = abs(repmat(aa,[1 size(aa,1)]) + repmat(aa',[size(aa,1) 1]) - 2*ab);
        d(1:size(d,1)+1:end) = 0;    
        d = sqrt(d);
    end
    % for column vectors use this
    %aa = sum(a.*a); bb = sum(b.*b); ab = a'*b; 
    %d = abs(repmat(aa',`[1 size(bb,2)]) + repmat(bb,[size(aa,2) 1]) - 2*ab);
end
