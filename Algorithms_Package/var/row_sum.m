%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% ROW_SUM   Sum for each row.
% A faster and more readable alternative to sum(x,2).
% source: Lightspeed package
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function s = row_sum(x)

    % unfortunately, this removes any sparseness of x.
    s = x*ones(cols(x),1);
end

