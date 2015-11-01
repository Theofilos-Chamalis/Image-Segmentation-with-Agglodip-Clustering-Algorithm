%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% seeks for a specific variable name in a par list (case-insesnsitive)
% source: Matlab File Exchange
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [found, val, vars] = parsepar(vars, param)

    isvar = cellfun(@(x) ischar(x) && strcmpi(x, param), vars);

    if (sum(isvar) > 1), error('Parameters can only be passed once'); end

    if any(isvar)
        found = true;
        idx = find(isvar);
        val = vars{idx+1};
        vars([idx idx+1]) = [];
    else
        found = false;
        val = [];
    end
end
