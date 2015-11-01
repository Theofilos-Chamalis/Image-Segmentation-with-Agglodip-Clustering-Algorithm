%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Perform kernel k-means clustering.
%   K: kernel matrix
%   init: k (1 x 1) or label (1 x n, 1<=label(i)<=k)
% Reference: [1] Kernel Methods for Pattern Analysis
% by John Shawe-Taylor, Nello Cristianini
%------------
% Written by Michael Chen (sth4nth@gmail.com)
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [label, val, energy] = knkmeans(K,init)

    n = size(K,1);
    if length(init) == 1
        label = ceil(init*rand(1,n));
    elseif length(init) == n
        if (~isrow(init))
            label = init';
        else
            label = init;
        end
    % elseif size(init,1) == 1 && size(init,2) < n
    %     k = size(init,2);
    %     
    else
        error('ERROR: init is not valid.');
    end
    
    last = 0;
    kprev = -1;
    while any(label ~= last)
        %disp('iter');
        [u,~,label] = unique(label);   % remove empty clusters
        k = length(u);
        if (kprev > 0 && k < kprev),  fprintf('Warning: %g Clusters removed!\n', kprev-k);  end
        kprev = k;
        E = sparse(label,1:n,1,k,n,n);
        E = bsxfun(@times,E,1./sum(E,2));
        T = E*K;
        Z = repmat(diag(T*E'),1,n)-2*T;
        last = label;
        [val, label] = min(Z,[],1);
    end
    
    [~,~,label] = unique(label);   % remove empty clusters
    energy = sum(val)+trace(K);
    label = label';
end