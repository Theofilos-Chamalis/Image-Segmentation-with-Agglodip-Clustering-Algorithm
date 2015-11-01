%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Tests a cluster for split based on the dip-dist criterion which uses 
% the dip-statistic on the distributions of similarities/distances of a 
% viewer point in cluster to all other members in that cluster.
%------------
% Input
%   D:       the dataset (row vectors)
%   nboot: the number of bootstrap uniform distributions to use
%
% Output
%   maxdip: the maximum dip value found in the data cluster
%   pmin:    the lowest probability found by the test (if pmin==0 cluster is definitely not unimodal)
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxdip, pmin] = test_unimodal_cluster2 (D, nboot, boot_dips)
    n = size(D,1);
    n1 = size(D,2);

    dip     = 100*ones(n, 1); % a large init value
    ext = exist('boot_dips','var');
    
    if(ext==1)
        boot_dip = zeros(1, nboot);
        %boot_dip = boot_dips(N,:);
        boot_dip = boot_dips(n1,:);
        boot_dip = boot_dip';
    else
        error('There was a problem loading offline bootstrap samples. Fix this and execute again. Quitting now...');
    end

    p_value = zeros(n,1);

    for i=1:n,
       dip(i) = HartigansDipTest(D(i,:));
       p_value(i) = sum(dip(i) < boot_dip) / nboot;
    end

    [maxdip,i] = max(dip);
    pmin = p_value(i);

end


