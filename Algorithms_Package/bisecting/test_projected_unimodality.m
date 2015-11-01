%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Tests a cluster for split based on the dip-dist criterion which uses 
% the dip-statistic on the projections of the data of a cluster.
%------------
% Input
%  projections : a cell containing the projections of the dataset (row vectors)
%  nboot : the number of bootstrap uniform distributions to use
%  nboot_accel : 1 or 0. If 1 we calculate one bootstrap sample (offline or online)
%  and then use it for all the projections' p_value calculations by dip test else we 
%  assign this job to the HartigansDipSignifTest routine to do these calculations
%  simultaneously. Also we make use of all the cores of the underlying machine 
%  using parallel parfor execution which is especially useful in large
%  datasets.
%
% Output
%   maxdip: the maximum dip value found in the projections
%   pmin  : the lowest probability found by the test (if pmin==0 cluster is definitely not unimodal)
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



function [maxdip, pmin] = test_projected_unimodality (projections, nboot, nboot_accel, boot_dips)
    
    n = size(projections,2);
    p_value = zeros(1,n);
    dip = 100*ones(1, n);

    if(nboot_accel == 1)
        N = length(projections{1});
        ext = exist('boot_dips','var');
        % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf 
        % of sample size N (the same as empirical pdf)
        % if ext = 1, we have calculated this offline
        if (ext == 0)
            boot_dip = zeros(nboot, 1);
            parfor i=1:nboot
               unifpdfboot = sort(unifrnd(0,1,1,N));
               unif_dip    = HartigansDipTest(unifpdfboot);
               boot_dip(i) = unif_dip;
            end;
            boot_dip = sort(boot_dip);
        elseif(ext == 1)
            boot_dip = boot_dips(N,:);
            boot_dip = boot_dip';
        end

        parfor i=1:n,
            [dip(i),p_value(i)] = HartigansDipSignifTest_no_boot(projections{i},nboot,boot_dip);
        end
    else
        parfor i=1:n,
            [dip(i),p_value(i)] = HartigansDipSignifTest(projections{i},nboot);
        end
    end

    [maxdip,i] = max(dip);
    pmin = p_value(i) ;

end
