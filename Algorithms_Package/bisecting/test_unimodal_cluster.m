%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Tests a cluster for split based on the dip-dist criterion which uses 
% the dip-statistic on the distributions of similarities/distances of a 
% viewer point in cluster to all other members in that cluster.
%------------
% Input
%   D:       the dataset (row vectors)
%   nboot: the number of bootstrap uniform distributions to use
%   exhaustive_search: whether all viewers should be checked or to use some
%                              other faster way to find the first viewer that
%                              rejects unimodality.
%   voting: the number of viewers that are required to indicate
%             multimodality in order to decide for a non-unimodal cluster 
%             that should be split.
%
% Output
%   maxdip: the maximum dip value found in the data cluster
%   pmin:    the lowest probability found by the test (if pmin==0 cluster is definitely not unimodal)
%   distr:    the similarity/distance distribution of the object with maxdip and minp
%   viewers_for_split: is a two-column matrix where in the first column all
%                           the ids of the objects that propose the cluster split 
%                           are stored in a descending order wrt the dip value. The 
%                           second column contains those dip values.
%------------
% Copyright (C) 2012-2015, Chamalis Theofilos, Kalogeratos Argyris.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [maxdip, pmin, distr, viewers_for_split] = test_unimodal_cluster (D, nboot, exhaustive_search, voting,boot_dips)
    n = size(D,1);

    dip     = 100*ones(n, 1); % a large init value
    p_THR   = 0; 
    ext = exist('boot_dips','var');
    
    if(ext==0)
        % calculate a bootstrap sample of size NBOOT of the dip statistic for a uniform pdf of sample size N (the same as empirical pdf)
        boot_dip = zeros(nboot, 1);
        unifpdfboot = sort(unifrnd(0, 1, n, nboot), 1, 'ascend');
        parfor i=1:nboot
           boot_dip(i) = HartigansDipTest(unifpdfboot(:,i)); % check the p_value of a uniform distribution
        end;
        boot_dip = sort(boot_dip, 'ascend');
    elseif(ext==1)
        % use the external bootstrap samples
        boot_dip = zeros(1, nboot);
        %boot_dip = boot_dips(N,:);
        boot_dip = boot_dips(n,:);
        boot_dip = boot_dip';
    end


    % enable exhaustive_search
    if (voting ~= 0), exhaustive_search = 1; end

    if (exhaustive_search == 0)
        p_value = 100*ones(n,1);  % a large init value

        % Setup the order in which the objects (view points) should be checked for
        % unimodality based on their distance to other cluster objects.
        % -> objects with higher distances lay out of the core mode of the cluster  
        % and are probably better viewpoints!
        [~, check_order] = sort(sum(D,1), 'descend');

        for i=check_order,
           dip(i) = HartigansDipTest(D(:,i));
           p_value(i) = sum(dip(i) < boot_dip);

           if (p_value(i) == 0),
               maxdip = dip(i); pmin = 0; distr = D(:,i);
               return;
           end
        end
       
        [maxdip,i] = max(dip);
        pmin       = p_value(i) / nboot;
        
    else % if exhaustive_search == 1
        p_value = zeros(n,1);

        parfor i=1:n,
           dip(i) = HartigansDipTest(D(:,i));
           p_value(i) = sum(dip(i) < boot_dip) / nboot;
        end
  
        if (voting == 0)
            [maxdip,i] = max(dip);
            pmin = p_value(i) ;
            distr = D(:,i);
        else
            viewers_for_split = find(p_value <= p_THR);
            split_votes = length(viewers_for_split);
            pmin = split_votes < voting*n; % check if the 'want to split' objects are the majority
            [maxdip_val,i] = max(dip);
            distr          = D(:,i);
            if (split_votes > 0)
                 maxdip = sum(dip(p_value <= p_THR)) / split_votes; % return average dip value between the points that vote for split
            else
                maxdip = maxdip_val;
            end

            if (nargout == 4)
                [sorted_viewers, ind] = sort(dip(viewers_for_split), 'descend');
                viewers_for_split = [viewers_for_split(ind), sorted_viewers, p_value(viewers_for_split(ind))];
            end
        end
    end
end


