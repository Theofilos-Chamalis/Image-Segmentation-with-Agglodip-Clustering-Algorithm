%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This function implements the bisecting k-means clustering algorithm.
% On top of it, x-means, g-means, dip-means and pdip-means
% algorithms are implemented as wrappers that each uses a different 
% statistical criterion to decide for a cluster split. 
%------------
% Input parameters
% X:          Data vectors (rows), 
% split_struct:information for the algorithms that learn the number of k.
%              For dip-means it must contain the following fields
%                - pval_threshold (default=0)       the probability <= to which the cluster must split
%                - exhaustive_search (default=1)  whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
%                - voting (default=0.1)                the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
%                - nboot  (default=1000)             number of uniform distribution samples to test each object with Hartigan test
%                - overall_distr (default=0)           if 1, Hartigan test is applied one time on the overall distribution of distances
%              For g-means it contains only one field: 
%                - pval_threshold (default=0.95)  which is the critical value for the internal normality statistic test 
%              When split_struct is not provided then default values are set for all fields.
% splittrials: the number of splits to try and then from them select the one 
%              with minimum error (default: 1)
% splitSELECT: this is the strategy according to which a cluster is selected to be split. 
%              (0) the largest error,
%              (1) the largest average error
%              (2) the largest cluster (default)
%              (3) dip-means: the cluster with lowest prob of unimodal sim/dist distribution (Hartigan)
%              (4) x-means with BIC criterion
%              (5) g-means with Anderson-Darling normality test
%              (6) pdip-means: the cluster with lowest prob of unimodal data distribution of its line projections(Hartigan)
% splitMODE:   two modes are available that implement different versions of bisecting k-means. 
%              (0) the one using one random object and the center of the group to compute the second initial prototype (default), 
%              (1) two randomly selected objects are the initial prototypes of the split.
%              (2) the split is based on Hartigan unimodality  test (can be used only in dip-means)
% refineMODE:  this is the refinement strategy
%              (0) no refinement (default)
%              (1) one k-means refinement at the end of the bisecting recursion
%              (2) refinement of all clusters after every bisection of one cluster
% kernelSpace:      true/false to solve the clustering problem in kernel space using kernel k-means
% attempts:          attempts to cluster the data (starting again from the single cluster case).
% smallest_cluster: the algorithm stops splitting clusters containing this number of objects (default=6) 
% rndseed:            initialization for the random number generator (dafault: random from cpu timer)
%
% Output
% R:              the partition of data
% min_err:      the numerical error corresponding to R partition
% R_ref:         the partition of data if refinement is applied 
% min_err_ref: the numerical error corresponding to the R_ref partition
%------------
% Copyright (C) 2012-2015, Chamalis Theofilos, Kalogeratos Argyris.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, min_err, R_ref, min_err_ref] = bisect_kmeans_default (X, varargin)

    % load the uniform distribution bootstrap samples that were generated
    % offline which will be used by the dip-means and pdip-means algorithms
    n = size(X,1);
    load('unifpdfbootext.mat','boot_dips');

    % -- parse input arguments --%
    [found, splittrials, varargin] = parsepar(varargin, 'split_trials');
    if (~found), splittrials = 1; end

    [found, splitSELECT, varargin] = parsepar(varargin, 'splitSELECT');
    if (~found), splitSELECT = 2; end

    [found, splitMODE, varargin] = parsepar(varargin, 'splitMODE');
    if (~found), splitMODE = 0; end

    [found, refineMODE, varargin] = parsepar(varargin, 'refineMODE');
    if (~found), refineMODE = 0; end

    [found, smallest_cluster, varargin] = parsepar(varargin, 'smallest_cluster');
    if (~found), smallest_cluster = 6; end   % the cluster must have at least this number of objects to be split candidate

    [found, attempts, varargin] = parsepar(varargin, 'attempts');
    if (~found), attempts = 1; end

    [~, split_struct, varargin] = parsepar(varargin, 'split_struct');

    [found, exps_group_id, varargin] = parsepar(varargin, 'rndseed');
    if (~found), exps_group_id = rand('state',sum(100*clock)); end

    [kernelSpace, Kernel, ~] = parsepar(varargin, 'kernelSpace');
    if (isempty(Kernel)), kernelSpace = 0; clear('Kernel'); end


    if (splitSELECT == 3 || splitSELECT == 6) % Hartigan unimodal test on similaritites/distances or on the projections of the data
        if (~isempty(split_struct))
            pval_threshold    = split_struct.pval_threshold;         % the probability <= to which the cluster must split
            exhaustive_search = split_struct.exhaustive_search; % whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
            voting            = split_struct.voting;                       % the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
            nboot             = split_struct.nboot;                       % number of uniform distribution samples to test each object with Hartigan test
            overall_distr     = split_struct.overall_distr;              % if 1, Hartigan test is applied one time on the overall distribution of distances
            clear('split_struct');
        else                        % default values
            pval_threshold    = 0;    
            exhaustive_search = 1;    
            voting            = 0.1;  
            nboot             = 1000; 
            overall_distr     = 0;    
        end

        if (kernelSpace == 0 && splitSELECT == 3)% || (kernelSpace == 1 && splitMODE == 2)) % no kernel space is to be used
            D = sqdist_rows(X);                                     % squareform(pdist(X, 'euclidean'));
        elseif(kernelSpace == 0 && splitSELECT == 6)        
            proj_Cell = {};  % proj_Cell is a cell containg all the projections of the data
        else
            splitMODE = 2; % force this type of split mode
            D = Kernel; clear('Kernel');
            selfSim = D(1:size(D,1)+1:end);
        end
    else
        % revert splitMode to default if not dip-means is requested by user
        if (splitMODE == 2), splitMODE = 0; end  

        if (splitSELECT == 4)

        elseif (splitSELECT == 5)
            if (~isempty(split_struct))
                 pval_threshold = split_struct.pval_threshold;
            else
                pval_threshold = 0.95;  % critical value for g-means
            end
        end
    end


    if (splitSELECT > 2), kmax = 1; end
    if (~kernelSpace && (splitSELECT == 3 || splitSELECT == 4 || splitSELECT == 6)), c_all = ComputeCentroids(X, ones(n,1), 1, 0); end

    min_err = flintmax;
    min_err_ref = flintmax;

    for exp_id=1:attempts,  % do a number of clustering attempts starting each time from the single cluster case
        % at first all data points are assigned to one cluster
        gIdx = ones(n,1);

        clmemberIDs = cell(kmax,1);    clmemberIDs{1} = (1:n);
        clmembers   = zeros(kmax,1);   clmembers(1) = n;

        Er = zeros(kmax,1);
        Er(1) = flintmax;       %the first cluster to be split
        k = 1;                     %the number of clusters up to now
        if (~kernelSpace && (splitSELECT == 3 || splitSELECT == 4 || splitSELECT == 6))
            c = c_all;               
        else
            c = [];
        end

        while (k < kmax || splitSELECT > 2)
            bestclmemberIDs = []; bestsplit_c = [];
            candidates = find(clmembers > smallest_cluster);

            %choose the cluster to split
            switch splitSELECT
            case 0
                [~, cltosplit] = max(Er(candidates));
                cltosplit = candidates(cltosplit);
            case 1
                [~, cltosplit] = max(Er(candidates)./double(clmembers(candidates)));
                cltosplit = candidates(cltosplit);
            case 2 
                [~, cltosplit] = max(clmembers);
            case 3
                if (isempty(candidates)), 
                    fprintf('-> Hartigan: Stopped splitting at k=%g  (smallest clusters reached)\n', k); 
                    break;
                end
                cluster_eval = ones(k,1); % a value >= (max prob = 1)
                pval         = ones(k,1);

                if (overall_distr == 0) % use individual object-viewers
                    if (splitMODE == 2)
                        viewers_for_split = cell(k,1);
                        for i=candidates,
                            [cluster_eval(i), pval(i), ~, viewers_for_split{i}] = test_unimodal_cluster (D(clmemberIDs{i},clmemberIDs{i}), nboot, exhaustive_search, voting);
                        end
                    else 
                        for i=candidates,
                            [cluster_eval(i), pval(i)] = test_unimodal_cluster (D(clmemberIDs{i},clmemberIDs{i}), nboot, exhaustive_search, voting,boot_dips);
                        end
                    end
                else
                    for i=candidates,
                        DD = tril(D(clmemberIDs{i},clmemberIDs{i}));
                        DD = DD(:);  DD = DD(DD~=0);
                        [cluster_eval(i), pval(i)] = HartigansDipSignifTest (DD, nboot); %cluster_eval = dip value
                    end
                end
                %[maxv, cltosplit] = max(cluster_eval);
                [~, cltosplit] = max(cluster_eval);            % select cluster with minimum p_value -> less probable to be unimodal wrt similarities/distributions (min(p_value) -> max(cluster_eval))

                minpv = pval(cltosplit);

                if (minpv > pval_threshold), 
                    fprintf('-> Hartigan: Stopped splitting at k=%g\n', k); 
                    break; 
                end
            case 4
                if (isempty(candidates)),
                    fprintf('-> BIC: Stopped splitting at k=%g  (smallest clusters reached)\n', k); 
                    break;
                end
                cluster_eval = -ones(k,1);
                max_improvement = 0;
                cltosplit = -1;

                for i=candidates,
                    % split cluster using function bisect(...)
                    [splitclmemberIDs, split_Idx, split_c, Error] = bisect(X(clmemberIDs{i},:), splitMODE, splittrials, exps_group_id*exp_id*k);
                    % BIC: Apply Bayesian Information Criterion to determine which whether the
                    % generated child clusters are more representative of the real distribution
                    parent_score = bic(X(clmemberIDs{i},:), ones(clmembers(i),1), c(i,:));
                    child_score  = bic(X(clmemberIDs{i},:), split_Idx, split_c);

                    if (child_score > parent_score)
                        cluster_eval(i) = abs(child_score - parent_score);
                        if (cluster_eval(i) > max_improvement)
                            max_improvement = cluster_eval(i);
                            bestclmemberIDs = splitclmemberIDs;
                            bestsplit_c     = split_c;
                            newEr           = Error;
                            cltosplit       = i;
                        end
                    end
                end

                if (cltosplit == -1), 
                    fprintf('-> BIC: Stopped splitting at k=%g\n', k); 
                    break; 
                end
            case 5
                if (isempty(candidates)), 
                    fprintf('-> And-Darl: Stopped splitting at k=%g  (smallest clusters reached)\n', k); 
                    break;
                end
                
                cluster_eval = -ones(k,1);
                pad = ones(k,1);
                maxad = 0;
                cltosplit = -1;
                
                for i=candidates,
                    % calculate the anderson-darling statistic
                    tcov = cov(X(clmemberIDs{i},:));  % first get the principal component of the group
                    if ( sum(sum(tcov ~= 0)) == 0 )    % deal with zero covariance matrix, consider this a group
                        pad(i) = 1; cluster_eval(i) = -1;
                        continue;
                    else
                        [pc, ~, ~] = pcacov(tcov);
                        projected  = pc(:,1)'*X(clmemberIDs{i},:)';                                                               % project the points on the principal component
                        [pad(i), cluster_eval(i)] = anderson_darling( projected, pval_threshold, smallest_cluster ); % calculate the statistic along the component
                    end

                    % check if this cluster is not accepted
                    if ( pad(i) == 0 )
                        % when tad is -2 it means that the distribution was far from normal (log(0) error in anderson_darling.m)
                        if ( cluster_eval(i) == -2 )
                            cltosplit = i;
                            break;
                        elseif ( cluster_eval(i) > maxad )
                            maxad     = abs(cluster_eval(i));
                            cltosplit = i;
                        end
                    end
                end

                if (cltosplit == -1), 
                    fprintf('-> And-Darl: Stopped splitting at k=%g\n', k); 
                    break; 
                end
              case 6
                    if (isempty(candidates)), 
                        fprintf('-> Hartigan: Stopped splitting at k=%g  (smallest clusters reached)\n', k); 
                        break;
                    end
                    
                    cluster_eval = ones(k,1); % a value >= (max prob = 1)
                    pval = ones(k,1);

                   if (splitMODE == 2)
                        fprintf('\nsplitMODE = 2 is not supported by pdip-means. Use 0 or 1 instead.\n\n');
                        break;
                    else 
                        for i=candidates,
                            % data_Proj outputs the projections of the dataset X row-wise in to the cell called
                            % proj_Cell. It takes as input: the data, 3 triggers(pca_proj, axis_proj, randproj
                            % with a value of 1 or 0) and randproj_number.
                             proj_Cell = data_Proj(X(clmemberIDs{i},:),1,1,0,18);
                             [cluster_eval(i), pval(i)] = test_projected_unimodality(proj_Cell, nboot, 1, boot_dips);  %cluster_eval = dip_value                                          
                        end
                    end
                    %[maxv, cltosplit] = max(cluster_eval); 
                    [~, cltosplit] = max(cluster_eval);   % select cluster with minimum p_value -> less probable to be unimodal wrt similarities/distributions (min(p_value) -> max(cluster_eval))

                    pvaltosplit = pval(cltosplit);

                    if (pvaltosplit > pval_threshold),  %if pval_threshold = 0, end the splits and continue to find the error
                        fprintf('-> Projected Hartigan: Stopped splitting at k=%g\n', k); 
                        break; 
                    end                
            end

            newcl = k+1; % the id of the new clusters (the pre-exists)

            % split cluster using function bisect(...)
            if (isempty(bestclmemberIDs)) % bestclmemberIDs exists only for BIC (x-means) case
                if (splitMODE == 2)
                     [bestclmemberIDs, ~, bestsplit_c, newEr] = bisect(X(clmemberIDs{cltosplit},:), splitMODE, splittrials, exps_group_id*exp_id*k, 'kernelSpace', kernelSpace, 'Kernel', D(clmemberIDs{cltosplit}, clmemberIDs{cltosplit}), 'viewers_for_split', viewers_for_split{cltosplit}, 'dmatrix', D(clmemberIDs{cltosplit}, viewers_for_split{cltosplit}(:,1)));
                else
                     [bestclmemberIDs, ~, bestsplit_c, newEr] = bisect(X(clmemberIDs{cltosplit},:), splitMODE, splittrials, exps_group_id*exp_id*k);
                end
            end

            % set correct cluster assignements for data points based on the 2-partition
            clmemberIDs{newcl}           = clmemberIDs{cltosplit}(bestclmemberIDs{2});
            clmemberIDs{cltosplit}       = clmemberIDs{cltosplit}(bestclmemberIDs{1});
            clmembers([cltosplit newcl]) = [length(clmemberIDs{cltosplit}) length(clmemberIDs{newcl})];
            k = k + 1;
            
            % compute the centroids of the bisecting k-means algorithm
            if (isempty(bestsplit_c)),   % bestsplit_c exists only for BIC (x-means) case        
                for t=[cltosplit newcl],
                    c(t,:) = sum(X(clmemberIDs{t},:),1) / clmembers(t); 
                end
            else
               c([cltosplit newcl],:) = bestsplit_c;
            end

            if (refineMODE == 2 && k > 2)
                fprintf('** Internal refinement...\n');
                % run a refinement phase with a k-sp version
                if (kernelSpace == 0)
                    [gIdx, c, ~, ~, ~, Ei] = ark_kmeans(X', [], c', 5, -1, 1, 0, 0, 0, 0);
                    c = c';
                else
                    for t=1:k, gIdx(clmemberIDs{t}) = t; end
                    [gIdx, Ei] = knkmeans(D, startIdx);
                     Ei = (Ei + Kernel(1:size(Kernel,1)+1:end))';
                     c = [];
                end
                clmemberIDs = find1DIndices(gIdx(:), 1:k);
                clmembers = cellfun(@length, clmemberIDs)';
                Er = cellfun(@(x) sum(Ei(x)), clmemberIDs);   % new error for each cluster
            else
                % new error for each cluster
                Er([cltosplit  newcl]) = newEr;
            end
        end

        % if no split occured, then calculate the overall error of the dataset
        if (k == 1)
            if (isempty(c)), c = ComputeCentroids(X, ones(n,1), 1, 0); end
            Er(1) = sum(sqrt( sum(bsxfun(@minus, X, c).^2, 2) ));
        end

        sumer = sum(Er);
        for j=1:k,  gIdx(clmemberIDs{j}) = j; end

        %%%%%%%%%%%%%%%%%%%%%%%%%%
        % REFINEMENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        if (refineMODE == 1)
            fprintf('** Refining final solution...\n');        
            % compute the centroids of the bisecting k-means algorithm
            if (kernelSpace == 0)
                c = ComputeCentroids(X, gIdx, k, 0);
                % run a refinement phase with a k-means version
                [gIdx_ref, c, sumd_ref, ~, ~] = ark_kmeans(X', [], c', 5, -1, 1, 0, 0, 0, 0);
                c = c';
            else
                [gIdx_ref, Ei, sumd_ref] = knkmeans(D, gIdx);
                Ei = (Ei + selfSim)';
                c = [];
            end
        else
            gIdx_ref = gIdx; sumd_ref = sumer;
        end

        if (sumer < min_err)
            min_err = sumer;
            R = gIdx;
        end
        
        if (sumd_ref < min_err_ref)
            min_err_ref = sumd_ref;
            R_ref = gIdx_ref;
        end 
        
     end  % end of attempts
 
end