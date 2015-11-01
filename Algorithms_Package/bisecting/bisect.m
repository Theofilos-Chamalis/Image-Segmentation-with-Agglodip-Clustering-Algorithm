%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This function implements the bisecting k-means clustering algorithm.
%------------
% Input parameters
% X: Data vectors (rows),
% splittrials: the number of splits to try and then ftom them select the one 
%              with minimum error
% splitMODE:   two modes are available that implement different versions of bisecting k-means. 
%              (0) the one using one random object and the center of the group to compute the second initial prototype, 
%              (1) two randomly selected objects are the initial prototypes of the split.
%              (2) this uses the unimodality test to find a good split
% kernelSpace: true/false to solve the clustering problem in kernel space using kernel k-means
%------------
% Algorithm:
%------------
% The algorithm is described by the following steps:
%     1.    Supposing we have a set M of data points 
%           (splitMODE=0): Randomly select a data point pl in space -> (left centroid candidate)
%           Compute the centroid c of data points in M
%           Compute the point pr as: pr = c-(pl-c) -> (right centroid candidate)
%           (splitMODE=1): Randomly select two data points
%           (splitMODE=2): uses the sim/dist distribution of a split voter based on Hartigan Unimodality Test
%     2.    Divide M into subclusters Ml and Mr by assigning the points
%           of M to the closest cluster candidate pl or pr
%     3.    Compute the centroids of Ml and Mr, cl and cr (as in 2)
%     4.    (convergence check) 
%           If cl == pl and cr == pr stop
%           else let cl = pl and cr = pr and goto step (2)
%
% Output parameters
% bestclmemberIDs: the ids of the objects assigned in each of the two clusters after the split
% bestgIdx:        the asobject-to-cluster assignment in a vector
% bestc:           the centroids of the two new clusters
% minEr:           the clustering error of the resulting data 2-partition
%------------
% Copyright (C) 2009-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [bestclmemberIDs, bestgIdx, bestc, minEr] = bisect (X, splitMODE, splittrials, seed, varargin)

[n, d]  = size(X);
D       = zeros(n, 2);
c       = zeros(2,d);
bestc   = zeros(2,d);
new_c   = zeros(2,d);
minclEr = flintmax;

if (splitMODE == 2)
    [found, ~, varargin] = parsepar(varargin, 'viewers_for_split');
    if (found)
        % read additional settings for split based on unimodality check
        [found, D, varargin] = parsepar(varargin, 'dmatrix');
        if (~found), error('Distance matrix not provided!'); end
    end
end

[~, kernelSpace, varargin] = parsepar(varargin, 'kernelSpace');
if (isempty(kernelSpace)), kernelSpace = 0; clear('Kernel'); end

[~, Kernel, ~] = parsepar(varargin, 'Kernel');
if (isempty(Kernel)), kernelSpace = 0; clear('Kernel'); end

clmemberIDs = cell(2,1);  % each bisecting results in a 2-partition of a cluster
clmembers   = zeros(2,1);

if (splitMODE == 0) 
    center = full( sum(X, 1) / n );
% elseif (splitMODE == 2)
% nothing
end

if (n == 2)  % singleton clusters produced
    bestclmemberIDs    = cell(2,1);
    bestclmemberIDs{1} = 1;   bestclmemberIDs{2} = 2;
    bestc = [X(1,:);  X(2,:)];
    minEr = [0; 0];
    return;
end

% the number of split trials might be smaller than the number of object in data subset (cluster)
splittrials = min(splittrials, n); 

rand ('state', seed+0);
randobjects = randperm(n);
randobjects = randobjects(1:splittrials);  % the number of split trials might be smaller than the number of object in data subset

for i=1:splittrials,    % try to split the cluster a number of times...
    % dbg: randobjects(i)
    % step 1
    if (splitMODE == 0)
        c(1,:) = full( X(randobjects(i), :) );
        c(2,:) = full( center - (c(1,:) - center) );
    elseif (splitMODE == 1) % select two different data points
        d1 = randobjects(i);
        d2 = floor(1 + rand*n);
        while (d1 == d2), d2= floor(1 + rand*n); end   % guarantee difference in seed objects
        c = full([ X(d1, :);     X(d2, :) ]);
    else%if (splitMODE == 2)
        % Choose a split viewer
        % Note: D columns are reordered wrt viewers_for_split(:,1) index ids
        % 1. choose a split viewer at random
            %distr = D(:, find(prob_viewers >= rand,1));
        % 2. choose the best split viewer
            %distr = D(:, 1);
        % 3. choose one of the splittrials best viewers
            distr = D(:, min(i,size(D,2)));
        
        mean_val = sum(distr) / length(distr);
        
        % Select a pair of objects having distance from split viewer
        % lower and greater than the average distance from split viewer
        % -> object spread lower and greater to mean by a "spread" value
        % and then optionally 1-dim k-means may applied to find better initial centroids
        % in sim/dist space and then proceed in the original space
            spread_val = std(distr);
            [startIdx, res] = ark_kmeans(distr', [], [(mean_val-spread_val) (mean_val+spread_val)], 5, -1, 1, 0, 0, 0, 0);
            [~,p] = min(abs( bsxfun(@minus, [distr distr], res) ), [], 1); % equally:  %[~,p1] = min(abs(distr - res(1)));  %[~,p2] = min(abs(distr - res(2)));
            p1 = p(1); p2 = p(2);
            
        c(1,:) = X(p1, :);
        c(2,:) = X(p2, :);
    end

    % use kernel k-means
    if (kernelSpace == 1)
        [gIdx, z] = knkmeans(Kernel,startIdx');
        z = (z + Kernel(1:size(Kernel,1)+1:end))';
        tmpind = (gIdx == 1);
        clmemberIDs{1} = find(tmpind);  clmemberIDs{2} = find(~tmpind);  % compute the 2-partition
        clmembers(1) = length(clmemberIDs{1});  clmembers(2) = length(clmemberIDs{2});
    else
        while (1)
            % step 2
            D = [sum(bsxfun(@minus, X, c(1,:)).^2, 2)    sum(bsxfun(@minus, X, c(2,:)).^2, 2)];  %alternative: D(:,1) = sum((bsxfun(@minus, X, pl)).^2, 2); D(:,2) = sum((bsxfun(@minus, X, pr)).^2, 2);
            
            [z,gIdx] = min(D,[],2);
            tmpind = (gIdx == 1);
            clmemberIDs{1} = find(tmpind);  clmemberIDs{2} = find(~tmpind);  % compute the 2-partition
            clmembers(1) = length(clmemberIDs{1});  clmembers(2) = length(clmemberIDs{2});

            % step 3: compute the mean vector
            new_c(1,:) = sum(X(clmemberIDs{1}, :),1) / clmembers(1);   new_c(2,:) = sum(X(clmemberIDs{2}, :),1) / clmembers(2); 

            % step 4
            if ( any(new_c(:) ~= c(:)) ),  c = new_c;
            else                           break;   end
        end    
    end
    
    % compute the cl. error of the two clusters and compare to the best found
    curEr{1} = z(clmemberIDs{1});
    curEr{2} = z(clmemberIDs{2});
    sum_curEr = sum(curEr{1}) + sum(curEr{2});
    
    if (sum_curEr < minclEr);  
         minEr = curEr; minclEr = sum_curEr; bestclmemberIDs = clmemberIDs; bestc = c; bestgIdx = gIdx;
    end                
end % end of trials

if (kernelSpace == 0),  minEr = [sum(sqrt(minEr{1}))  sum(sqrt(minEr{2}))];  
else                    minEr = [sum(minEr{1})        sum(minEr{2})]; 
end