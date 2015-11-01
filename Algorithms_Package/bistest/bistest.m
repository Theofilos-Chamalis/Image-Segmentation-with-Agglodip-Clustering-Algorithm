%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Agglodip hierarchical algorithm image segmentation
% demonstration script.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% remove all warnings if any, found at old matlab versions
w = warning ('off','all');

% find all windows of type figure (if any), and close them
delete(findall(0,'Type','figure'))

% define the RNG seed  
rseed = 10;

% the real number of cluster/segments is not known beforehand
real_k = -1;
[N,d] = size(X);  
DATASIZE = N;

% set the number of initial clusters for agglodip 
numberOfInitialClusters = numberOfInitialClustersEXT;

% set the minimum number of elements a cluster can have
smallestCl = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run algorithm
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

methods = 1;

method_name = '   Agglodip';

split_struct = cell(1,1);

split_struct{1} = struct;
split_struct{1}.pval_threshold    = 0.00;    % the probability <= to which the cluster must split
split_struct{1}.exhaustive_search = 1;      % whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
split_struct{1}.voting            = 0;           % the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
split_struct{1}.nboot             = 1000;      % number of normal distribution samples to test each object with Hartigan test

merge_trials = 10;                                  % times to try a merge (we use 10 for more stable results)

result = zeros(max(methods), 7);
m = methods;

fprintf('\n++++++++++++++++++++++++++++++AGGLODIP-SEGMENTATION++++++++++++++++++++++++++++++\n');
tic;
[R, sumer, R_ref, sumer_ref, realClustFinal] = bisect_agglodip_segmentation(X, gIdx_ref_Init,c_all, Am1, numberOfInitialClusters,'merge_struct', split_struct{m}, 'merge_trials', merge_trials, 'mergeSELECT', 1, 'splitMODE', 0, 'refineMODE', 0, 'smallest_cluster', smallestCl, 'attempts', 1, 'rndseed', 0+rseed);
toc;
disp(' ');  

k = length(unique(R_ref));
result(m, 1:3) = [k, sumer, sumer_ref];

if (real_k > 0)
    [pq, RI, ARI, conf_matrix, conf_matrix_probC, conf_matrix_probR] =  partition_quality(C,R_ref);
    VI = varinfo(C,R_ref);
    result(m, 4:6) = [RI, ARI, VI];
end       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Show results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('----------------------------------------------------------------------------------------------------------------------------------------\n\t\t\t\tClustering results for real_k = %g\n----------------------------------------------------------------------------------------------------------------------------------------\n', real_k);
fprintf('%g. %10s -- k: %3g, error: %5.5f, (supervised measures N/A)\n', m, method_name, result(m,1), result(m,3));

fprintf('RNG seed used: %f\n\n\n', rseed);