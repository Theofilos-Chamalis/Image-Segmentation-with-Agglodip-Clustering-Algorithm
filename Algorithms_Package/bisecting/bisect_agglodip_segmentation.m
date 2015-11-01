%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This function implements the agglodip (agglomerative dip-means) hierarchical clustering algorithm.
% This version of the algorithm is changed to accomodate the needs of image segmentation by
% using the neighbors of each superpixel for merging and other minor changes.
%------------
% Input parameters
% X:            Data vectors (rows), 
% merge_struct:information for the algorithms that learn the number of k.
%                - pval_threshold (default=0)     the probability <= to which the cluster would split
%                - exhaustive_search (default=1)  whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
%                - voting (default=0.1)           the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
%                - nboot  (default=1000)          number of normal distribution samples to test each object with Hartigan test
%               When merge_struct is not provided then default values are set for all fields.
% mergetrials: the number of merges to try and then from them select the one 
%               with minimum error (default: 1)
% mergeSELECT: this is the strategy according to which a cluster is selected to be merged. 
%               (1) agglodip: the cluster with lower prob of unimodal sim/dist distribution on the line projections(Hartigan)
% splitMODE:   one mode is available in bisect_agglodip. 
%               (0) the one using one random object and the center of the group to compute the second initial prototype (default), 
% refineMODE:  this is the refinement strategy
%               (0) no refinement (default)
%               (1) one k-means refinement at the end of the merging recursion
% attempts:          attempts to cluster the data (starting again from the initial clusters that k-means created in the beginning)
% smallest_cluster: the algorithm stops splitting clusters containing this number of objects (default=6) 
% rndseed:            initialization for the random number generator (dafault: random from cpu timer)
%
% Output
% R:              the partition of data
% min_err:      the numerical error corresponding to R partition
% R_ref:         the partition of data if refinement is applied 
% min_err_ref: the numerical error corresponding to the R_ref partition
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [R, min_err, R_ref, min_err_ref, realClustFinal] = bisect_agglodip_segmentation (X, gidxrefinit, c_all, Am1, numberOfInitialClusters, varargin)

    % load the uniform distribution bootstrap samples that were generated
    % offline which will be used by the agglodip algorithm
    n = size(X,1);
    load('unifpdfbootext.mat','boot_dips');

    initclNo = numberOfInitialClusters;
    fprintf('\nThe initial number of clusters is: %d\n',initclNo);
    
    % find the neighbors of each superpixel
    Am1Mat = zeros(length(Am1));
    for j1=1:length(Am1)
       neighbors = find(Am1(:,j1) == 1); 
       for j2=1:length(neighbors)
            Am1Mat(j1,neighbors(j2)) = 1;
            Am1Mat(neighbors(j2),j1) = 1; 
       end
    end
    
    GRAPH_FLAG = 0;                   % Use the connected components of the graph produced by the pval values of the first iteration to cluster the data
    CENTROID_TO_ALL_FLAG = 0;   % Use only the centroids of each of the 2 clusters to calculate the distances between them and the rest of the 2 clusters data
	LOCAL_SEARCH_FLAG = 1;        % Use to perform local search for merging a superpixel only with one of its neighbors
        
    gIdx_ref_Init = gidxrefinit;
    kmax = initclNo + 1;

    % -- parse input arguments -- %
    [found, mergetrials, varargin] = parsepar(varargin, 'merge_trials');
    if (~found), mergetrials = 1; end

    [found, mergeSELECT, varargin] = parsepar(varargin, 'mergeSELECT');
    if (~found), mergeSELECT = 1; end

    [found, splitMODE, varargin] = parsepar(varargin, 'splitMODE');
    if (~found), splitMODE = 0; end

    [found, refineMODE, varargin] = parsepar(varargin, 'refineMODE');
    if (~found), refineMODE = 0; end    
    if(refineMODE == 2)
        fprintf('\nrefineMODE can be set only to 1 in agglomerative versions.\nSetting refineMODE = 1...\n\n');
        refineMODE = 1;
    end

    [found, smallest_cluster, varargin] = parsepar(varargin, 'smallest_cluster');
    if (~found), smallest_cluster = 8; end   % the cluster must have at least this number of objects

    [found, attempts, varargin] = parsepar(varargin, 'attempts');
    if (~found), attempts = 1; end

    [~, merge_struct, varargin] = parsepar(varargin, 'merge_struct');

    [found, exps_group_id, varargin] = parsepar(varargin, 'rndseed');
    if (~found), exps_group_id = rand('state',sum(100*clock)); end

    [kernelSpace, Kernel, ~] = parsepar(varargin, 'kernelSpace');
    if (isempty(Kernel)), kernelSpace = 0; clear('Kernel'); 
    else
        kernelSpace = 0; clear('Kernel'); fprintf('\n\nKernel cannot be used in dip-means reverse. Clearing Kernel now and continuing... \n\n');
    end


    if (~isempty(merge_struct))
        pval_threshold    = merge_struct.pval_threshold;            % the probability <= to which the cluster must split
        exhaustive_search = merge_struct.exhaustive_search;    % whether Hartigan test must be done for all cluster objects or to stop when first prob=0 is found
        voting            = merge_struct.voting;                          % the spliting criterion is based on a voting (0<voting<=1), or on the worst indication from the objects of the cluster (voting=0)
        nboot             = merge_struct.nboot;                          % number of normal distribution samples to test each object with Hartigan test
        clear('merge_struct');
    else                        % default values
        pval_threshold    = 0;    
        exhaustive_search = 1;    
        voting            = 0.1;  
        nboot             = 1000; 
    end

    % calculate variables that will be used later to find
    % which cluster has moved to which of the final ones
    clmemberIDsStart = cell(kmax,1);
    for clindexStart = 1:n
        if isempty(clmemberIDsStart{gIdx_ref_Init(clindexStart)})
           clmemberIDsStart{gIdx_ref_Init(clindexStart)} = clindexStart;
        else
           clmemberIDsStart{gIdx_ref_Init(clindexStart)} = horzcat(clmemberIDsStart{gIdx_ref_Init(clindexStart)},clindexStart);
        end
    end

    clmemberIDsStart = clmemberIDsStart';
    firstClusterMember = zeros(1,initclNo);

   for i=1:initclNo
           if(~isempty(clmemberIDsStart{:,i}))
               firstClusterMember(i) = clmemberIDsStart{:,i}(1);
           end
   end

   clear clmemberIDsStart clindexStart; 
       
   Am1Matrix = Am1Mat;
   kmax = initclNo;
    
    % Unimodality tests are done on each cluster we have up to now
    % to check if we need to split any of them
    %[gIdx_ref_Init,c_all,initclNo,kmax] = UnimodalityChecking(X, gIdx_ref_Init, initclNo, smallest_cluster, mergeSELECT);
    
    if CENTROID_TO_ALL_FLAG == 0
        D = sqdist_rows(X);
    end
    
    min_err = flintmax;
    min_err_ref = flintmax;
    
    for exp_id=1:attempts,  % do a number of clustering attempts starting with the clusters that the initial k-means created

        gIdx = gIdx_ref_Init;
        clmemberIDs = cell(kmax,1);    
        
        for clindex = 1:n
            if isempty(clmemberIDs{gIdx(clindex)})
               clmemberIDs{gIdx(clindex)} = clindex;
            else
               clmemberIDs{gIdx(clindex)} = horzcat(clmemberIDs{gIdx(clindex)},clindex);
            end
        end

        clmemberIDs = clmemberIDs';
        clmembers   = zeros(1,kmax);   
        
        for clMembNo = 1:initclNo
            clmembers(clMembNo) = length(clmemberIDs{clMembNo});
        end

        Er = zeros(1,kmax);
        k = initclNo;                    % the number of clusters up to now
        
        if (mergeSELECT == 1)
            c = c_all;
        else
            c = [];
        end

        FIRST_RUN_FLAG = 1;       % Flag used to distinguish the first run of the merging part from the others

        while (k < kmax || mergeSELECT < 2)
                candidates = find(clmembers >= 1);

            % choose the cluster to split
            switch mergeSELECT
              case 1
                    if (isempty(candidates)), 
                        fprintf('-> Hartigan: Stopped merging at k=%g  (smallest clusters reached)\n', k); 
                        break;
                    end

                   if (splitMODE == 1 || splitMODE == 2)
                       fprintf('\nsplitMODE = 1 or 2 is not supported by dip-means reverse. Use 0 instead.\n\n');
                       break;
                   else                  
                       clear tempclmemberIDs;
                       if (FIRST_RUN_FLAG==1)
                           cluster_eval = ones(k); % a value >= (max prob = 1)
                           pval = zeros(k);
                           PropMatrix = zeros(k);
                         
                           for i=candidates
                                for j= i + 1 : length(candidates)
                                    
                                    tempclmemberIDs = [clmemberIDs{i},clmemberIDs{j}]; %Create a temporary cluster made from the members of cluster i & j
                                    
                                    if(clmembers(i)>=clmembers(j))
                                        N1 = clmembers(j);
                                        N2 = clmembers(i);
                                    else
                                        N1 = clmembers(i);
                                        N2 = clmembers(j);
                                    end 
                                    
                                    if CENTROID_TO_ALL_FLAG == 1 %Centroids to all case
                                        Xtemp = X(tempclmemberIDs(:),:);
                                        Xtemp1 = X(clmemberIDs{i},:);
                                        Xtemp2 = X(clmemberIDs{j},:);
                                        Xtlength = size(Xtemp,1);
                                        Xtlength1 = size(Xtemp1,1);
                                        Xtlength2 = size(Xtemp2,1);
                                        gIdXtemp1 = ones(Xtlength1,1);
                                        gIdXtemp2 = ones(Xtlength2,1);
                                        ctemp1 = ComputeCentroids(Xtemp1, gIdXtemp1, 1, 0);
                                        ctemp2 = ComputeCentroids(Xtemp2, gIdXtemp2, 1, 0);
                                        Dtemp1 = [];

                                        for ti = 1:Xtlength
                                            Dtemp1 = [Dtemp1 pdist([ctemp1; Xtemp(ti,:)])];
                                        end

                                        Dtemp2 = [];
                                        for ti2 = 1:Xtlength
                                            Dtemp2 = [Dtemp2 pdist([ctemp2; Xtemp(ti2,:)])];
                                        end

                                        Dtemp = [Dtemp1;Dtemp2];
                                        if(LOCAL_SEARCH_FLAG == 1)
                                            if(Am1Matrix(i,j)==1)
                                                [cluster_eval(i,j), pval(i,j)] = test_unimodal_cluster2 (Dtemp, nboot, boot_dips);
                                            end
                                        else
                                           [cluster_eval(i,j), pval(i,j)] = test_unimodal_cluster2 (Dtemp, nboot, boot_dips); 
                                        end
                                        clear Xtemp Xtemp1 Xtemp2 Xtlength Xtlength1 Xtlength2 gIdXtemp gIdXtemp1 gIdXtemp2 ctemp ctemp1 ctemp2 Dtemp Dtemp1 Dtemp2;                                   
                                    else
                                        if(LOCAL_SEARCH_FLAG == 1)
                                            if(Am1Matrix(i,j)==1)
                                                [cluster_eval(i,j), pval(i,j)] = test_unimodal_cluster (D(tempclmemberIDs,tempclmemberIDs), nboot, exhaustive_search, voting, boot_dips);
                                            end
                                        else
                                            [cluster_eval(i,j), pval(i,j)] = test_unimodal_cluster (D(tempclmemberIDs,tempclmemberIDs), nboot, exhaustive_search, voting, boot_dips);
                                        end
                                    end
                                    
                                    if pval(i,j)>0
                                        PropMatrix(i,j) = N1/N2; %Keep the number of points in the cluster if pvalue is positive
                                    else
                                        PropMatrix(i,j) = 0;
                                    end
                                end
                                FIRST_RUN_FLAG = 0;
                           end
                           
                           if GRAPH_FLAG == 1 %Graph case
                               InitialPvalMatrix = pval;
                               InitialDipMatrix = cluster_eval;
                               GraphMatrix = zeros(k);
                               Gedges = 0;
                               
                               for i = 1:k
                                   for j = i:k
                                       if(InitialPvalMatrix(i,j)~=0) % if pvalue is positive create an edge in the GraphMatrix
                                           GraphMatrix(i,j) = 1;
                                           GraphMatrix(j,i) = 1;
                                           Gedges = Gedges +1;
                                           InitialPvalMatrix(j,i) = InitialPvalMatrix(i,j);
                                       end
                                   end
                                   InitialPvalMatrix(i,i) = 1;
                               end
                               
                               SparseGraph = sparse(GraphMatrix);
                               [NoOfFinalClusters,CompAssignment] = graphconncomp(SparseGraph);
                
                               %Assignments of the new members to one cluster and
                               %elimination of the other one that is merged
                                for i=1:length(CompAssignment)
                                        if(CompAssignment(i)~=i)
                                            clmemberIDs{CompAssignment(i)} = [clmemberIDs{CompAssignment(i)},clmemberIDs{i}];
                                            clmemberIDs{i} = [];
                                            clmembers(CompAssignment(i)) = clmembers(CompAssignment(i)) + clmembers(i);
                                        end
                                end
                               
                               clmemberIDs = clmemberIDs(~cellfun('isempty',clmemberIDs));
                               clmembers = clmembers(1:NoOfFinalClusters);
                               k = NoOfFinalClusters;
                               
                               for j2=1:k,  gIdx(clmemberIDs{j2}) = j2; end
                               c = ComputeCentroids(X, gIdx, k, 0);
                               break;
                           end                           
                       else %2nd stage of algorithm, in 1st we compute all pairs of clusters to get the pvalue, now just the new (merged) with the old ones
                              cluster_eval2 = ones(1,length(candidates));
                              pval2 = zeros(1,length(candidates));
                              PropMatrix2 = zeros(1,length(candidates));
                              
                              for i=candidates
                                  
                                    tempclmemberIDs = [clmemberIDs{best_cl},clmemberIDs{i}];
                                    
                                    if(clmembers(best_cl)>=clmembers(i))
                                        N1 = clmembers(i);
                                        N2 = clmembers(best_cl);
                                    else
                                        N1 = clmembers(best_cl);
                                        N2 = clmembers(i);
                                    end
                                    
                                    if CENTROID_TO_ALL_FLAG == 1
                                        Xtemp = X(tempclmemberIDs(:),:);
                                        Xtemp1 = X(clmemberIDs{i},:);
                                        Xtemp2 = X(clmemberIDs{best_cl},:);
                                        Xtlength = size(Xtemp,1);
                                        Xtlength1 = size(Xtemp1,1);
                                        Xtlength2 = size(Xtemp2,1);
                                        gIdXtemp1 = ones(Xtlength1,1);
                                        gIdXtemp2 = ones(Xtlength2,1);
                                        ctemp1 = ComputeCentroids(Xtemp1, gIdXtemp1, 1, 0);
                                        ctemp2 = ComputeCentroids(Xtemp2, gIdXtemp2, 1, 0);
                                        Dtemp1 = [];

                                        for ti = 1:Xtlength
                                            Dtemp1 = [Dtemp1 pdist([ctemp1; Xtemp(ti,:)])];
                                        end

                                        Dtemp2 = [];
                                        for ti2 = 1:Xtlength
                                            Dtemp2 = [Dtemp2 pdist([ctemp2; Xtemp(ti2,:)])];
                                        end

                                        Dtemp = [Dtemp1;Dtemp2];
                                        if(LOCAL_SEARCH_FLAG == 1)
                                            if(Am1Matrix(i,best_cl)==1)
                                                [cluster_eval2(i), pval2(i)] = test_unimodal_cluster2 (Dtemp, nboot, boot_dips);
                                            end
                                        else
                                            [cluster_eval2(i), pval2(i)] = test_unimodal_cluster2 (Dtemp, nboot, boot_dips);
                                        end
                                       clear Xtemp Xtemp1 Xtemp2 Xtlength Xtlength1 Xtlength2 gIdXtemp gIdXtemp1 gIdXtemp2 ctemp ctemp1 ctemp2 Dtemp Dtemp1 Dtemp2; 
                                    else
                                        if(LOCAL_SEARCH_FLAG == 1)
                                            if(Am1Matrix(i,best_cl)==1)
                                                [cluster_eval2(i), pval2(i)] = test_unimodal_cluster (D(tempclmemberIDs(:),tempclmemberIDs(:)), nboot, exhaustive_search, voting, boot_dips);
                                            end
                                        else
                                            [cluster_eval2(i), pval2(i)] = test_unimodal_cluster (D(tempclmemberIDs(:),tempclmemberIDs(:)), nboot, exhaustive_search, voting, boot_dips);
                                        end
                                        clear tempclmemberIDs;
                                    end              
                                    
                                    if pval2(i)>0
                                        PropMatrix2(i)=N1/N2;
                                    else
                                        PropMatrix2(i)=0;
                                    end
                                    
                                    if i>best_cl
                                        cluster_eval(best_cl,i) = cluster_eval2(i);
                                        pval(best_cl,i) = pval2(i);
                                        PropMatrix(best_cl,i) = PropMatrix2(i);
                                        cluster_eval(i,best_cl)=1;
                                        pval(i,best_cl)=0;
                                        PropMatrix(i,best_cl)=0;
                                    elseif i==best_cl
                                        cluster_eval(i,i)=1;
                                        pval(i,i)=0;
                                        PropMatrix(i,i)=0;
                                    else
                                        cluster_eval(i,best_cl) = cluster_eval2(i);
                                        pval(i,best_cl) = pval2(i);
                                        PropMatrix(i,best_cl) = PropMatrix2(i);
                                        cluster_eval(best_cl,i) = 1;
                                        pval(best_cl,i) = 0;
                                        PropMatrix(best_cl,i) = 0;
                                    end
                              end
                      end
                   end                       
                    
                   [propvalues,prIndex1] = max(PropMatrix.*pval); %get the two clusters with the biggest product of pvalue * (size1/size2)
                   [maxprop,prIndex2] = max(propvalues);
                   cltomerge1 = prIndex1(prIndex2);
                   cltomerge2 = prIndex2;
                   
                   if(cltomerge2>cltomerge1(1))
                      temp = cltomerge2;
                      cltomerge2 = cltomerge1(1) ;
                      cltomerge1(1) = temp;
                      cltomerge1(cltomerge2) = temp;
                   else
                      cltomerge1(cltomerge2) = cltomerge1;
                   end
                   
                   dipval = cluster_eval(cltomerge2,cltomerge1(1));

                    if (maxprop == pval_threshold) %if maxprop reaches the threshold 0
							fprintf('-> Hartigan: Stopped merging at k=%g\n', k); 
							break;						
                    end
            end

            %Assignments of the new members to one cluster and
            %elimination of the other one that is merged
            clmemberIDs{cltomerge2}=[clmemberIDs{cltomerge2},clmemberIDs{cltomerge1(cltomerge2)}];
            clmembers(cltomerge2) = clmembers(cltomerge2) + clmembers(cltomerge1(cltomerge2));

            k = k-1;

            c(cltomerge2,:) = sum(X(clmemberIDs{cltomerge2},:),1) / clmembers(cltomerge2);
            c(cltomerge1(cltomerge2),:) = [];
            clmemberIDs{cltomerge1(cltomerge2)} = [];
            clmembers(cltomerge1(cltomerge2)) = [];
            cluster_eval(cltomerge1(cltomerge2),:) = [];
            cluster_eval(:,cltomerge1(cltomerge2)) = [];
            pval(cltomerge1(cltomerge2),:) = [];
            pval(:,cltomerge1(cltomerge2)) = [];
            PropMatrix(cltomerge1(cltomerge2),:) = [];
            PropMatrix(:,cltomerge1(cltomerge2)) = [];
            clmemberIDs = clmemberIDs(~cellfun('isempty',clmemberIDs));
            
            if(LOCAL_SEARCH_FLAG == 1)
                %Inheriting neighbors to the new merged cluster
                for i3=1:size(Am1Matrix,2)
                    if(Am1Matrix(cltomerge1(cltomerge2),i3)==1)
                        Am1Matrix(cltomerge1(cltomerge2),i3)=0;
                        Am1Matrix(i3,cltomerge1(cltomerge2))=0;
                        Am1Matrix(cltomerge2,i3)=1;
                        Am1Matrix(i3,cltomerge2)=1;
                    end
                end

               %Delete old cluster 
               Am1Matrix(cltomerge1(cltomerge2),:) = [];
               Am1Matrix(:,cltomerge1(cltomerge2)) = [];
           end

           best_cl = cltomerge2;       
           clear cltomerge1 cltomerge2 pvalue dipval dipvalue maxpval cluster_eval2 pval2 gIdx2 pvalues pIndex1 pIndex2 PropMatrix2 propvalues maxprop prIndex1 prIndex2;

        end %end of an iteration of the algorithm
        
        for j2=1:k,  gIdx(clmemberIDs{j2}) = j2; end

        % calculate the overall error of the dataset
        if (k == 1)
            if (isempty(c)), c = ComputeCentroids(X, ones(n,1), 1, 0); end        
            Er(1) = sum(sqrt( sum(bsxfun(@minus, X, c).^2, 2) ));
        else
            for k1 =1:k
                Er(k1)= sum(sqrt( sum(bsxfun(@minus, X(clmemberIDs{k1},:), c(k1,:)).^2, 2) ));
            end
        end
        
        Er(find(Er==0)) = [];
        sumer = sum(Er);

        if(k>0)
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % REFINEMENT
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            if (refineMODE == 1)
                fprintf('** Refining final solution...\n');        
                % compute the centroids of the bisecting k-means algorithm
                c = ComputeCentroids(X, gIdx, k, 0);

                % run a refinement phase with a k-means version
                [gIdx_ref, c, sumd_ref, ~, ~] = ark_kmeans(X', [], c', 5, -1, 1, 0, 0, 0, 0);
                c = c';
            else
                gIdx_ref = gIdx; sumd_ref = sumer;
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
     
    %Find which cluster has moved where in the final clustering
    realClustFinal = [];
    
    for i1 = 1:size(clmemberIDs,2)
        for i2 = 1:size(clmemberIDs{:,i1},2)
            for i3 = 1:size(firstClusterMember,2)
                if(firstClusterMember(i3) == clmemberIDs{:,i1}(i2))
                    realClustFinal = [realClustFinal;i3 i1]; 
                end
            end
        end
    end 
end