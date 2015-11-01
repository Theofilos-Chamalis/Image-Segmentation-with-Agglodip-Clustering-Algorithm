%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Unimodality tests are done on each cluster produced by the global kmeans 
% or the classic kmeans algorithm to check if we need to split them based
% on if they are unimodal or not. The algorithm that decides for each split
% is pdip-means.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gIdx_ref_Init, c_all, k2, kmax] = UnimodalityChecking(X, gIdx_ref_Init, initclNo, smallest_cluster, mergeSELECT)
    check_timer = tic;
    
    fprintf('=====================Start of unimodality check=====================\n\n');
	% define the RNG seed
	%rseed = sum(100*clock);    
	rseed = 10;
	% rand('state', rseed);  
	% randn('state', rseed);
    split_struct = cell(5,1);
    split_struct{2} = struct;
    split_struct{2}.pval_threshold    = 0.00; 
    split_struct{2}.exhaustive_search = 1; 
    split_struct{2}.voting            = 0.01;
    split_struct{2}.nboot             = 1000;
    split_struct{2}.overall_distr     = 0;

    Y = zeros(1,length(X(1)));
    Y(1) = 99999;
    y = 1;

    while(true)
        for t=1:length(X)
            if(gIdx_ref_Init(t)==y)
                if(Y(1,1)==99999)
                    Y = X(t,:);
                end
                Y = [Y;X(t,:)];
            end        
        end

        % use pdip-means algorithm for the splitting
        [R_temp, sumer_temp, R_ref_temp, sumer_ref_temp] = bisect_kmeans_default(Y, 'split_struct', split_struct{2}, 'split_trials', 10, 'splitSELECT', 6, 'splitMODE', 0, 'refineMODE', 1, 'smallest_cluster', smallest_cluster, 'attempts', 1, 'rndseed', 0+rseed);
        k_ref_temp = length(unique(R_ref_temp)); 

        if(k_ref_temp>1)
            for t=1:length(X)
                if(gIdx_ref_Init(t)>y)
                    gIdx_ref_Init(t) = gIdx_ref_Init(t) + k_ref_temp - 1;
                end
            end

            for t=1:length(X)
                if(gIdx_ref_Init(t)==y)
                    gIdx_ref_Init(t) = gIdx_ref_Init(t) + R_ref_temp(1) - 1;
                    R_ref_temp(1) = [];
                end
            end

            y = y + k_ref_temp - 1;
            initclNo = initclNo + k_ref_temp - 1;
        end

        Y = zeros(1,length(X(1)));
        Y(1) = 99999;
        if(y>=initclNo)
            break;
        end
        y = y+1;
        clear R_temp sumer_temp R_ref_temp sumer_ref_temp k_ref_temp;
    end

    k2 = length(unique(gIdx_ref_Init));
    if (mergeSELECT > 2), kmax = k2; end
    
    c_all = ComputeCentroids(X, gIdx_ref_Init, k2, 0);
    fprintf('\nAfter the unimodality check the number of clusters is: %d',k2);
    if(initclNo == k2)
        fprintf(' and no split occured!');
    end
    
    check_timer_stop = toc(check_timer);
    fprintf('\nThe elapsed time for the unimodality check is %f seconds\n\n',check_timer_stop);
    fprintf('=====================End of unimodality check=====================\n\n\n');

end