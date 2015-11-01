%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% k-means clustering
%
% [Idx, c, E, D, rnd_num] = kmeans(X, k, init_type, iters, attempts_num, store_best_infile, save_res, keep_best_D, exp_id)
%
% Input:
% X             - (d x n) with column vectors
% W             - weigths for each vector instance
% k             -  number of clusters
% init_type     - various initializations
% iters         - maximum number of iterations
% attempts_num  - number of restarts
% store_best_infile - store best solution in a output file 
% save_res          - save results in files (true/false)
% keep_best_D       - keep the best distance matrix (true/false)
% exp_id            - a unique identifier for the output files


% Output:
% Idx   - (n x 1) index of nearest center for each data point
% c     - (k x d) cluster centers in row vectors
% E     - sum of squared distances to nearest centroid
% D     - the distance matrix
% rnd_num - a random identifier
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
% Aknowledgement:
% the anchors() function is taken from fastkmeans package available from
% Matlab File Exchange. Lightspeed package is also required
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Idx, c, E, D, rnd_num, Ei] = ark_kmeans(X, W, k, init_type, iters, attempts_num, store_best_infile, save_res, keep_best_D, exp_id)

    best_E = flintmax;
    [d, n] = size(X);

    if (n < k), fprintf('Error: more clusters(%g) than data instances(%g)!\n', k, n); end;

    % consider uniform unit weights
    if (isempty(W) || sum(W) == n),   %W = ones(1, n);
                                     WEIGHTED_KMEANS = 0;
    else                             WEIGHTED_KMEANS = 1;
    end

    %in a different case for unormalized data:   
    SqSum_X = sum(X.^2, 1);
    % SqSum_X = ones(1,n);
    SqSum_Xall = sum(X(:).^2); % this is x^2 for all vector instances

    rnd_num = randi([1 1000000],1,1); % this is an id for the file output
    if (~exist('exp_id', 'var'));
         rand ('state', sum(100*clock));
         fid_str = ['best_res' num2str(rnd_num) '.mat'];
    else rand ('state', 0 + exp_id);
         fid_str = ['best_res' num2str(exp_id) '_' num2str(rnd_num) '.mat'];
    end    

    DATA_PREC = 'double';

    % set no limitation to the maximum iterations of the k-sp procedure 
    if (iters <= 0),  
        iters = 10000; %intmax('int32');  
    end
    % if iters in [0, 1) then the algorithm stops when less than MAXITER objects reassigned

    if (issparse(X)), SPARSE_X = 1;
    else              SPARSE_X = 0;
    end

    for attempt=1:attempts_num,
        done    = 0;
        iter    = 0;
        E       = zeros(iters, 1);
        Ei      = zeros(1, n);
        Ei_new  = zeros(1, n);
        Idx_old = zeros(1, n);
        Idx     = zeros(n, 1);
        zerocls = 0;

        if (numel(k) > 1);     % initialize with given means, compute assignment only
            c = full(k);
            k = size(c,2);
        else
            if (init_type == 1), 
                r = randperm(n); r = r(1:k);
                c = full( X(:, r) );       
            elseif (init_type == 2),  % k-means++
                % for the first centroid
                c = zeros(d, k,  DATA_PREC);
                r = randi([1 n],1,1);
                c(:, 1) = full( X(:, r) );   % the first centroid is uniform randomly selected
                valid_ids = ([1:r-1, r+1:n]);

                D(valid_ids, 1) = sum( sqrt( bsxfun(@minus, X(:,valid_ids), c(:, 1)).^2 ) );

                % for the rest k-1 centroids
                for t=2:k,
                    %dbg: disp (t)
                    % only for the new cluster center
                    D(valid_ids, t) = sqrt( sum(bsxfun(@minus, X(:,valid_ids), c(:, t)).^2,1) ); %sum((X - repmat(c(:,1), [n 1])).^2, 2);

                    % Partition data to closest centroids
                    v = min(D(valid_ids, 1:t), [], 2);
                    v = cumsum(v / sum(v));

                    cid = find(v > rand, 1);
                    c(:,t) = full( X(:, valid_ids(cid)) );
                    valid_ids(cid) = [];
                end
            elseif (init_type == 3), % furthest-first method, mean vector as first center
                c = anchors(sum(X,2)', k, X');
                c = full( c' );
            elseif (init_type == 4), % furthest-first method, random first center from data
                c = anchors(X(:,randi([1 n],1,1))', k, X');
                c = full( c' );
            end

            clchanged = true(k,1);
            changed_clusters = 1:k;
        end 

        D = zeros(k,n);
        clchanged = true(k,1);
        changed_clusters = 1:k;
        c_sum = zeros(d, k, DATA_PREC);
        cluster_ids = 1:k;
        clmemberIDs = cell(k,1);

        clmembers = zeros(1,k);
        sum_cl_weights = zeros(1,k);

    %%%%%%%%%%%%%%%%%%%%%%%%%
    % loop until converge
    %%%%%%%%%%%%%%%%%%%%%%%%%
        while ~done
            iter = iter + 1;

            SqSum_c = sum(c.^2, 1) * 0.5;  % this is c^2 for each centroid (divided by 2 in order to multiply it in the next step)

            % only for the changed clusters
                partial_computations = 1;
                %partial_computations
                if (partial_computations > 1) 
                    % if n >> number of clusters
                        j = 0; csz = ceil(n/partial_computations);
                        for i=1:partial_computations,
                            obj_id = 1+j : min(j+csz, n);
                            D(changed_clusters, obj_id) = bsxfun(@minus, c(:,changed_clusters)'*X(:,obj_id), SqSum_c(changed_clusters)');
                            j = j + length(obj_id);
                        end
                else D(changed_clusters, :) = bsxfun(@minus, c(:,changed_clusters)'*X, SqSum_c(changed_clusters)');
                end


                [Ei Idx] = max(D, [], 1);              
                Idx = Idx';

            %find reassigned objects
            if (iter > 1)
                reassigned = find(Idx ~= Idx_old); % note the objects that reassigned to different clusters
                n_reassigned = length(reassigned); 

                clchanged(:) = false;

                tmpmat = false(k, 1); 
                tmpmat(Idx(reassigned)) = true;  tmpmat(Idx_old(reassigned)) = true;

                clchanged(tmpmat) = true;
                changed_clusters  = cluster_ids(clchanged);
                k_changed = length(changed_clusters);
            else
                reassigned = 1:n;
                n_reassigned = n;
                k_changed = k;
                changed_clusters = 1:k;
            end

            E(iter) = sum( SqSum_X - 2*Ei ); % avoid sqrt computation

            if ((iters >=1 && iter >= iters) || (iters >  1 && n_reassigned == 0) || (iters < 1 && n_reassigned > iters*DATASIZE))  % terminate after given # of iterations
                done = 1; zerocls = 0;
            else
                old_members = clmemberIDs(changed_clusters);
                clmemberIDs(changed_clusters) = find1DIndices(Idx(:), changed_clusters);            
                clmembers(changed_clusters) = cellfun(@length, clmemberIDs(changed_clusters));

                k_to_check = 1:k_changed;
                zerocls = sum(clmembers(changed_clusters) == 0);
                if (zerocls > 0), k_to_check = k_to_check( clmembers(changed_clusters(k_to_check)) > 0 ); end
                c_sum_update = zeros(d, length(k_to_check));
                q_ind = 1;

                if (WEIGHTED_KMEANS == 0)
                    sum_cl_weights(changed_clusters) = clmembers(changed_clusters); % clmembers is mutable and sum_cl_weights is its reference

                    for q=k_to_check,
                        Idx_q = changed_clusters(q);
                        left = old_members{q}(Idx(old_members{q}) ~= Idx_q);   
                        came = clmemberIDs{Idx_q}(Idx_old(clmemberIDs{Idx_q}) ~= Idx_q);
                        c_sum_update(:,q_ind) = full(row_sum(X(:,came))) - full(row_sum(X(:,left)));
                        q_ind = q_ind + 1;
                    end
                else 
                    % compute the sum of weights of objects in each cluster
                    sum_cl_weights(changed_clusters) = cellfun(@(x) sum(W(x)), clmemberIDs(changed_clusters), 'UniformOutput', true);

                    for q=k_to_check,
                       Idx_q = changed_clusters(q);

                       left = old_members{q}(Idx(old_members{q}) ~= Idx_q);   
                       came = clmemberIDs{Idx_q}(Idx_old(clmemberIDs{Idx_q}) ~= Idx_q);
                       if (~isempty(came)),  c_sum_update(:,q_ind) = full(row_sum(bsxfun(@times, X(:,came), W(came)))); end
                       if (~isempty(left)),  c_sum_update(:,q_ind) = c_sum_update(:,q_ind) - full(row_sum(bsxfun(@times, X(:,left), W(left))));  end                       
                       q_ind = q_ind + 1;
                    end

                end

                c_sum(:,changed_clusters(k_to_check)) = c_sum(:,changed_clusters(k_to_check)) + c_sum_update;
            end

                % handle empty clusters
                    if (zerocls > 0)
                        clzeromembers = find(clmembers == 0, zerocls);
                        % find least represented objects, without making empty any cluster
                        [~, new_objects] = sort(SqSum_X - 2*Ei, 'descend'); % SqSum_X will be all the same for normalized data

                         q = 1; q_ind = 1;                   
                        while (q <= zerocls),
                            obj = new_objects(q_ind);
                            Idx_q = Idx(obj);                        
                            if (clmembers(Idx_q) == 1), q_ind = q_ind + 1;
                                                        continue;
                            end
                            new_objects(q) = obj;
                            clmembers(Idx_q) = clmembers(Idx_q) - 1; % this should stay here
                            if (WEIGHTED_KMEANS == 0), sum_cl_weights(Idx_q) = clmembers(Idx_q);
                            else                       sum_cl_weights(Idx_q) = sum_cl_weights(Idx_q) - W(obj);
                            end
                            clmemberIDs{Idx_q}(clmemberIDs{Idx_q} == obj) = [];
                            q = q + 1;  q_ind = q_ind + 1;
                        end
                        new_objects = new_objects(1:zerocls);

                        clmembers(clzeromembers) = 1;

                        if (WEIGHTED_KMEANS == 0)
                             c_sum(:,clzeromembers)        = full( X(:, new_objects) );
                             c_sum(:,Idx(new_objects))     = c_sum(:,Idx(new_objects)) - full( X(:, new_objects) );
                             sum_cl_weights(clzeromembers) = clmembers(clzeromembers);
                        else new_vector = full( bsxfun(@times, X(:, new_objects), W(new_objects)) );
                             c_sum(:,clzeromembers)        = new_vector;
                             c_sum(:,Idx(new_objects))     = c_sum(:,Idx(new_objects)) - new_vector;
                             sum_cl_weights(clzeromembers) = W(new_objects);
                        end
                        clchanged(Idx(new_objects)) = true;  
                        clchanged(clzeromembers) = true;
                        changed_clusters = find(clchanged == true);

                        Idx(new_objects) = clzeromembers;
                        clmemberIDs(clzeromembers) = num2cell(new_objects');

                        % update clustering error 
                        E_update = sum( Ei(new_objects) ); % the previous error that must be substracted
                        Ei(new_objects) = 0;

                        E_update = sum(Ei(new_objects)) - E_update;
                        E(iter) = E(iter) + E_update;  % avoid sqrt
                        clear('new_vector');
                    end               


                % IF WEIGHTS ARE USED THEN THE NORMALIZATION IS THE SUM OF CLUSTER WEIGHTS
                % find mean centroid vector
                if (done ~= 1)
                    c(:,changed_clusters) = bsxfun(@times, c_sum(:,changed_clusters), 1 ./ sum_cl_weights(changed_clusters));
                end

            Idx_old = Idx;
        end

        Ei = sqrt(SqSum_X - 2*Ei);
        if (WEIGHTED_KMEANS == 1),   E = sum( Ei .* W ); % compute the sqrt on the error function
        else                         E = sum( Ei      ); % compute the sqrt on the error function
        end

        if (E < best_E)
            best_E = E;        

            if (store_best_infile)
                if (keep_best_D), save(fid_str, 'E', 'Idx', 'c', 'D', '-v7.3');
                else              save(fid_str, 'E', 'Idx', 'c');
                end
            else       
                best_Idx = Idx;
                best_c   = c;
                if (keep_best_D), best_D = D; end
            end
        end
    end

    % Manage the output info
    if (store_best_infile) % the file has already been written on disk
        if (attempts_num > 1 && best_E ~= E),
            if (keep_best_D == true), load(fid_str, 'E', 'Idx', 'c', 'D');                      
            else                      load(fid_str, 'E', 'Idx', 'c');    
            end        
        %else
        end
    else                
        if (attempts_num > 1 && best_E ~= E),
            E   = best_E; 
            Idx = best_Idx;
            c   = best_c;
            clear('best_E', 'best_Idx', 'best_c');
        end

        if (save_res)
            if (keep_best_D),
                 D = D_best;
            %else %save(fid_str, 'E', 'Idx', 'c');
            end
        end
    end 
end




