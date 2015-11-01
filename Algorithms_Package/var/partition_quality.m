%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Input
%   R: the clustering assignment
%   C: ground truth labels
% Output
%   pq:          is a clustering quality measure related to Rand Index
%   conf_matrix: the confusion matrix of labels 
%                (rows: classes, cols: clusters of the result)
%   conf_matrix_probR: the confusion matrix with probabilities per cluster size
%   conf_matrix_probC: the confusion matrix with probabilities per class size
%------------
% Copyright (C) 2012-2013, Argyris Kalogeratos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pq, RI, ARI, conf_matrix, conf_matrix_probC, conf_matrix_probR] =  partition_quality(C, R)

    M = max(C);   % the number of classes in dataset
    k = max(R);     % the number of clusters in the result
    N = length(C);

    indC = find1DIndexMaps(C, 1:M);
    indR = find1DIndexMaps(R, 1:k);
    ndocs_C = sum(indC, 1)';
    ndocs_R = sum(indR, 1)';
    conf_matrix = double(indC') * double(indR);

    % sort clusters wrt to the class they represent better (most members come
    % from that data class)
    [~, best_class] = max(conf_matrix, [], 1);
    [~, ind] = sort(best_class, 'ascend');
    indR = indR(:,ind);
    ndocs_R = ndocs_R(ind);
    conf_matrix = conf_matrix(:, ind);


    conf_matrix_probC = bsxfun(@rdivide, conf_matrix, ndocs_C);
    conf_matrix_probR = bsxfun(@rdivide, conf_matrix, ndocs_R');

    p_i = ndocs_C / N;

    p_ij = conf_matrix / N;
    pq = sum(p_ij(:).^2) / sum(p_i.^2);

    [ARI,RI] = valid_RandIndex(C,R);
end