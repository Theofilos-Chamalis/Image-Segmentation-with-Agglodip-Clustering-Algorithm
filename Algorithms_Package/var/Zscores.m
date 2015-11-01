%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Function used to transform the sampling distribution into a standard normal distribution
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Z = Zscores(X)
    
    means_col = mean(X,1);
    stdev_col = std(X,1,1);
    [m,n] = size(X);
    Z = zeros(m,n);
      
    for i = 1:m
        for j = 1:n
            Z(i,j) = (X(i,j)-means_col(j))/stdev_col(j);
        end        
    end
    
end