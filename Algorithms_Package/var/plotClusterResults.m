%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This script was created to plot the final clustering that
% an algorithm has found as well as the real clusters (if provided)
% so we can compare their plots if we wish. Note that
% this script plots only the first 2 dimensions (first
% along with second) of the data and that it has a 
% maximum limit of 56 clusters that it can plot with 
% different colors. Changes to aleviate the above 
% limitations are easy to be made.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Here new color/shape combinations can be added to plot the data
colorCell = {'r^','g^','b^','m^','c^','k^','y^','rs','gs','bs','ms','cs','ks','ys','rv','gv','bv','mv','cv','kv','yv','r+','g+','b+','m+','c+','k+','y+','ro','go','bo','mo','co','ko','yo','r.','g.','b.','m.','c.','k.','y.','rh','gh','bh','mh','ch','kh','yh','rp','gp','bp','mp','cp','kp','yp'};

% Final Clusters
figure(1);
hold on;
for clustk = 1:k_plot
    for leng = 1:length(R_ref_plot)
        if R_ref_plot(leng) == clustk

                plot(X(leng,1),X(leng,2),colorCell{clustk});
            
        end
    end
end
title('Final Clustering');

% Real Clusters provided by variable C
 if (real_k > 0)
    figure(2);
    hold on;
    for clustk = 1:real_k
        for leng = 1:length(C)
            if C(leng) == clustk

                    plot(X(leng,1),X(leng,2),colorCell{clustk});
            
            end
        end
    end
     title('Real Clusters');
 end


