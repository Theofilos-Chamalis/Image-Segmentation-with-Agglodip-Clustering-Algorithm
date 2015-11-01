%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This function was created to plot each step of the
% clustering process to see which cluster is merged
% in every iteration of the agglomerative algorithms.
% These clusters are pointed with yellow color.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [] = plotEachStep(X, gIdx_ref_Init, iterationNo, algorithmName, mergedCluster)

% plot the resulting clusters
k2 = length(unique(gIdx_ref_Init));
colorCell = {'r^','g^','b^','m^','c^','k^','rs','gs','bs','ms','cs','ks','rv','gv','bv','mv','cv','kv','r+','g+','b+','m+','c+','k+','ro','go','bo','mo','co','ko','r.','g.','b.','m.','c.','k.','rh','gh','bh','mh','ch','kh','rp','gp','bp','mp','cp','kp'};
colorCell2 = horzcat(colorCell,colorCell,colorCell,colorCell,colorCell);
figure(iterationNo+1);
hold on;
for clustk = 1:k2
    for leng = 1:length(gIdx_ref_Init)
        if gIdx_ref_Init(leng) == clustk
                if(mergedCluster == 0 || clustk~=mergedCluster)
                    plot(X(leng,1),X(leng,2),colorCell2{clustk});
                else
                    plot(X(leng,1),X(leng,2),'yo');
                end     
        end
    end
end

if(iterationNo == 0)
    iter = ' initial clusters';
elseif(iterationNo == 99999)
    iter = ' final clustering';
else
    iter = strcat(' iteration No ',num2str(iterationNo));
end

s = strcat(algorithmName,iter);
title(s);
% pause to see the results
pause;

end