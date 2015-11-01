%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This script prints the final segments that agglodip
% algorithm has found on the image as well as the
% initial oversegmentation of the image for a easy
% visual comparison
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

for i=1:size(l2,1)
       for j=1:size(l2,2)
          for m1=1:size(realClustFinal,1)
               if realClustFinal(m1,1) == l2(i,j)
                   l2(i,j) = realClustFinal(m1,2);
               end
          end 
       end
end
 
str1 = sprintf('Number Of Clusters Initially : %d',C1s);
str2 = sprintf('Number Of Clusters Finally : %d',k);
show(drawregionboundaries(l1,im1, [255 0 144]),110,str1);
show(drawregionboundaries(l2,im1, [255 0 144]),111,str2);
clear str1 str2;