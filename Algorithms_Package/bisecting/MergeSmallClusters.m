%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Merges the small clusters to the clusters with the nearest centroid.
%------------
% Input
%   X:                    the dataset (row vectors)
%   gIdx_ref_Init:     the matrix that contains the object-to-cluster assingment
%   c_all:                the centroids set
%   smallest_cluster: a threshold that defines which clusters are considered small
% Output
%   gIdx_ref_Init:     the new matrix that contains the object-to-cluster assingment
%   c_all:                the new centroids set
%   initclNo:            the number of clusters
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [gIdx_ref_Init,c_all,initclNo,Am1Mat] = MergeSmallClusters(gIdx_ref_Init,c_all,smallest_cluster,X, Am1Mat, LOCAL_SEARCH_FLAG)

    repz = 1;
    az = unique(gIdx_ref_Init);
    out = [az,histc(gIdx_ref_Init(:),az)];
    smallcl=[0];
    rept = 1;

    for i1=1:length(out)
        if(out(i1,2)<smallest_cluster)
            if(smallcl(1)==0)
                smallcl(1)=i1;
            else
                smallcl=[smallcl;i1];
            end        
        end
    end
    
    while(smallcl(1,1)~=0)
        % Check if there are any small clusters 
        if(smallcl(1,1)~=0)
            smallcl(:,2)=0;
            smallclsize = size(smallcl);

            for i1=1:smallclsize(1,1)
                    q = pdist2(c_all(smallcl(i1,1),:),c_all);
                    q(smallcl(i1)) = 999;
                    [~,smalldist] = min(q);
                    smallcl(i1,2) = smalldist;
                    smallcl2(i1,1) = smallcl(i1,2);
                    smallcl2(i1,2) = smallcl(i1,1);
            end

            % Remove the duplicates from smallcl
             [La,Lb]=ismember(smallcl,smallcl2,'rows');
             samepairs = [];

             for i1=1:length(Lb)
                 if(Lb(i1)~=0)
                     if(Lb(i1)<=Lb(Lb(i1)))
                         samepairs = [samepairs; Lb(i1),Lb(Lb(i1))];
                     end
                 end
             end

            if(~isempty(samepairs))
               samepairs=sort(samepairs,1,'descend');
               samepairs(:,1) = samepairs(:,2);
               samepairs(:,2) = [];
               for i1=1:length(samepairs)
                smallcl(samepairs(i1),:)=[];
               end
               smallcl = flipud(smallcl);
            else
                smallcl = flipud(smallcl);
            end

            smallclsize = [];
            smallclsize = size(smallcl);      
            
            for i1=1:smallclsize(1,1)
                for i2=1:length(gIdx_ref_Init)
                    if(gIdx_ref_Init(i2) == smallcl(i1,1))
                        gIdx_ref_Init(i2) = smallcl(i1,2);
                        c_all(smallcl(i1,1)) = c_all(smallcl(i1,2));
                    end
                end
                
                for i3=1:size(Am1Mat,2)
                    if(Am1Mat(smallcl(i1,1),i3)==1)
                        Am1Mat(smallcl(i1,1),i3)=0;
                        Am1Mat(i3,smallcl(i1,1))=0;
                        Am1Mat(smallcl(i1,2),i3)=1;
                        Am1Mat(i3,smallcl(i1,2))=1;
                    end
                end
            end
            
			clear az;
            az = unique(gIdx_ref_Init);
            out2 = [az,histc(gIdx_ref_Init(:),az)];
            out2 = flipud(out2);
            
            %Eliminate empty clusters and reduce their number
            for i1 = 1:size(out2,1)
                if(out2(i1,2)==0)
                   for i2 = 1:length(gIdx_ref_Init)
                      if gIdx_ref_Init(i2) > out2(i1,1)
                          gIdx_ref_Init(i2) = gIdx_ref_Init(i2) - 1;
                      end
                   end
                   Am1Mat(out2(i1,1),:) = [];
                   Am1Mat(:,out2(i1,1)) = [];
                   c_all(out2(i1,1),:)=[];
                end
            end
        end
                
        if(LOCAL_SEARCH_FLAG == 0)
%             [gIdx_ref_Init, c_all, ~, ~, ~] = ark_kmeans(X', [], c_all', 20, -1, 1, 0, 0, 0, 0);
%             c_all = c_all';
        end

        
        %Fill again smallcl variable to check for small clusters
        clear az out out2 smallcl;

        az = unique(gIdx_ref_Init);
        out = [az,histc(gIdx_ref_Init(:),az)];
        smallcl=[0];

        for i1=1:length(out)
            if(out(i1,2)<smallest_cluster)
                if(smallcl(1)==0)
                    smallcl(1)=i1;
                else
                    smallcl=[smallcl;i1];
                end        
            end
        end
        
        initclNo = length(unique(gIdx_ref_Init));
        repz = repz+1;
    end
	
	initclNo = length(unique(gIdx_ref_Init));
	fprintf('\n\nAfter merging the small clusters, the number of clusters is now: %d\n\n',initclNo);
	
end
