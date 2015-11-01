%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% This functions converts the image im1 to grayscale and
% for every pixel we represent it by its neighbors within
% a window of size windowSize. Then we apply PCA on the
% set of these pixels to reduce the dimensions and output
% the matrix Xgray_PCA_Init2 that contains them.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Xgray_PCA_Init2 = PCA_Mat_Create(im1,windowSize)

    im4 = rgb2gray(im1);
    im4length1 = size(im4,1);
    im4length2 = size(im4,2);
    
    % A window size of 2 means that we use pixels within of distance 2 from the current one
    % to describe it, i.e. 24 pixels around the current (center) one and itself.
    winSize = windowSize;
    bRows = im4length1*im4length2;
    bCols = (2*winSize+1)^2;
    Xgray_PCA_Init = ones(bRows,bCols)*9999;

    %We are inside the accepted pixel range of the image so we fill
    %Xgray_PCA_Init cells with the appropriate neighboring values
    for l1 = 1:im4length1
        for l2 = 1:im4length2
            if(l1>winSize && l1<=im4length1-winSize && l2>winSize && l2<=im4length2-winSize)
                  row2BeInserted = [];
                  for ii2 = l1-winSize:l1+winSize
                      for ii3 = l2-winSize:l2+winSize
                          row2BeInserted = [row2BeInserted im4(ii2,ii3)];
                      end
                  end
                  Xgray_PCA_Init(((l1-1)*im4length2)+l2 , :) = row2BeInserted; 
            end

        end
    end


    %Initialize and crop the Xgray_PCA_Init2 matrix to correct dimensions
    dimensionsToKeep = windowSize*2+1;
    Xgray_PCA_Init2 = Xgray_PCA_Init;
    Xgray_PCA_Init2(:,dimensionsToKeep+1:size(Xgray_PCA_Init,2)) = [];
    
    %Compute the PCA of the whole Xgray_PCA_Init2 matrix
    [~,score,~] = pca(Xgray_PCA_Init); 
    Xgray_PCA_Init2 = score(:,1:dimensionsToKeep);
end