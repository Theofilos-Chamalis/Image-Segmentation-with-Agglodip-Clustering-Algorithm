%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% The preprocessing phase of the segmentation is executed here.
% First an image is given as input and the SLIC (Simple Linear
% Iterative Algorithm) algorithm produces an oversegmentation of
% the image, with the number of superpixels (groups of pixels) we
% instruct it to. Then there are 2 approaches; the first converts
% the pixels from RGB color space to the CIE 1976 (L*,A*,B*)
% color space and the second one transforms the image to grayscale
% and then every pixel is represented by itself and its neighbors that
% lie within a window of size w. PCA is applied after that to the whole
% set to reduce the dimensions. After that we choose the central 
% pixels of each superpixel to represent each superpixel/cluster and
% give them as input to agglodip hierarchical clustering algorithm to
% execute the segmentation by merging neighbor superpixels.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
%PCA approach (enter 1 to use or 0 to use (L*,A*,B*) approach)
PCA_FLAG = 0;

%Input of the image to be used for segmentation
im1 = imread('a481x321(110).jpg');
im2 = rgb2lab(im1); %convert RGB to LAB
im3 = im1;
im4 = im2;


if PCA_FLAG == 1
    %Input params are image and window size (default = 2)
    Xpca = PCA_Mat_Create(im1,2);
end

%execute the SLIC algorithm that produces the superpixels of the image
[l1, Am1, C1] = slic(im1,100,30,2,'median'); 
C1s = size(C1,2);

%Use of the 'central' pixels of the superpixel and ignore the rest
%A central window size of 5 means that out of each row and column in the superpixel, 
%we get the 5 central pixels and then we use the intersection of all these to use as 
%data of each cluster for input in agglomerative dip means. The lowest value 4
CENTRAL_WINDOW_FLAG = 1; 
centralWindowSize = 5; 

ltemp = zeros(1,size(l1,2));
ltemp1 = [];

%Length and width of the whole image
sizel11 = size(l1,1);
sizel12 = size(l1,2);

if CENTRAL_WINDOW_FLAG == 1
    
    im3With2550Values = zeros(size(im2,1),size(im2,2),size(im2,3));
    im3With2550Values(:,:,:) = 2550; % create a matrix with very big values as initialization so then we can distinguish these values from normal ones of the image
    rowIndicesOfPixelsToKeep{1,1} = {};
    colIndicesOfPixelsToKeep{1,1} = {};

    for r1 = 1:sizel11
        uniqueElementsOnRows{r1} = unique(l1(r1,:),'stable');
    end
    uniqueElementsOnRows = uniqueElementsOnRows';
    
    for s1 = 1:sizel12
        uniqueElementsOnCols{s1} = unique(l1(:,s1),'stable');
    end
    
    %Find the indices of pixels to keep on rows and columns respectively
    for r2 = 1:size(uniqueElementsOnRows,1)
        for r3 = 1:size(uniqueElementsOnRows{r2,1},2)
            frequencyOfElementsOnRows{r2,1}(r3) = numel(l1(l1(r2,:) == uniqueElementsOnRows{r2,1}(r3)));

            if(r3>1)
                cumulativeFrequencyOfElementsOnRows{r2,1}(r3) = cumulativeFrequencyOfElementsOnRows{r2,1}(r3-1) + frequencyOfElementsOnRows{r2,1}(r3);
            else
                cumulativeFrequencyOfElementsOnRows{r2,1}(r3) = frequencyOfElementsOnRows{r2,1}(r3);
            end

            halfNumber = floor((frequencyOfElementsOnRows{r2,1}(r3))/2);
            halfBound = floor(centralWindowSize/2);

            if(halfNumber>centralWindowSize+2)
                if(r3>1)
                    for r4 = - halfBound: halfBound
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (cumulativeFrequencyOfElementsOnRows{r2,1}(r3-1) + halfNumber + r4)];               
                    end
                else
                    rowIndicesOfPixelsToKeep{r2,1} = {};
                    for r4 = - halfBound: halfBound
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (halfNumber + r4)];
                    end
                end
            elseif(halfNumber<=centralWindowSize+2 && halfNumber>3)
                if(r3>1)
                    for r4 = - halfBound+1: halfBound-1
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (cumulativeFrequencyOfElementsOnRows{r2,1}(r3-1) + halfNumber + r4)];               
                    end
                else
                    rowIndicesOfPixelsToKeep{r2,1} = {};
                    for r4 = - halfBound+1: halfBound-1
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (halfNumber + r4)];
                    end
                end
            elseif(halfNumber == 3)
                if(r3>1)
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (cumulativeFrequencyOfElementsOnRows{r2,1}(r3-1) + halfNumber)];               
                else
                        rowIndicesOfPixelsToKeep{r2,1} = {};
                        rowIndicesOfPixelsToKeep{r2,1} = [rowIndicesOfPixelsToKeep{r2,1}, (halfNumber)];
                end
            else
                if(r3<=1)
                    rowIndicesOfPixelsToKeep{r2,1} = {};
                end
            end
        end
    end 
    
    clear halfNumber halfBound;
    
     for s2 = 1:size(uniqueElementsOnCols,2)
        for s3 = 1:size(uniqueElementsOnCols{1,s2},1)
            frequencyOfElementsOnCols{1,s2}(s3) = numel(l1(l1(:,s2) == uniqueElementsOnCols{1,s2}(s3)));

            if(s3>1)
                cumulativeFrequencyOfElementsOnCols{1,s2}(s3) = cumulativeFrequencyOfElementsOnCols{1,s2}(s3-1) + frequencyOfElementsOnCols{1,s2}(s3);
            else
                cumulativeFrequencyOfElementsOnCols{1,s2}(s3) = frequencyOfElementsOnCols{1,s2}(s3);
            end

            halfNumber = floor((frequencyOfElementsOnCols{1,s2}(s3))/2);
            halfBound = floor(centralWindowSize/2);

            if(halfNumber>centralWindowSize+2)
                if(s3>1)
                    for s4 = - halfBound: halfBound
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (cumulativeFrequencyOfElementsOnCols{1,s2}(s3-1) + halfNumber + s4)];               
                    end
                else
                    colIndicesOfPixelsToKeep{1,s2} = {};
                    for s4 = - halfBound: halfBound  
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (halfNumber + s4)];
                    end
                end
            elseif(halfNumber<=centralWindowSize+2 && halfNumber>3)
                if(s3>1)
                    for s4 = - halfBound+1: halfBound-1
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (cumulativeFrequencyOfElementsOnCols{1,s2}(s3-1) + halfNumber + s4)];               
                    end
                else
                    colIndicesOfPixelsToKeep{1,s2} = {};
                    for s4 = - halfBound+1: halfBound-1
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (halfNumber + s4)];
                    end
                end
            elseif(halfNumber == 3)
                if(s3>1)
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (cumulativeFrequencyOfElementsOnCols{1,s2}(s3-1) + halfNumber )];                   
                else
                        colIndicesOfPixelsToKeep{1,s2} = {};
                        colIndicesOfPixelsToKeep{1,s2} = [colIndicesOfPixelsToKeep{1,s2}, (halfNumber)];
                end
            else
                if(s3<=1)
                    colIndicesOfPixelsToKeep{1,s2} = {};
                end
            end
        end
    end
    
    rowIndicesOfPixelsToKeep2 = {};
    rowIndicesOfPixelsToKeep3 = {};

    for a1=1:size(rowIndicesOfPixelsToKeep,1)
        rowIndicesOfPixelsToKeep2{a1,1} = {};
        rowIndicesOfPixelsToKeep3{a1,1} = {};
    end

    for w1 = 1:size(colIndicesOfPixelsToKeep,2)
        for w2 = 1:size(colIndicesOfPixelsToKeep{1,w1},2)
            g1 = cell2mat(colIndicesOfPixelsToKeep{1,w1}(w2));
            rowIndicesOfPixelsToKeep2{g1,1} = [ rowIndicesOfPixelsToKeep2{g1,1},w1];
        end
    end

    for w3 = 1:size(rowIndicesOfPixelsToKeep,1)
        if(isempty(rowIndicesOfPixelsToKeep2{w3,1}) || isempty(intersect(cell2mat(rowIndicesOfPixelsToKeep{w3,1}),cell2mat(rowIndicesOfPixelsToKeep2{w3,1}))))
            rowIndicesOfPixelsToKeep3{w3,1} = {0};
        else
            rowIndicesOfPixelsToKeep3{w3,1} = [rowIndicesOfPixelsToKeep3{w3,1}, intersect(cell2mat(rowIndicesOfPixelsToKeep{w3,1}),cell2mat(rowIndicesOfPixelsToKeep2{w3,1}))];
        end
    end
    
end

%Write on the im3with2550values matrix the LAB values of the pixels to keep,
%and reformat the ltemps so that they can be used as the gIdx matrix later
%(assignments of pixels to clusters)
for y1=1:sizel11
    for y2=1:sizel12
        ltemp(y2) = l1(y1,y2);        
        if(y2 == sizel12)
            ltemp2 = ltemp';
            ltemp1 = [ltemp1;ltemp2];
            ltemp = zeros(1,sizel12);
        end
              
        if CENTRAL_WINDOW_FLAG == 1
            if(cell2mat(rowIndicesOfPixelsToKeep3{y1,1}) ~= 0)
                for u1=1:size(rowIndicesOfPixelsToKeep3{y1,1}{1},2)
                    if(rowIndicesOfPixelsToKeep3{y1,1}{1}(u1) == y2)
                        im3With2550Values(y1,y2,1) = im2(y1,y2,1);
                        im3(y1,y2,1) = 110;
                        im3With2550Values(y1,y2,2) = im2(y1,y2,2);
                        im3(y1,y2,2) = 0;
                        im3With2550Values(y1,y2,3) = im2(y1,y2,3);
                        im3(y1,y2,3) = 255;
                    end
                end                
            end
        end
        
    end
end

im2temp = zeros(size(im2,2),3);
im2temp1 = [];
sizeim21 = size(im2,1);
sizeim22 = size(im2,2);

for y1=1:sizeim21
    for y2=1:sizeim22
        if CENTRAL_WINDOW_FLAG == 0
            im2temp(y2,1) = im2(y1,y2,1);
            im2temp(y2,2) = im2(y1,y2,2);
            im2temp(y2,3) = im2(y1,y2,3);
            if(y2 == sizeim22)
                im2temp1 = [im2temp1;im2temp];
                im2temp = zeros(sizeim22,3);
            end
        else
            im2temp(y2,1) = im3With2550Values(y1,y2,1);
            im2(y1,y2,1) =  im3With2550Values(y1,y2,1);
            
            im2temp(y2,2) = im3With2550Values(y1,y2,2);
            im2(y1,y2,2) =  im3With2550Values(y1,y2,2);
            
            im2temp(y2,3) = im3With2550Values(y1,y2,3);
            im2(y1,y2,3) =  im3With2550Values(y1,y2,3);
            if(y2 == sizeim22)
                im2temp1 = [im2temp1;im2temp];
                im2temp = zeros(sizeim22,3);
            end
        end
    end
end

%Use the central pixels of the generated superpixels
if CENTRAL_WINDOW_FLAG == 1
    removedNo = 0;
    sizeim2temp11 = size(im2temp1,1);
    
    for h=1:sizeim2temp11
        if(im2temp1(h,1)==2550)
            removedNo = removedNo+1;
        end
    end
    
    finalSize = sizeim2temp11 - removedNo;
    removedIndices = zeros(removedNo,1);
    o1 = 1;
    
    for h=1:sizeim2temp11
        if(im2temp1(h,1)==2550)
            removedIndices(o1) = h;
            o1 = o1 + 1;
        end
    end
  
    allIndeces = zeros(sizeim2temp11,1);
    
    for v1 = 1:sizeim2temp11
        allIndeces(v1,1) = v1;
    end
    
    diffIndices = setdiff(allIndeces,removedIndices);
    
    im2temp2 = zeros(finalSize,3);
    ltemp3 = zeros(finalSize,1);
    
    for q1 = 1:finalSize
        im2temp2(q1,1) = im2temp1(diffIndices(q1),1);
        im2temp2(q1,2) = im2temp1(diffIndices(q1),2);
        im2temp2(q1,3) = im2temp1(diffIndices(q1),3);
        ltemp3(q1,1) = ltemp1(diffIndices(q1),1);
    end

    %Final variable assignments to be used by the bistest script
    gIdx_ref_Init = ltemp3;
    Xfinal = im2temp2;
    %Uncomment this to use normalized values
    % X = normc(Xfinal); 
    X = Xfinal;
else %Use the whole dataset as input to agglodip if central window flag is 0
    gIdx_ref_Init = ltemp1;
    Xfinal = im2temp1;
    %Uncomment this to use normalized values
    % X = normc(Xfinal);
    X = Xfinal;
end

if PCA_FLAG == 1
    finalSize2 = size(Xpca,2);
    Xpca_grayscale = zeros(finalSize,finalSize2);
    %Keep only the central pixels from the grayscale image
    for q2 = 1:finalSize
        for q3 = 1:finalSize2
           Xpca_grayscale(q2,q3) = Xpca(diffIndices(q2),q3); 
        end
    end
    
    %Prepare the output variable X
    clear c_all X;
    X = Xpca_grayscale;

end


kmax = C1s + 1;
numberOfInitialClustersEXT = C1s;
l2 = l1;

%Look for empty clusters
valuez = unique(gIdx_ref_Init);
vlzLength = length(valuez);
mxVlz = C1s;
indxMatt = ones(mxVlz,1);

for i=1:mxVlz
   indxMatt(i) = i; 
end

emptyClusters = setdiff(indxMatt,valuez);

%Look for elements of the empty clusters
AdditionalElementsIndx = [];
for i = 1:length(emptyClusters)
    [rowz,colz] = find(l1 == emptyClusters(i));
    Lrowz = length(rowz);

    if Lrowz>=16
        for i = -floor(Lrowz/4):1:floor(Lrowz/4)
            AdditionalElementsIndx = [AdditionalElementsIndx; rowz(floor(Lrowz/2)+i) colz(floor(Lrowz/2)+i)];
        end
    elseif Lrowz>=9
        for i = -floor(Lrowz/3)+1:1:floor(Lrowz/3)-1
            AdditionalElementsIndx = [AdditionalElementsIndx; rowz(floor(Lrowz/2)+i) colz(floor(Lrowz/2)+i)];
        end
    elseif Lrowz>=5 
        for i = -1:1:1
            AdditionalElementsIndx = [AdditionalElementsIndx; rowz(floor(Lrowz/2)+i) colz(floor(Lrowz/2)+i)];
        end
    else
        AdditionalElementsIndx = [AdditionalElementsIndx; rowz(1) colz(1)];
    end
end

%Take the central pixels from superpixels that have
%very small width or height
if PCA_FLAG == 1
    
    if sizel11>sizel12
        temp = sizel11;
        sizel11 = sizel12;
        sizel12 = temp;
    end
    
    for i=1:size(AdditionalElementsIndx,1)
        X = [X; Xpca(AdditionalElementsIndx(i,1)*sizel11 + AdditionalElementsIndx(i,2),1), Xpca(AdditionalElementsIndx(i,1)*sizel11 + AdditionalElementsIndx(i,2),2), Xpca(AdditionalElementsIndx(i,1)*sizel11 + AdditionalElementsIndx(i,2),3), Xpca(AdditionalElementsIndx(i,1)*sizel11 + AdditionalElementsIndx(i,2),4), Xpca(AdditionalElementsIndx(i,1)*sizel11 + AdditionalElementsIndx(i,2),5)];
        gIdx_ref_Init = [gIdx_ref_Init; l1(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2))];
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),1) = 110;
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),2) = 0;
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),3) = 255;
    end
    
    [c_all,~,~] = ComputeCentroids(X,gIdx_ref_Init,C1s,0);
    c_allBlack = c_all;
    Xpca_grayscale = X;
else
    for i=1:size(AdditionalElementsIndx,1)
        X = [X; im4(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),1), im4(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),2), im4(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),3)];
        gIdx_ref_Init = [gIdx_ref_Init; l1(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2))];
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),1) = 110;
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),2) = 0;
        im3(AdditionalElementsIndx(i,1),AdditionalElementsIndx(i,2),3) = 255;
    end
    
    [c_all,~,~] = ComputeCentroids(X,gIdx_ref_Init,C1s,0);
    c_allColor = c_all;
    XColor = X;
end

%Show the initial oversegmentation with magenta lines as boundaries
%and with blue color for the central pixels that are chosen for each
%superpixel
str1 = sprintf('Number Of Clusters Initially : %d',C1s);
%Show with RGB color
show(drawregionboundaries(l1,im3, [255 0 144]),120,str1);
%Show with LAB color
% show(drawregionboundaries(l1,im2, [255 0 144]),110,str1);


clearvars -except Am1 BOUNDARY_FLAG C1 C1s c_all gIdx_ref_Init im1 im2 kmax l1 l2 numberOfInitialClustersEXT X XBlackPCA XColor c_allBlack c_allColor;

%Run the whole segmentation process by uncommenting the next lines
%that uses the agglodip algorithm with a local search for merging the
%superpixels and then print the resulting segmentation of the image
pause;
bistest
Print_Resulting_Image