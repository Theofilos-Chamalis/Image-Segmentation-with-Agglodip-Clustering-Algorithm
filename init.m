%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%------------
% Initialization for Agglodip Image Segmentation package. 
% Run RunTheSlic.m for a demonstration.
%------------
% Copyright (C) 2014-2015, Chamalis Theofilos.
%------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear;
fclose('all');
clc;

% copyright/copyleft info
fprintf('Agglodip Image Segmentation package v.0.1. Copyright (C) 2014-2015, Chamalis Theofilos.\n'); 
fprintf('This is free software distributed under the GNU General Public License; for details see LICENSE.txt.\n\n');

% add folders and their subfolders to the path
addpath(genpath('Algorithms_Package'));
addpath(genpath('External_Bootstrap'));
addpath(genpath('Images'));
            
% open the files that are required
% to run the agglodip segmentation.
% the whole segmentation process (Preprocessing,
% Segmentation with agglodip and print of the
% final segments) can be run separately or just
% via the Image_Segmentation_Preprocessing.m
% file
edit('Image_Segmentation_Preprocessing.m');
edit('bistest.m');
edit('Print_Resulting_Image');

% main file that call the agglodip algorithm
edit('bisect_agglodip_segmentation.m');

% misc files for image segmentation
%edit('PCA_Mat_Create.m');
%edit('slic.m');
%edit('spdbscan.m');

% files for Hartigan's test
%edit('HartigansDipSignifTest_no_boot.m');
%edit('HartigansDipSignifTest.m');
%edit('HartigansDipTest.m');
%edit('hartigansdiptestdemo.m');

% some of the functions used
%edit('test_projected_unimodality.m');
%edit('test_unimodal_cluster.m');
%edit('test_unimodal_cluster2.m');
%edit('data_Proj');
%edit('MergeSmallClusters.m');
%edit('UnimodalityChecking.m');
%edit('plotClusterResults');
%edit('compute_unifpdfboot');

% pgmeans specific files
%edit('pgmeans_init_params.m');


