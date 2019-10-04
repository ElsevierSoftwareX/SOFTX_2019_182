
function  SegmentationMethod2AMPAO(handles)
if ~strcmp(handles.data.format,'vol')
    errordlg('Method II works only on VOL data','File Error');
    return
end

folderPath = handles.PathName;
% folderPath = 'D:\kafieh\AMPAO';
filterindex = 2;
% folderPath = strcat(folderPath,'\');

% image = double(handles.data.BScan);


%% --------------------Layer Segmentation Lang (Seg_script.m)

% Script showing how to run the segmentation algorithm. Contains a brief
% description of what inputs are necessary to run the algorithm

% clear all; close all;

if ~exist('OCTLayerSegmentation','file')
    addpath(fullfile(pwd,'Source'),...
        fullfile(pwd,'Source/anigaussm'),...
        fullfile(pwd,'Source/graph_code'),...
        fullfile(pwd,'Source/graph_code/maxflow_mex'),...
        fullfile(pwd,'Source/randomforest_base'));
end

% ---------------------------------------- %
% Input data formats (three possibilities)
% ---------------------------------------- %

% 1) string containing the path to single file
PathName = handles.PathName;

FileName = handles.FileName;

% 'C:\Users\user\Documents\MATLAB\VOLdata/voldata.vol';
% [FileName,PathName] =uigetfile( '*.vol')
% 2) cell array containing path to multiple files
%'C:\Users\user\Documents\MATLAB\OCTLayerSegmentation2.11\1\voldata.vol';
% filenames = {'C:/OCTData/OCT_001.vol',...
%              'C:/OCTData/OCT_002.vol'};

% 3) a text file containing a separate file on each line (can be generated
%    using the script 'make_file_list_(scanner).m')

% filenames = 'file_list.txt';


% ---------------------------------------------------------------------- %
% Set algorithm parameters
%
%   Type 'doc OCTLayerSegmentation' for more information about the input
%   parameters
% ---------------------------------------------------------------------- %

% Algorithm specific parameters
params.minseg =true; % if true -> runs faster, slightly less accurate
params.smooth = true; % smooth the output segmentation boundaries
params.segmethod = 1; % use 1 for best results
params.resizedata = true; % resize images to the size used for training

% Display options
params.printtoscreen = true; % display algorithm status in command window
params.displayresult = true; % show result in interactive octViewer GUI
params.displaygrid = false; % display a figure with the ETDRS grid results

% Results options
params.resultfolder = uigetdir(handles.PathName,'Select location to save results');
% params.resultfolder = 'D:\kafieh\students\bahare salafian\codes\Codes\Internship\results'; % folder to save results
params.logfile = true; % save a log file - can be found in the results folder
params.skip_completed = false; % do not run algorithm if result exists already
params.overwrite_results = false; % overwrite results if they already exist
params.saveXMLfiles = false; % save xml results viewable in MIPAV
params.gridradii = [500 1500 2500]; % radius of the ETDRS grid circles (um)

% ------------- %
% Run algorithm
% ------------- %

% stats = thicknessStatistics(boundaryPoints,header,circ_rad,showgrid)

tic
OCTLayerSegmentation2(handles,params);
% OCTLayerSegmentation1(strcat([handles.PathName,handles.FileName]),params);
toc



% uisave( 'imageLayer','SegmentedData')

