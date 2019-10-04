% Build all of the required libraries

addpath(fullfile(pwd,'Source'),...
        fullfile(pwd,'Source/anigaussm'),...
        fullfile(pwd,'Source/graph_code'),...
        fullfile(pwd,'Source/graph_code/maxflow_mex'),...
        fullfile(pwd,'Source/randomforest_base'));

% Graph cut
run('Source/graph_code/maxflow_mex/make_maxflow.m')

% Random forest
if ispc
    run('Source/randomforest_base/compile_windows.m')
else
    run('Source/randomforest_base/compile_linux.m')
end

% Anigauss
currdir = pwd;
cd Source/anigaussm
mex -v -largeArrayDims anigauss_mex.c anigauss.c 

cd(currdir)