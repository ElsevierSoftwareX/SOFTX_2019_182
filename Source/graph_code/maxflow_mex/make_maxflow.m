%MAKE   Compiles the maxflowmex library.
%   
%   (c) 2008 Michael Rubinstein, WDI R&D and IDC
%   $Revision: 1.1 $
%   $Date: 2013/11/15 00:19:29 $
%

mex -largeArrayDims maxflowmex.cpp maxflow-v3.0/graph.cpp maxflow-v3.0/maxflow.cpp
