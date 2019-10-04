function A = buildGraphWithSmoothnessConstraints(Z,X,Y,N,param,mask)
% Build a sparse graph A where A(i,j) is non-zero if the graph has an edge
% from vertex i to vertex j. The weights are infinity as stated in LI06 -
% Section 4.2. Note that according to convention in LI06, this function is
% building the graph G = (V,E), not G_st. The source and sink nodes are
% not included in this graph, they are added by the min cut max flow code
% (maxflow.m). 
%
% Inputs:
%   Z - Number of rows in the volume
%   X - Number of columns (A-scans) in the volume
%   Y - Number of B-scans in the volume
%   N - Number of surfaces to segment
%   param - struct with parameters DxL, DxU, DyL, DyU, dl, du representing
%           smoothness constraints (see the functions createErEdges.m and
%           createEsEdges.m for more details)
% Outputs:
%   A - Adjacency matrix for the graph

% Building the graph depends on the size of our volume. Note following
% convention in the paper, a MATLAB volume of dimensions IxJxK corresponds
% to ZxXxY in the paper

% References:
%   LI06 - Kang Li, Xiaodong Wu, Danny Z. Chen, and Milan Sonka. "Optimal
%       Surface Segmentation in Volumetric Images - A Graph-Theoretic
%       Approach." IEEE Transactions on Pattern Analysis and Machine
%       Intelligence (TPAMI). 28(1). January 2006. pp. 119-134.
%   GAR09 - Mona Garvin, Michael Abramoff, Xiaodong Wu, Stephen Russell,
%       Trudy Burns, and Milan Sonka. "Automated 3-D Intraretinal Layer
%       Segmentation of Macular Spectral-Domain Optical Coherence
%       Tomography Images." IEEE Transactions on Medical Imaging. 28(9).
%       September 2009. pp. 1436-1447

% Infinity definition (set weights for all non-zero edges to infinity)
% **Note that many graph cut solvers assume integer data, and therefore,
%   using a double infinity definition will cause an overflow and the code
%   will not work. Boykov's code has no restrictions on input data type.
DEF_INF = inf;
% DEF_INF = double(int32(inf));

if nargin < 6
    mask = true(Z,X,Y,N);
end

%% Create Intra Column Edges Ea

A = createEaEdges(Z,X,Y,N,DEF_INF,mask);

%% Create inter column edges Er
if Y == 1
%     Aer = createErEdges([],Z,X,N,param,DEF_INF);
    
    Aer = createErEdges3D_X([],Z,X,Y,N,param,DEF_INF);
    Aer(~mask(:),:) = [];
    Aer(:,~mask(:)) = [];
    A = A + Aer;
    clear Aer
else
    A = A + createErEdges3D_X([],Z,X,Y,N,param,DEF_INF,mask);
    A = A + createErEdges3D_Y([],Z,X,Y,N,param,DEF_INF,mask);
end

%% Create Inter Surface Edges Es
if N > 1
    if Y == 1
        if N > 1
            Aes = createEsEdges([],Z,X,N,param,DEF_INF);
        end
    else
        A = A + createEsEdges3D([],Z,X,Y,N,param,DEF_INF,mask);
        
%         Aes = createEsEdges3D([],Z,X,Y,N,param,DEF_INF);
%         Aes(~mask(:),:) = [];
%         Aes(:,~mask(:)) = [];
%         A = A + Aes;
%         clear Aes
    end
end

%% Add smoothness penalty between adjacent pixels
% Vertical edges weight 1, horizontal edges sqrt(2)
% Vertical edges weight 0, horizontal edges sqrt(2)-1

param.DxL = zeros(size(param.DxL));
param.DxU = zeros(size(param.DxU));

if ~isfield(param,'wts')
    wt = 0;
else
    wt = param.wts;
end

if wt == 0
    %% Remove zero voxels
    A(~mask(:),:) = [];
    A(:,~mask(:)) = [];
    return
end

A = A + createErEdges3D_X([],Z,X,Y,N,param,wt,mask);

% As = createErEdges3D_X([],Z,X,Y,N,param,wt);
% As(~mask(:),:) = [];
% As(:,~mask(:)) = [];
% A = A + As;

if Y > 1
    param.DyL = zeros(size(param.DyL));
    param.DyU = zeros(size(param.DyU));

%     A = A + createErEdges3D_Y([],Z,X,Y,N,param,wt,mask);
    A = A + createErEdges3D_Y([],Z,X,Y,N,param,wt*Y/X,mask);
    
%     As = createErEdges3D_Y([],Z,X,Y,N,param,wt);
% %     As = createErEdges3D_Y([],Z,X,Y,N,param,wt*Y/X);
%     As(~mask(:),:) = [];
%     As(:,~mask(:)) = [];
%     A = A + As;
end

%% Remove zero voxels
A(~mask(:),:) = [];
A(:,~mask(:)) = [];