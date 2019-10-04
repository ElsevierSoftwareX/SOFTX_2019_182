function A = createErEdges(A,Z,X,N,param,DEF_INF)
% Create edges in the graph between A-scans (see Equation 3 in LI06). Edges
% are created between each vertex and the vertices DxL pixels down and one
% pixel to the right and DxU pixels down and one to the left. These edges
% act as smoothness constraints where DxL is the maximum allowed decrease
% in the surface as you move in the X direction, while DxU is the maximum
% allowed increase in the surface as you move in the X direction. For a
% visual reference, see Figure 3(b) in GAR09.
%
% Inputs: 
%   A - Input adjacency matrix to add Er edges to
%   Z - Number of rows in the image
%   X - Number of columns in the image
%   N - Number of surfaces to find
%   param - struct with entries DxL and DxU. DxL and DxU should be either
%           1x1 or XxN matrices.
%   DEF_INF - definition of infinity for edge weights (this matters if the
%             min flow max cut code requires integer costs or not
% Output:
%   A - The updated graph with non-zero entries as defined by the Er edges

% References:
%   LI06 - Kang Li, Xiaodong Wu, Danny Z. Chen, and Milan Sonka. "Optimal
%       Surface Segmentation in Volumetric Images - A Graph-Theoretic
%       Approach." IEEE Transactions on Pattern Analysis and Machine
%       Intelligence. 28(1). January 2006. pp. 119-134.
%   GAR09 - Mona Garvin, Michael Abramoff, Xiaodong Wu, Stephen Russell,
%       Trudy Burns, and Milan Sonka. "Automated 3-D Intraretinal Layer
%       Segmentation of Macular Spectral-Domain Optical Coherence
%       Tomography Images." IEEE Transactions on Medical Imaging. 28(9).
%       September 2009. pp. 1436-1447

if nargin < 6
    DEF_INF = inf;
end

if numel(param.DxL) == 1
    param.DxL = repmat(param.DxL,X,N);
end
DxL = param.DxL;
if numel(param.DxU) == 1
    param.DxU = repmat(param.DxU,X,N);
end
DxU = param.DxU;

%% DxL edges
% Vertices send flow DxL pixels down and one pixel right

% List of all vertices
vtxSet = 1:N*X*Z;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,X,N);
for i = 1:N
    for j = 1:max(DxL(:,i))
        vtxRemove(end-j+1,DxL(:,i)>=j,i) = true;
    end
end
% Vertices to remove on the right edge of the image
vtxRemove(:,end,:) = true;
vtxSet(vtxRemove(:)) = [];

DxL = repmat(shiftdim(DxL,-1),[Z 1 1]);
DxL(vtxRemove(:)) = [];

A = A + sparse(vtxSet,vtxSet+Z+DxL,DEF_INF*ones(1,numel(vtxSet)),N*X*Z,N*X*Z);

%% DxU edges 
% Vertices send flow DxU pixels down and one pixel to the left

% List of all vertices
vtxSet = 1:N*X*Z;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,X,N);
for i = 1:N
    for j = 1:max(DxU(:,i))
        vtxRemove(end-j+1,DxU(:,i)>=j,i) = true;
    end
end
% Vertices to remove on the right edge of the image
vtxRemove(:,1,:) = true;
vtxSet(vtxRemove(:)) = [];

DxU = repmat(DxU,[Z 1 1]);
DxU(vtxRemove(:)) = [];

A = A + sparse(vtxSet,vtxSet-Z+DxU,DEF_INF*ones(1,numel(vtxSet)),N*X*Z,N*X*Z);
