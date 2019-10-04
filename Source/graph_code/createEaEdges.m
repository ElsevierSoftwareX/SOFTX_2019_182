function A = createEaEdges(Z,X,Y,N,DEF_INF,mask)
% Create edges in the graph between each vertex and the vertex below it.
% See Equation 2 in LI06.
%
% Inputs:
%   Z - Number of rows in the volume
%   X - Number of columns in the volume
%   Z - Number of B-scans in the volume
%   N - Number of surfaces to find
%   DEF_INF - definition of infinity for edge weights (this matters if the
%             min flow max cut code requires integer costs or not
% Output:
%   A - The updated graph with non-zero entries as defined by the Ea edges

% Reference:
%   LI08 - Kang Li, Xiaodong Wu, Danny Z. Chen, and Milan Sonka. "Optimal
%       Surface Segmentation in Volumetric Images - A Graph-Theoretic
%       Approach." IEEE Transactions on Pattern Analysis and Machine
%       Intelligence. 28(1). January 2006. pp. 119-134.  

if nargin < 6
    mask = true(Z,X,Y,N);
end
if nargin < 5
    DEF_INF = inf;
end

% All vertices
vtxSet = 1:N*X*Z*Y;
% Remove base set vertices, they have no edges
vtxSetRemov = Z:Z:N*X*Z*Y;
% vtxSet(vtxSetRemov) = [];

mask(vtxSetRemov) = false;
vtxSet(~mask(:)) = [];

% Create edges between each vertex and the vertex below it
n = N*X*Z*Y;
A = sparse(vtxSet,vtxSet+1,DEF_INF,n,n);