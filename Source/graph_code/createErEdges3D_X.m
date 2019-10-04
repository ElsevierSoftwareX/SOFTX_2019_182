function A = createErEdges3D_X(A,Z,X,Y,N,param,DEF_INF,mask)

%%Inputs:Z,X are image dimensions. N is the number of surfaces, param is
%%the parameter listing.
%%Outputs:Updated smoothness matrix containing Er edges.

if nargin < 8
    mask = true(Z,X,Y,N);
end
if nargin < 7
    DEF_INF = inf;
end

if isempty(A)
    A = sparse(N*X*Y*Z,N*X*Y*Z);
end

if numel(param.DxL) == 1
    param.DxL = repmat(param.DxL,[X,Y,N]);
end
DxL = param.DxL;
if numel(param.DxU) == 1
    param.DxU = repmat(param.DxU,[X,Y,N]);
end
DxU = param.DxU;

%% DxL edges
% Vertices send flow DxL pixels down and one pixel right

% List of all vertices
vtxSet = 1:N*X*Y*Z;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,N,X,Y);
for i = 1:N
    for j = 1:max(max(DxL(:,:,i),[],1),[],2)
        vtxRemove(end-j+1,i,DxL(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% Vertices to remove on the right edge of the image
vtxRemove(:,end,:,:) = true;
% vtxSet(vtxRemove(:)) = [];

mask2 = mask;
mask2(vtxRemove) = false;
vtxSet(~mask2(:)) = [];

DxL = repmat(shiftdim(DxL,-1),[Z 1 1 1]);
% DxL(vtxRemove(:)) = [];
DxL(~mask2(:)) = [];

A = A + sparse(vtxSet,vtxSet+Z+DxL,DEF_INF,N*X*Y*Z,N*X*Y*Z);

%% DxU edges 
% Vertices send flow DxU pixels down and one pixel to the left

% List of all vertices
vtxSet = 1:N*X*Z*Y;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,N,X,Y);
for i = 1:N
    for j = 1:max(max(DxU(:,:,i),[],1),[],2)
        vtxRemove(end-j+1,i,DxU(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% Vertices to remove on the right edge of the image
vtxRemove(:,1,:,:) = true;
% vtxSet(vtxRemove(:)) = [];

mask(vtxRemove) = false;
vtxSet(~mask(:)) = [];

% DxU = repmat(DxU,[Z 1 1 1]);
DxU = repmat(shiftdim(DxU,-1),[Z 1 1 1]);
% DxU(vtxRemove(:)) = [];
DxU(~mask(:)) = [];

A = A + sparse(vtxSet,vtxSet-Z+DxU,DEF_INF,N*X*Y*Z,N*X*Y*Z);
