function A = createErEdges3D_Y(A,Z,X,Y,N,param,DEF_INF,mask)

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

if numel(param.DyL) == 1
    param.DyL = repmat(param.DyL,[X,Y,N]);
end
DyL = param.DyL;
if numel(param.DyU) == 1
    param.DyU = repmat(param.DyU,[X,Y,N]);
end
DyU = param.DyU;

%% DyL edges 
% Vertices send flow DyL pixels down and one pixel back

% List of all vertices
vtxSet = 1:N*X*Z*Y;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,N,X,Y);
for i = 1:N
    for j = 1:max(max(DyL(:,:,i),[],1),[],2)
        vtxRemove(end-j+1,i,DyL(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% Remove last B-scan from vertices
vtxRemove(:,:,end,:) = true;
% vtxSet(vtxRemove(:)) = [];

mask2 = mask;
mask2(vtxRemove) = false;
vtxSet(~mask2(:)) = [];

DyL = repmat(shiftdim(DyL,-1),[Z 1 1 1]);
% DyL(vtxRemove(:)) = [];
DyL(~mask2(:)) = [];

A = A + sparse(vtxSet,vtxSet+X*Z+DyL,DEF_INF,N*X*Y*Z,N*X*Y*Z);

%% DyU edges 
% Vertices send flow DyU pixels down and one pixel forward

% List of all vertices
vtxSet = 1:N*X*Z*Y;
% Vertices to remove at the bottom of the image
vtxRemove = false(Z,N,X,Y);
for i = 1:N
    for j = 1:max(max(DyU(:,:,i),[],1),[],2)
        vtxRemove(end-j+1,i,DyU(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% Remove first B-scan from vertices
vtxRemove(:,:,1,:) = true;
% vtxSet(vtxRemove(:)) = [];

mask(vtxRemove) = false;
vtxSet(~mask(:)) = [];

DyU = repmat(shiftdim(DyU,-1),[Z 1 1 1]);
% DyU(vtxRemove(:)) = [];
DyU(~mask(:)) = [];

A = A + sparse(vtxSet,vtxSet-X*Z+DyU,DEF_INF,N*X*Y*Z,N*X*Y*Z);
