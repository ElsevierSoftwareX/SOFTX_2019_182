function A = createEsEdges3D(A,Z,X,Y,N,param,DEF_INF,mask)

if nargin < 8
    mask = true(Z,X,Y,N);
end
if nargin < 7
    DEF_INF = inf;
end

if isempty(A)
    A = sparse(N*X*Y*Z,N*X*Y*Z);
end

if numel(param.dl) == 1
%     param.dl = repmat(param.dl,X,N-1);
    param.dl = repmat(param.dl,[X,Y,N-1]);
end
dl = param.dl;
if numel(param.du) == 1
%     param.du = repmat(param.du,X,N-1);
    param.du = repmat(param.du,[X,Y,N-1]);
end
du = param.du;

mask2 = mask;

%% dl edges
% Vertices send flow from surface i+1 to surface i dl voxels upward 

% List of all vertices (first surface removed)
vtxSet = (X*Z*Y+1):N*X*Z*Y;
% Remove vertices at the top of the image
vtxRemove = false(Z,N-1,X,Y);
for i = 1:N-1
    for j = 1:max(max(dl(:,:,i),[],1),[],2)
        vtxRemove(j,i,dl(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% vtxSet(vtxRemove(:)) = [];

mask2 = mask;
mask2(:,:,:,1) = [];
mask2(vtxRemove) = false;
vtxSet(~mask2(:)) = [];

dl = repmat(shiftdim(dl,-1),[Z 1 1 1]);
% dl(vtxRemove(:)) = [];
dl(~mask2(:)) = [];

A = A + sparse(vtxSet,vtxSet-Z*X*Y-dl,DEF_INF,N*X*Y*Z,N*X*Y*Z);

%% du edges
% Vertices send flow from surface i to surface i+1 du voxels downward 

% List of all vertices (last surface removed)
vtxSet = 1:(N-1)*X*Y*Z;
% Remove vertices at the bottom of the image
vtxRemove = false(Z,N-1,X,Y);
for i = 1:N-1
    for j = 1:max(max(du(:,:,i),[],1),[],2)
        vtxRemove(end-j+1,i,du(:,:,i)>=j) = true;
    end
end
vtxRemove = permute(vtxRemove,[1 3 4 2]);
% vtxSet(vtxRemove(:)) = [];

mask(:,:,:,end) = [];
mask(vtxRemove) = false;
vtxSet(~mask(:)) = [];

du = repmat(shiftdim(du,-1),[Z 1 1 1]);
% du(vtxRemove(:)) = [];
du(~mask(:)) = [];

A = A + sparse(vtxSet,vtxSet+Z*X*Y+du,DEF_INF,N*X*Z*Y,N*X*Z*Y);
