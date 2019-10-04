function A = createEsEdges(A,Z,X,N,param,DEF_INF)
% Create edges in the graph between surfaces (See Equation 4 in LI06).
% Edges are created between each vertex and the corresponding vertex in the
% adjacent surface. These edges act as distance constraints between each
% surface, defining the minimum (dl(i,i+1)) and maximum (du(i,i+1))
% distance between surfaces i and i + 1. For a visual reference image, see
% Fig. 3(d) in GAR09.
%
% Inputs:
%   A - Input adjacency matrix to add Es edges to
%   Z - Number of rows in the image
%   X - Number of columns in the image
%   N - Number of surfaces to find
%   param - struct with entries du and dl. DxL and DxU should be either
%           1x1 or Xx(N-1) matrices.
%   DEF_INF - definition of infinity for edge weights (this matters if the
%             min flow max cut code requires integer costs or not
% Output:
%   A - The updated graph with non-zero entries as defined by the Es edges

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

if numel(param.dl) == 1
    param.dl = repmat(param.dl,X,N-1);
end
dl = param.dl;
if numel(param.du) == 1
    param.du = repmat(param.du,X,N-1);
end
du = param.du;

%% dl edges
% Vertices send flow from surface i+1 to surface i dl pixels upward 

% List of all vertices (first surface removed)
vtxSet = (X*Z+1):N*X*Z;
% Remove vertices at the top of the image
vtxRemove = false(Z,X,N-1);
for i = 1:(N-1)
    for j = 1:max(dl(:,i))
        vtxRemove(j,dl(:,i)>=j,i) = true;
    end
end
vtxSet(vtxRemove(:)) = [];

dl = repmat(shiftdim(dl,-1),[Z 1 1]);
dl(vtxRemove(:)) = [];

A = A + sparse(vtxSet,vtxSet-Z*X-dl,DEF_INF*ones(1,numel(vtxSet)),N*X*Z,N*X*Z);

%% du edges
% Vertices send flow from surface i to surface i+1 du pixels downward 

% List of all vertices (last surface removed)
vtxSet = 1:(N-1)*X*Z;
% Remove vertices at the bottom of the image
vtxRemove = false(Z,X,N-1);
for i = 1:N-1
    for j = 1:max(du(:,i))
        vtxRemove(end-j+1,du(:,i)>=j,i) = true;
    end
end
vtxSet(vtxRemove(:)) = [];

du = repmat(shiftdim(du,-1),[100 1 1]);
du(vtxRemove(:)) = [];

A = A + sparse(vtxSet,vtxSet+Z*X+du,DEF_INF*ones(1,numel(vtxSet)),N*X*Z,N*X*Z);

% for p = 1 : N-1
%     i = ijPairs(p,1);   %%Surface i
%     j = ijPairs(p,2);   %%Surface j
%    
%     du = param.du{i,j}; %%Max allowed distance between i,j at x
%     dl = param.dl{i,j}; %%Min allowed distance between i,j at x
%     
%     fromVtxSeti = (p*X*Z+1):((p+1)*X*Z);   %%dl edges are meant to be sent from surface j to surface i in the upward direction
%     topVtcsi = p*X*Z+1:Z:((p+1)*X*Z);      %%Top Set of vertices above which flow is not to be sent
%     topVtxSet = [];                          %%Initialize collective set of vertices above which flow is not to be sent. 
%                                              %%This includes, the topVtcs as well as all vertices upto dl below it 
%     for x = 1 : X
%         topVtxSet = [topVtxSet topVtcsi(x):(topVtcsi(x)+dl(x)-1)];  %#ok<AGROW>
%     end    
%     fromVtxSeti = setdiff(fromVtxSeti,topVtxSet);
%     dl = repmat(dl,Z,1);dl = dl(:)';
%     dl = dl(fromVtxSeti-p*X*Z-1);
% 
%     A = A + sparse(fromVtxSeti,fromVtxSeti-Z*X-dl,DEF_INF*ones(1,numel(fromVtxSeti)),N*X*Z,N*X*Z);
%     
%     fromVtxSetj = ((p-1)*X*Z+1):(p*X*Z);   %%du edges are meant to be sent from surface i to surface j in the downward direction
%     botVtcsj = (p-1)*X*Z+Z:Z:(p*X*Z);      %%Bottom set of vertices below which flow is not to be sent
%     botVtxSet = [];                          %%Initialize collective set of vertices below which flow is not to be sent. 
%                                              %%This includes, the botVtcs as well as all vertices upto du above it 
%     for x = 1 : X
%        botVtxSet = [botVtxSet (botVtcsj(x)-du(x)+1):botVtcsj(x)];  %#ok<AGROW>
%     end  
%     fromVtxSetj = setdiff(fromVtxSetj,botVtxSet);
%     du = repmat(du,Z,1);du = du(:)';
%     du = du(fromVtxSetj-(p-1)*X*Z);
%     A = A + sparse(fromVtxSetj,fromVtxSetj+Z*X+du,DEF_INF*ones(1,numel(fromVtxSetj)),N*X*Z,N*X*Z);
% %     A = A + sparse(fromVtxSetj-Z*X+dl,fromVtxSetj,DEF_INF*ones(1,numel(fromVtxSetj)),N*X*Z+2,N*X*Z+2);    
% end

end