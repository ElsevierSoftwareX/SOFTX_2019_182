function T = createDataTerm(c)
% Create the data term for input into Boykov max flow min cut code
% (maxflow.m). T will be an Nx2 matrix where N is the number of nodes in
% the graph. Column 1 is non-zero at the nodes connected to the source and
% column 2 is non-zero at the nodes connected to the sink.
%
% c is a matrix of costs at each pixel or voxel, which are converted into
% weights in the graph according to equation 1 in LI08

% Note: works for 2D and 3D graphs

% Reference:
%   LI08 - Kang Li, Xiaodong Wu, Danny Z. Chen, and Milan Sonka. "Optimal
%       Surface Segmentation in Volumetric Images - A Graph-Theoretic
%       Approach." IEEE Transactions on Pattern Analysis and Machine
%       Intelligence (TPAMI). 28(1). January 2006. pp. 119-134.   

n = numel(c);
T = zeros(n,2);

%% Make w term  (LI08 - Equation 1)

w = c(1:end-1,:,:,:) - c(2:end,:,:,:);

% Translate the cost of the base set to be negative (see LI08 section 4.1)
bc = sum(c(end,:)); % Sum of cost of base nodes
if bc > 0
    % Subtract bc+1 from cost of any base node, this avoids an empty
    % zero-cost cut (Note: this doesnt seem to work? Instead just set the
    % weight of all base nodes to -1) 
    c(end) = c(end) - (bc+1);
end
% w = cat(1,w,c(end,:,:,:));

% Negative cost for base set by setting all base nodes to -1
w = cat(1,w,-1*ones(size(c(end,:,:,:))));


%% Create source sink edges (LI08 - Section 4.2)

% Get positive and negative weights
wneg = w < 0;
wpos = ~wneg;
indneg = find(wneg);
indpos = find(wpos);
% Connect to source or sink
T(indneg,1) = -w(indneg);
T(indpos,2) = w(indpos);
% Make sparse for maxflow.m
T = sparse(T);
              
end