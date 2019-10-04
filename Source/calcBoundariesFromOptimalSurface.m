function bd_points = calcBoundariesFromOptimalSurface(prob_vol,params)
% Calculate the optimal set of surfaces from a 4D volume of probability
% values. We fit the surfaces in a sequential manner to save time and
% memory. Since it is the easiest boundary to find, the ILM is fit first.
% Next, we find the 2nd to 5th boundaries (RNFL-GCL, IPL-INL, INL-OPL,
% OPL-ONL) since they are inseparable at the fovea. Finally, we fit the
% last 4 boundaries (ELM, IS-OS, OS-RPE, BM) since they are close together
% in the outer retina.
%
% Inputs:
%   prob_vol - an MxNxKxB volume containing the probability of voxel
%     (m,n,k) belonging to boundary b. Right now a fixed set of
%     boundaries must be included (B = 9), but that may change in the
%     future.
%   params - a structure containing the constraint parameters of the
%     algorithm
% Output: 
%   bd_points - an NxKxB volume of points representing the surface of each
%     boundary. For example, the location of surface b at lateral position
%     n and B-scan k is prob_vol(bd_points(n,k,b),n,k,b).

if size(prob_vol,4) ~= 9
    error('Probability volume must have 9 values!')
end

% Transform values so that the max value becomes the min and vice versa
% (necessary for min-cut/max-flow optimization) - note that the max is not
% necessarily 1 since we do smoothing before hand, but will still be
% non-negative 
prob_vol = max(prob_vol(:)) - prob_vol;

% Get the three regions we will segment separately
prob_ilm = prob_vol(:,:,:,1);
prob_vol(:,:,:,1) = [];

prob_inner = prob_vol(:,:,:,1:4);
prob_vol(:,:,:,1:4) = [];

prob_outer = prob_vol(:,:,:,1:4);
clear prob_vol

bd_points = zeros(size(prob_ilm,2),size(prob_ilm,3),9);

%% Start with the top boundary (ILM)

fprintf('Finding ILM boundary...')

% - Get mask for non-zero values
th = 0.9*max(prob_ilm(:)); % threshold
ilm_mask = get_mask(prob_ilm,th);

% - Construct sparse matrix of node costs
fprintf('creating data term...')
t1 = tic;
T = createDataTerm2(prob_ilm,ilm_mask);
tn = toc(t1);

sz = size(prob_ilm);
if length(sz) == 3
    sz(4) = 1;
end

% - Set graph parameters
params_n = params;
params_n.DxL = params_n.DxL(:,:,1);
params_n.DxU = params_n.DxU(:,:,1);
params_n.DyL = params_n.DyL(:,:,1);
params_n.DyU = params_n.DyU(:,:,1);
params_n.du = [];
params_n.dl = [];

% - Construct graph
fprintf('constructing graph edges...\n\t')
t1 = tic;
A = buildGraphWithSmoothnessConstraints(sz(1),sz(2),sz(3),sz(4),params_n,ilm_mask);
tg = toc(t1);

s = checkMemory(A);
if ~s
    return
end

% - Solve for optimal surface
fprintf('solving min cut/max flow...')
t1 = tic; 
[~, labels] = maxflow(A,T);
tm = toc(t1);

% - Extract Boundary from min closed set
fprintf('extracting boundary...')
lbls = true(size(ilm_mask));
lbls(ilm_mask) = labels;
ilm_pts = extractBoundaryMulti3D_new(lbls);
fprintf('done!\n')

% Remove any potential missing surface points (can this happen due to
%   problems with the mask?) 
if any(isnan(ilm_pts(:)))
    ilm_pts = inpaint_nans(ilm_pts);
end

bd_points(:,:,1) = ilm_pts;

%% Solve for the next 4 boundaries (RNFL-GCL, IPL-INL, INL-OPL,
%    OPL-ONL)

fprintf('Finding inner boundaries...')

% - Mask out areas by looking at the constraints relative to the ILM (by
%     looking at minimum and maximum thickness of each layer)
inner_mask = false(size(prob_inner));
% Sum of lower bounds is the minimum distance from ILM
dlcs = repmat(bd_points(:,:,1),[1 1 4]) + cumsum(params.dl(:,:,1:4),3) - 1;
% Sum of upper bounds is the maximum distance from ILM
ducs = repmat(bd_points(:,:,1),[1 1 4]) + cumsum(params.du(:,:,1:4),3) + 1;
for ii = 1:4
    imi = inner_mask(:,:,:,ii);
    for i = 1:size(inner_mask,1)
        imi(i,dlcs(:,:,ii) >= i) = false;
        imi(i,(dlcs(:,:,ii) < i) & (ducs(:,:,ii) > i)) = true;
        imi(i,ducs(:,:,ii) <= i) = false;
    end
    inner_mask(:,:,:,ii) = imi;
end

% - Construct sparse matrix of node costs
fprintf('creating data term...')
t1 = tic;
T = createDataTerm2(prob_inner,inner_mask);
tn = tn + toc(t1);

sz = size(inner_mask);
if length(sz) == 3
    sz(4) = 1;
end

% - Set graph parameters
params_n = params;
params_n.DxL = params_n.DxL(:,:,2:5);
params_n.DxU = params_n.DxU(:,:,2:5);
params_n.DyL = params_n.DyL(:,:,2:5);
params_n.DyU = params_n.DyU(:,:,2:5);
params_n.du = params_n.du(:,:,2:4);
params_n.dl = params_n.dl(:,:,2:4);

% - Construct graph
fprintf('constructing graph edges...\n\t')
t1 = tic;
A = buildGraphWithSmoothnessConstraints(sz(1),sz(2),sz(3),sz(4),params_n,inner_mask);
tg = tg + toc(t1);

s = checkMemory(A);
if ~s
    return
end

% - Solve for optimal surface
fprintf('solving min cut/max flow...')
t1 = tic; 
[~, labels] = maxflow(A,T);
tm = tm + toc(t1);

% - Extract Boundary from min closed set
fprintf('extracting boundary...')
lbls = ones(size(inner_mask));
lbls(inner_mask) = labels;
inner_pts = extractBoundaryMulti3D_new(lbls);
fprintf('done!\n')

% Remove any potential missing surface points (can happen due to problems
%   with the mask)
for i = 1:4
    ipi = inner_pts(:,:,i);
    if any(isnan(ipi(:)))
        inner_pts(:,:,i) = inpaint_nans(ipi);
    end
end

bd_points(:,:,2:5) = inner_pts;

%% Solve for the last 4 boundaries (ELM, IS-OS, OS-RPE, BM)

fprintf('Finding outer boundaries...')

% - Mask out areas by looking at the constraints relative to the OPL-ONL
%     (by looking at minimum and maximum thickness of each layer)
outer_mask = false(size(prob_inner));
% Sum of lower bounds is the minimum distance from OPL-ONL
dlcs = repmat(bd_points(:,:,5),[1 1 4]) + cumsum(params.dl(:,:,5:8),3) - 1;
% Sum of upper bounds is the maximum distance from OPL-ONL
ducs = repmat(bd_points(:,:,5),[1 1 4]) + cumsum(params.du(:,:,5:8),3) + 1;
for ii = 1:4
    omi = outer_mask(:,:,:,ii);
    for i = 1:size(outer_mask,1)
        omi(i,dlcs(:,:,ii) >= i) = false;
        omi(i,(dlcs(:,:,ii) < i) & (ducs(:,:,ii) > i)) = true;
        omi(i,ducs(:,:,ii) <= i) = false;
    end
    outer_mask(:,:,:,ii) = omi;
end

% - Construct sparse matrix of node costs
fprintf('creating data term...')
t1 = tic;
T = createDataTerm2(prob_outer,outer_mask);
tn = tn + toc(t1);

sz = size(outer_mask);
if length(sz) == 3
    sz(4) = 1;
end

% - Set graph parameters
params_n = params;
params_n.DxL = params_n.DxL(:,:,6:9);
params_n.DxU = params_n.DxU(:,:,6:9);
params_n.DyL = params_n.DyL(:,:,6:9);
params_n.DyU = params_n.DyU(:,:,6:9);
params_n.du = params_n.du(:,:,6:8);
params_n.dl = params_n.dl(:,:,6:8);

% - Construct graph
fprintf('constructing graph edges...\n\t')
t1 = tic;
A = buildGraphWithSmoothnessConstraints(sz(1),sz(2),sz(3),sz(4),params_n,outer_mask);
tg = tg + toc(t1);

s = checkMemory(A);
if ~s
    return
end

% - Solve for optimal surface
fprintf('solving min cut/max flow...')
t1 = tic; 
[~, labels] = maxflow(A,T);
tm = tm + toc(t1);

% - Extract Boundary from min closed set
fprintf('extracting boundary...')
lbls = ones(size(outer_mask));
lbls(outer_mask) = labels;
outer_pts = extractBoundaryMulti3D_new(lbls);
fprintf('done!\n')

% Remove any potential missing surface points (can happen due to problems
%   with the mask)
for i = 1:4
    opi = outer_pts(:,:,i);
    if any(isnan(opi(:)))
        outer_pts(:,:,i) = inpaint_nans(opi);
    end
end

bd_points(:,:,6:9) = outer_pts;

disp(['Graph constructed in ' num2str(tn+tg) ' seconds']);
disp(['Max flow calculated in ' num2str(tm) ' seconds']);

function s = checkMemory(A)

s = true;
if ~isunix
    % Make sure we wont run out of memory
    s = getfield(whos('A'),'bytes'); % size of A
    % Maxflow requires memory of about 4 times the size of A
    [~,sys] = memory;
    mem_available = sys.PhysicalMemory.Available;
    % Probably want to leave at least 1GB to keep it from running sluggish
    if (mem_available-s*4) < 1e9
        answ = questdlg('Running maxflow will exceed system memory, continue?',...
            'Low memory','Yes','No','No'); 
        if strcmp(answ,'No')
            s = false;
        end
    end
end

function mask = get_mask(vol,th)
% Create a mask for the volume where we want to mask out connected low
% probability voxels at the top and bottom of the volume (threshold defined
% by th). We finish by doing a dilation so that the mask is not too close
% to boundary pixels

vol_th = vol > th;
vol_th = min(vol_th,[],4);

% Remove small points by closing
vol_th = imclose(vol_th,strel('line',5,0));

mask_top = false(size(vol_th));
mask_bot = false(size(vol_th));

% Look for zero voxels starting at the top row the volume
fn = vol_th(1,:,:);
mask_top(1,:,:) = fn;
j = 2;
while any(fn(:)) && j <= size(vol_th,1)
    % Get points that were previously and are still zero
    fn = vol_th(j,:,:) & fn;
    mask_top(j,:,:) = fn;
    j = j + 1;
end

% Fill in empty spots
[xe, ye] = ind2sub([size(vol_th,2) size(vol_th,3)],find(squeeze(all(mask_top==1,1))));
ye_u = unique(ye);
for i = 1:length(ye_u)
    inds_bsc = find(ye==ye_u(i));
    ye_b = ye(inds_bsc);
    xe_b = xe(inds_bsc);
    
    xp = 1:size(mask_top,2);
    xp(xe_b) = [];
    for j = 1:length(xe_b)
        [~,min_ind] = min(abs(xp-xe_b(j)));
        tv = find(mask_top(:,xp(min_ind),ye_b(j)),1,'last');
        mask_top((tv+1):end,xe_b(j),ye_b(j)) = 0;
    end
end

% Look for zero voxels starting at the bottom row the volume
nr = size(vol_th,1);
fn = vol_th(nr,:,:);
mask_bot(nr,:,:) = fn;
j = nr-1;
while any(fn(:)) && j >= 1
    % Get points that were previously and are still zero
    fn = vol_th(j,:,:) & fn;
    mask_bot(j,:,:) = fn;
    j = j - 1;
end

% Fill in empty spots
[xe, ye] = ind2sub([size(vol_th,2) size(vol_th,3)],find(squeeze(all(mask_bot==1,1))));
ye_u = unique(ye);
for i = 1:length(ye_u)
    inds_bsc = find(ye==ye_u(i));
    ye_b = ye(inds_bsc);
    xe_b = xe(inds_bsc);
    
    xp = 1:size(mask_top,2);
    xp(xe_b) = [];
    for j = 1:length(xe_b)
        [~,min_ind] = min(abs(xp-xe_b(j)));
        tv = find(mask_bot(:,xp(min_ind),ye_b(j)),1,'first');
        mask_bot(1:(tv-1),xe_b(j),ye_b(j)) = 0;
    end
end

% Combine top and bottom masks
mask = ~(mask_top | mask_bot);

% Dilate across A-scans for safety (we don want breaks in the graph)
st = strel('rectangle',[5 20]);
mask = imdilate(mask,st);

% Final check to make sure all A-scans have a probability
ma = all(mask==0,1);
if any(ma(:))
    error('Mask generation failed! Detected a break in the mask.')
end
