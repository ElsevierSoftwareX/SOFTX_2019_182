function feat_vol = calculateFeatures3D(feat_list,img_vol,retina_mask,header)

% Note that this file is separate from calculateFeatures2D since we need
%   need the full volume to calculate 3D features, while we only need
%   individual images from the volume to calculate 2D features
% An alternative to this is to input only which slices we want to calculate
%   features on and combine the 2 files

if nargin < 3
    header = [];
end

feat_names = feat_list(:,1);

% Volume size
sz = size(img_vol);
if length(sz) == 2
    sz(3) = 1;
end

feat_inds = calc_feat_inds(feat_list);
num_feat = feat_inds(end)-1;
feat_vol = zeros([sz, num_feat],class(img_vol));


% -- Z Position of each point (slice number) --
ind = find(strcmp('zpos',feat_names));
if ~isempty(ind)
    feat_vol(:,:,:,feat_inds(ind)) = bsxfun(@times,ones(sz),shiftdim(1:sz(3),-1));
end


% -- Distance to the auto-detected fovea center (as detected by the retina
%    mask), can be signed or unsigned depending on which side of the retina
%    you are on -- 
ind = find(strncmp('fovea_dist',feat_names,10));
if ~isempty(ind)
    retina_mask = retina_mask>0;
    
    % Central retina mask
    cr = round([0.3 0.7]*size(retina_mask,3));
    cr2 = round([0.3 0.7]*size(retina_mask,2));
    rm = retina_mask(:,cr2(1):cr2(2),cr(1):cr(2));

    % ILM, ISOS, BM Boundaries
    pts = convertLabelImageToPoints(rm);
    % Fovea thickness
    d = pts(:,:,2) - pts(:,:,1);
    [~,ptm] = min(d(:));
    [ptm1 ptm2] = ind2sub(size(d),ptm);

    % Fovea center (adding back the crop)
    f_c = [ptm1 + cr2(1)-1, ptm2 + cr(1)-1];
    
    [~,X,Z] = ndgrid(1:size(img_vol,1),1:size(img_vol,2),1:size(img_vol,3));
    X = (f_c(1) - X)*header.ScaleX; 
    Z = (f_c(2) - Z)*header.Distance;
    
%     X = (f_c(1) - X)./size(img_vol,2);
%     Z = (f_c(2) - Z)./size(img_vol,3);
    
    if strcmp('fovea_dist_signed',feat_names)
        % Take the sign of which side you are on (left or right of fovea)
        feat_vol(:,:,:,feat_inds(ind)) = sign(X).*sqrt(X.^2 + Z.^2);
    elseif strcmp('fovea_dist_xz',feat_names)
        feat_vol(:,:,:,feat_inds(ind):feat_inds(ind)+1) = cat(4,X,Z);
    else
        feat_vol(:,:,:,feat_inds(ind)) = sqrt(X.^2 + Z.^2);
    end
end


% -- Average intensity around a pixel --
ind = find(strncmp('int_mean',feat_names,8));
if ~isempty(ind)
    for i = 1:length(ind)
        f_params = feat_list{ind(i),2};
        n = f_params;

        h = fspecial3('average',n);
        ave_vol = imfilter(img_vol,h,'same','symmetric');

        feat_vol(:,:,:,feat_inds(ind(i))) = ave_vol;
    end
end


% -- Median intensity around a pixel --
ind = find(strcmp('int_med',feat_names));
if ~isempty(ind)
    f_params = feat_list{ind,2};
    feat_vol(:,:,:,feat_inds(ind)) = medfilt3(img_vol,f_params);
end


% -- Standard deviation of voxels in an nxnxn region around each voxel --
ind = find(strcmp('stdev',feat_names));
if ~isempty(ind)
    n = feat_list{ind,2};
    h = fspecial3('average',n);
    std_vol = sqrt(imfilter(img_vol.^2,h,'same','symmetric') - ...
              (imfilter(img_vol,h,'same','symmetric')).^2);
    feat_vol(:,:,:,feat_inds(ind)) = std_vol;
end


% -- Sobel Gradient --
ind = find(strcmp('sobel_grad',feat_names));
if ~isempty(ind)
    % Possible parameters 'x', 'y', 'mag'
    f_params = feat_list{ind,2};
    
    hs = [1 2 1;2 4 2;1 2 1];
    h = zeros(3,3,3);
    h(:,1,:) = hs; h(:,3,:) = -hs;
    gx = imfilter(img_vol,h,'same','symmetric');
    h = zeros(3,3,3);
    h(1,:,:) = hs; h(3,:,:) = -hs;
    gy = imfilter(img_vol,h,'same','symmetric');
    h = zeros(3,3,3);
    h(:,:,1) = hs; h(:,:,3) = -hs;
    gz = imfilter(img_vol,h,'same','symmetric');
    
    feat_vol(:,:,:,feat_inds(ind)) = sqrt(gx.^2 + gy.^2 + gz.^2);
%     feat_vol(:,:,:,feat_inds(ind)) = sign(gy).*sqrt(gx.^2 + gy.^2 + gz.^2);
end

% -- Block based gradient --
%       Difference of mean intensity in nxn blocks m pixels above and below
%       the current pixel
ind = find(strncmp('blk_grad',feat_names,8));
if ~isempty(ind)
    n_prev = -1;
    for i = 1:length(ind)
        % Required parameters {n,m}
        f_params = feat_list{ind(i),2};

        if length(f_params) ~= 2
            error('blk_grad requires 2 parameters as input')
        end

        n = f_params{1};
        m = f_params{2};

        % n should be an odd number
        if mod(n,2) == 0
            error('blk_grad first parameter must be odd')
        end
        
        % Only recalculate if size of kernel is different
        if n ~= n_prev
            % Get average of nxnxn region around each pixel
            h = fspecial3('average',n);
    %         h = ones(n,n);
            ave_vol = imfilter(img_vol,h,'same','symmetric');
        end
        
        % pad with replicate endpoints
        padvals = repmat(ave_vol(1,:,:),[m 1 1]);        
        above = cat(1,padvals,ave_vol(1:end-m,:,:));
        padvals = repmat(ave_vol(end,:,:),[m 1 1]);
        below = cat(1,ave_vol((1+m):end,:,:),padvals);
        
        if strcmp(feat_names{ind(i)},'blk_grad_above')
%             below = ave_vol;
            
            nhalf = round((n-1)/2);
            padvals = repmat(ave_vol(end,:,:),[nhalf 1 1]);
            below = cat(1,ave_vol((1+nhalf):end,:,:),padvals);
        elseif strcmp(feat_names{ind(i)},'blk_grad_below')
%             above = ave_vol;
            
            nhalf = round((n-1)/2);
            padvals = repmat(ave_vol(1,:,:),[nhalf 1 1]);        
            above = cat(1,padvals,ave_vol(1:end-nhalf,:,:));
        end
        
        feat_vol(:,:,:,feat_inds(ind(i))) = above - below;
        n_prev = n;
    end
    clear above below ave_vol
end


function feat_inds = calc_feat_inds(feat_list)
% Note: assume valid parameters

feat_inds = 1:(size(feat_list,1)+1);

for i = 1:size(feat_list,1)
    if strcmp(feat_list{i,1},'fovea_dist_xz')
        feat_inds(i+1:end) = feat_inds(i+1:end)+1;
    end
end

