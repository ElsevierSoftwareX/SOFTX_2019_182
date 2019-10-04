function feat_vol = calculateFeatures2D(feat_list,img_vol,retina_mask,header)

if nargin < 4
    header = [];
end
if nargin < 3
    header = [];
    retina_mask = [];
end
retina_mask = uint8(retina_mask);

feat_names = feat_list(:,1);

if ~isempty(header)
    % The scale of the original data, to be used when calculating smoothing
    % parameters
    scale_orig = [3.87167 6.07151 129.44131];
    scale_new = (1000*[header.ScaleZ header.ScaleX header.Distance]);
    scale_ratio = scale_orig./scale_new;
    if all(abs(scale_ratio-1)) < 0.3
        % If the scale is close to the same, probably spectralis data, so don't
        % change anything (for consistency with BMOE paper work, with a fixed
        % sigma for all data)
        scale_ratio = [1 1 1];
        scale_new = scale_orig;
    end
else
    scale_ratio = [1 1 1];
end

% Volume size
sz = size(img_vol);
if length(sz) == 2
    sz(3) = 1;
end

feat_inds = calc_feat_inds(feat_list);
num_feat = feat_inds(end)-1;
feat_vol = zeros([sz, num_feat],class(img_vol));


% -- X Position of each point --
ind = find(strcmp('xpos',feat_names));
if ~isempty(ind)
    feat_vol(:,:,:,feat_inds(ind)) = repmat(1:sz(2),[sz(1) 1 sz(3)]);
end


% -- Y Position of each point --
ind = find(strcmp('ypos',feat_names));
if ~isempty(ind)
    feat_vol(:,:,:,feat_inds(ind)) = repmat((1:sz(1))',[1 sz(2) sz(3)]);
end


% -- Intensity --
ind = find(strcmp('intensity',feat_names));
if ~isempty(ind)
    feat_vol(:,:,:,feat_inds(ind)) = img_vol;
end

% All intensities withing a neighborhood
ind = find(strcmp('nhood',feat_names));
if ~isempty(ind)
    nhs = feat_list{ind,2};
    hs = ceil((nhs-1)/2);

    colpts = round((-hs:hs)*scale_ratio(2));
    rowpts = round((-hs:hs)*scale_ratio(1));
    if any(diff(colpts) == 0) || any(diff(rowpts) == 0)
        warning('neighborhood intensity feature does not scale correctly with pixel size')
    end
    
    nhv = zeros([sz,nhs^2]);
    k = 1;
    for i = colpts % col
        for j = rowpts % row
            if i < 0
                pv = cat(2,repmat(img_vol(:,1,:),[1 -i 1]),img_vol(:,1:end-(-i),:));
            elseif i > 0
                pv = cat(2,img_vol(:,(i+1):end,:),repmat(img_vol(:,end,:),[1 i 1]));
            else
                pv = img_vol;
            end
                
            if j < 0
                pv = cat(1,repmat(pv(1,:,:),[-j 1 1]),pv(1:end-(-j),:,:));
            elseif j > 0
                pv = cat(1,pv((j+1):end,:,:),repmat(pv(end,:,:),[j 1 1]));
            else
%                 pv = pv;
            end                
            
            nhv(:,:,:,k) = pv; k = k + 1;
        end
    end
    
    feat_vol(:,:,:,feat_inds(ind):feat_inds(ind)+nhs^2-1) = nhv;
    clear nhv pv
end


% -- Average intensity around a pixel --
ind = find(strcmp('int_mean',feat_names));
if ~isempty(ind)
    for i = 1:length(ind)
        img_vol_tf = img_vol;
        if length(feat_list{ind(i),2}) == 2
            n = feat_list{ind(i),2}{1};
            tf = feat_list{ind(i),2}{2};
            if strcmp(tf,'log')
                img_vol_tf = log(img_vol);
                img_vol_tf(isinf(img_vol_tf)) = nan;
                img_vol_tf = impute_nan_col(img_vol_tf);
            elseif strncmp(tf,'exp',3)
                expval = str2double(tf(4:end))/100;        
                img_vol_tf = img_vol.^expval;
            end
        else
            n = feat_list{ind(i),2};
        end
        h = fspecial('average',[n n]);
        feat_vol(:,:,:,feat_inds(ind(i))) = imfilter(img_vol_tf,h,'same','symmetric');
    end
    clear img_vol_tf
end


% -- Sobel Gradient --
ind = find(strcmp('sobel_grad',feat_names));
if ~isempty(ind)
    % Possible parameters 'x', 'y', 'mag'
    f_params = feat_list{ind,2};
    
    h = fspecial('sobel');
    i = 0;
    if any(strcmp('x',f_params))
        feat_vol(:,:,:,feat_inds(ind)) = imfilter(img_vol,h','same','symmetric');
        i = 1;
    end
    if any(strcmp('y',f_params))
        feat_vol(:,:,:,feat_inds(ind)+i) = imfilter(img_vol,h,'same','symmetric');
        i = i + 1;
    end
    if any(strcmp('mag',f_params))
        if i == 2
            feat_vol(:,:,:,feat_inds(ind)+i) = sqrt(feat_vol(:,:,:,feat_inds(ind)+i-1).^2 + feat_vol(:,:,:,feat_inds(ind)+i-2).^2);
        else
            feat_vol(:,:,:,feat_inds(ind)+i) = sqrt(imfilter(img_vol,h','same','symmetric').^2 + imfilter(img_vol,h,'same','symmetric').^2);
        end
        i = i + 1;
    end
    if any(strcmp('mag_y',f_params))
        yg = imfilter(img_vol,h,'same','symmetric');
        feat_vol(:,:,:,feat_inds(ind)+i) = sign(yg).*sqrt(imfilter(img_vol,h','same','symmetric').^2 + yg.^2); 
        clear yg
    end
end

% -- Laplacian --
ind = find(strcmp('laplacian',feat_names));
if ~isempty(ind)
    h = fspecial('laplacian',0);
    feat_vol(:,:,:,feat_inds(ind)) = imfilter(img_vol,h,'same','symmetric');
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
            % Get average of nxn region around each pixel
            h = fspecial('average',n);
    %         h = ones(n,n);
            ave_vol = imfilter(img_vol,h,'same','symmetric');
        end
        
        if strcmp(feat_names{ind(i)},'blk_grad_horiz')
            % pad with replicate endpoints
            padvals = repmat(ave_vol(:,1,:),[1 m 1]);        
            left = cat(2,padvals,ave_vol(:,1:end-m,:));
            padvals = repmat(ave_vol(:,end,:),[1 m 1]);
            right = cat(2,ave_vol(:,(1+m):end,:),padvals);
            
            feat_vol(:,:,:,feat_inds(ind(i))) = right - left;         
        else        
            % pad with replicate endpoints
            padvals = repmat(ave_vol(1,:,:),[m 1 1]);        
            above = cat(1,padvals,ave_vol(1:end-m,:,:));
            padvals = repmat(ave_vol(end,:,:),[m 1 1]);
            below = cat(1,ave_vol((1+m):end,:,:),padvals);

            if strcmp(feat_names{ind(i)},'blk_grad_above')
                below = ave_vol;

    %             nhalf = round((n-1)/2);
    %             padvals = repmat(ave_vol(end,:,:),[nhalf 1 1]);
    %             below = cat(1,ave_vol((1+nhalf):end,:,:),padvals);
            elseif strcmp(feat_names{ind(i)},'blk_grad_below')
                above = ave_vol;

    %             nhalf = round((n-1)/2);
    %             padvals = repmat(ave_vol(1,:,:),[nhalf 1 1]);        
    %             above = cat(1,padvals,ave_vol(1:end-nhalf,:,:));
            end

            feat_vol(:,:,:,feat_inds(ind(i))) = above - below;
        end
        n_prev = n;
    end
    clear above below left right ave_vol
end


% -- Neighborhood average of gradient values
ind = find(strncmp('grad_ave',feat_names,8));
if ~isempty(ind)
    n_prev = -1;
    grad_vol = imfilter(img_vol,fspecial('sobel'),'replicate');
    for i = 1:length(ind)
        % Required parameters {n,m}
        f_params = feat_list{ind(i),2};

        if length(f_params) ~= 2
            error('blk_grad requires 2 parameters as input')
        end
        
        n = f_params{1};
        n1 = round(n*scale_ratio(1));
        n2 = round(n*scale_ratio(2));
        m = round(f_params{2}*scale_ratio(1));
        
        % n should be an odd number
        if mod(n1,2) == 0
            n1 = n1 + 1;
        end
        if mod(n2,2) == 0
            n2 = n2 + 1;
        end

%         % n should be an odd number
%         if mod(n,2) == 0
%             error('blk_grad first parameter must be odd')
%         end
        
        % Only recalculate if size of kernel is different
        if n ~= n_prev
            % Get average of nxn region around each pixel
            h = fspecial('average',[n1 n2]);
    %         h = ones(n,n);
            ave_vol = imfilter(grad_vol,h,'same','symmetric');
        end
        
        % pad with replicate endpoints
        if m == 0
            feat_vol(:,:,:,feat_inds(ind(i))) = ave_vol;
        elseif m > 0
            % Take region above current voxel
            padvals = repmat(ave_vol(1,:,:),[m 1 1]);        
            feat_vol(:,:,:,feat_inds(ind(i))) = cat(1,padvals,ave_vol(1:end-m,:,:));
        else
            m = -m;
            padvals = repmat(ave_vol(end,:,:),[m 1 1]);
            feat_vol(:,:,:,feat_inds(ind(i))) = cat(1,ave_vol((1+m):end,:,:),padvals);
        end
        
        n_prev = n;
    end
    clear grad_vol ave_vol
end


% -- Anisotropic gaussian gradient
ind = find(strcmp('anigauss',feat_names));
if ~isempty(ind)
    for i = 1:length(ind)
        % Required parameters {u,v,d}
        f_params = feat_list{ind(i),2};
        
        if length(f_params) < 3
            error('anigauss requires at least 3 parameters as input')
        end
        
        % angle (relative to image coordinate system, y-axis pointing down)
        d = f_params{3};
        % minor-direction
        u = f_params{1}; % approximate scaling for small angles only
        % major-direction
        v = f_params{2};
        % Scale adjustments - convert u and v to microns and then back to
        %   the new scale. Since u and v are rotated, their values are 
        %   relative to a mixture of the pixel size
        if ~all(scale_ratio == 1)
            t = sqrt((cosd(d)*u*scale_orig(1))^2 + (sind(d)*u*scale_orig(2))^2);
            w = sqrt((cosd(d)*v*scale_orig(2))^2 + (sind(d)*v*scale_orig(1))^2);
            d2 = atand(scale_orig(1)/scale_orig(2)*tand(d)); % angle microns
            u = sqrt((cosd(d2)*t/scale_new(1))^2 + (sind(d2)*t/scale_new(2))^2);
            v = sqrt((cosd(d2)*w/scale_new(2))^2 + (sind(d2)*w/scale_new(1))^2);
            d = atand(scale_new(2)/scale_new(1)*tand(d2));
        end
        
        % Old way
%         d = atand(scale_ratio(1)/scale_ratio(2)*tand(d));        
%         u = scale_ratio(1)*u; % approximate scaling for small angles only?
%         v = scale_ratio(2)*v;
        
        % derivative in the y-direction
        if length(f_params) >= 4
            dy = f_params{4};
        else
            dy = 1;
        end
        % Take absolute value or not
        if length(f_params) < 5
            absv = false;
        else
            if strcmp(f_params(5),'abs')
                absv = true;
            else
                absv = false;
            end
        end
        
        dx = 0;
        
        for j = 1:size(feat_vol,3)
            af = zeros(size(img_vol,1),size(img_vol,2),length(d));
            for k = 1:length(d)
                af(:,:,k) = anigauss_mex(double(img_vol(:,:,j)),u,v,d(k),dy,dx);
            end
            if ~absv
                % Signed max
                [~, ci] = max(abs(af),[],3);
                inds = (1:size(img_vol,1)*size(img_vol,2))' + (ci(:)-1)*size(img_vol,1)*size(img_vol,2);
                feat_vol(:,:,j,feat_inds(ind(i))) = reshape(af(inds),size(img_vol,1),size(img_vol,2));
            else
                % Max of absolute value
                mv = max(abs(af),[],3);
                feat_vol(:,:,j,feat_inds(ind(i))) = mv;
            end
        end
    end
end

% -- Distance to top or bottom of retina mask --
% Parameters: 'top','bot','norm'
ind = find(strcmp('retina_dist',feat_names));
if ~isempty(ind)
    if isempty(retina_mask)
        error('retina_dist requires a retina mask volume as input!')
    end
    retina_mask = int8(retina_mask); % allow diff(.) to be negative
    for i = 1:length(ind)        
        f_params = feat_list{ind(i),2};

        retina_mask_diff = (retina_mask+1).*cat(1,diff(retina_mask,1),zeros(1,size(retina_mask,2),size(retina_mask,3))); 
        dist_mask = zeros(size(retina_mask));
        for j = 1:size(retina_mask,3)
            rm = retina_mask_diff(:,:,j);
            [x1 y1] = find(rm==1);
    %         [x2 y2] = find(rm==2); % distance to center
            [x2 y2] = find(rm==-6); % distance to bottom

            if length(y1)~=size(rm,2)
                % Not enough ILM points
                error('Size mismatch for retinal mask')
            end
            if length(y2)~=size(rm,2)
                % Missing BM points, maybe ISOS went down too far
                
                % First try at fix, extend BM across if it is flattened
                if all(x2 == x2(1))
                    x2(end+1:size(rm,2)) = x2(1);
                else
                    % Retina mask is bad
                    error('Size mismatch for retinal mask')
                end                
            end
            dists = x2-x1;
            if strcmp(f_params,'bot')
                % Signed distance to bottom of mask
                rm = -bsxfun(@minus,repmat((1:sz(1))',[1 sz(2)]),x2');
            elseif strcmp(f_params,'top')
                rm = bsxfun(@minus,repmat((1:sz(1))',[1 sz(2)]),x1');
            else 
                % Relative distance from bottom of mask to the top
                rm = bsxfun(@rdivide,-bsxfun(@minus,repmat((1:sz(1))',[1 sz(2)]),x2'),dists');
%                 rm(rm>1) = -1;
%                 rm(rm<0) = -1;
            end
            
            dist_mask(:,:,j) = rm;
        end 
        feat_vol(:,:,:,feat_inds(ind(i))) = dist_mask;
    end
    clear dist_mask retina_mask_diff
end


% -- Binary retina mask --
ind = find(strcmp('retina_binary',feat_names));
if ~isempty(ind)
    if isempty(retina_mask)
        error('retina_binary requires a retina mask volume as input!')
    end
    feat_vol(:,:,:,feat_inds(ind)) = retina_mask > 0;
end


% -- Standard deviation if the pixels in an nxn region around each pixel --
ind = find(strcmp('stdev',feat_names));
if ~isempty(ind)
    n = feat_list{ind,2};
    h = fspecial('average',n);
    std_vol = sqrt(imfilter(img_vol.^2,h,'same','symmetric') - ...
              (imfilter(img_vol,h,'same','symmetric')).^2);
    feat_vol(:,:,:,feat_inds(ind)) = std_vol;
end

function feat_inds = calc_feat_inds(feat_list)
% Note: assume valid parameters

feat_inds = 1:(size(feat_list,1)+1);

for i = 1:size(feat_list,1)
    if strcmp(feat_list{i,1},'sobel_grad')
        f_params = feat_list{i,2};
        feat_inds(i+1:end) = feat_inds(i+1:end)+(length(f_params)-1);
    elseif strcmp(feat_list{i,1},'nhood')
        f_params = feat_list{i,2};
        feat_inds(i+1:end) = feat_inds(i+1:end)+f_params.^2-1;
    end            
end