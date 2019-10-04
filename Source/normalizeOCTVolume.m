function out_vol = normalizeOCTVolume(in_vol,method,header)

if nargin < 2
    method = 2;
    header = [];
elseif nargin < 3
    header = [];
end

if ~isempty(header)
    % The scale of the original data, to be used when calculating smoothing
    % parameters
    scale_orig = [3.87167 6.07151 129.44131];
    scale_ratio = scale_orig./...
                    (1000*[header.ScaleZ header.ScaleX header.Distance]);
    if mean(abs(scale_ratio-1)) < 0.3
        % If the scale is close to the same, probably spectralis data, so don't
        % change anything (for consistency with BMOE paper work, with a fixed
        % sigma for all data)
        scale_ratio = [1 1 1];
    end
else
    scale_ratio = [1 1 1];
end

% Median filter kernel
h = [round(15*scale_ratio(1)) 1];

out_vol = zeros(size(in_vol),class(in_vol));
if method == 1
    % Median filter standard deviation normalization
    
    % Get standard deviation of whole volume, we want each slice to have
    % approximately this standard deviation (which also keeps the
    % intensity range of 0-1 meaningful)
    stdv = nanstd(in_vol(:));
    
    % Loop through each slice, median filter the image and normalize the
    % slice by the standard deviation of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol(:,:,j),h);
        out_vol(:,:,j) = in_vol(:,:,j)./nanstd(med(:))*stdv;
    end
elseif method == 2
    % Median filter contrast stretching normalization
    if isinteger(in_vol)
        mv = double(intmax(class(in_vol))); % imadjust requires double
    else
        mv = 1;
    end
    
    maxOffset = 0.05*mv;
    
    % Loop through each slice, median filter the image and normalize by
    % contrast stretching to the max of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol(:,:,j),h);
        ms = double(max(med(:)))+maxOffset;
        if ms > mv
            ms = mv;
        end
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[0 ms]/mv,[0 1]);
    end    
elseif method == 3
    % Quantile contrast stretching
    qle1 = 0.3;
    qle2 = 0.999;
    
    q1 = quantile(in_vol(:),qle1);
    q2 = quantile(in_vol(:),qle2);
    
    for j = 1:size(in_vol,3)
        sl = in_vol(:,:,j);
%         q1 = quantile(sl(:),qle1);
%         q2 = quantile(sl(:),qle2);
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[q1 q2],[0 1]);
    end    
elseif method == 4
    % Threshold based on otsu
    in_vol = in_vol.^0.25;
    for j = 1:size(in_vol,3)
        sl = in_vol(:,:,j);
        t = graythresh(sl);
        h = hist(sl(sl>t),255)/sum(sl(:)>t);
    end
    
elseif method == 5
    % Method of Girard et al. (2011), Shadow removal and contrast
    % enhancement in OCT images of the human optic nerve head
    cs = flipdim(cumtrapz(flipdim(in_vol,1),1),1);
    cs = imfilter(cs,fspecial('gaussian',19,3),'symmetric');
    out_vol = in_vol./cs;
    out_vol(isnan(out_vol)) = min(out_vol(:)); % 0/0
    
elseif method == 6
    % Median filter contrast stretching normalization with lower threshold
    % Median filter contrast stretching normalization
    if isinteger(in_vol)
        mv = double(intmax(class(in_vol))); % imadjust requires double
    else
        mv = 1;
    end
    
    maxOffset = 0.05*mv;
    minOffset = 0.05;
    
    % Loop through each slice, median filter the image and normalize by
    % contrast stretching to the max of the median filtered image
    for j = 1:size(in_vol,3)
        med = medfilt2(in_vol(:,:,j),h);
        ms = double(max(med(:)))+maxOffset;
        if ms > mv
            ms = mv;
        end
        mns = double(min(med(:)))+minOffset;
        out_vol(:,:,j) = imadjust(in_vol(:,:,j),[mns ms]/mv,[0 1]);
    end    
end

% out_vol(out_vol > 1) = 1;
% out_vol(out_vol < 0) = 0;
out_vol(isnan(in_vol)) = nan;