function [retinaMask, shifts, boundaries, nbpt] = retinaDetector2_scale(img_vol,header,params,doplots,bsc_indep)
% Detect the retina boundaries. The retina mask will have a value of 0
% outside of the retina, a value of  1 between the ilm and inner segment,
% and a value of 2 between the inner segment and Bruch's membrane (bottom
% of the RPE)
%
% Note, if you want a more accurate ISOS boundary, use retinaDetector()
% code, not retinaDetector2
%
% Shifts is an YxZ vector containing the number of pixels to shift each
% a-scan to flatten the image. Use retinaFlatten.m to flatten the data
% using these shifts.
%
% Boundaries is an YxZx3 array containing the value of the ILM, ISOS and
% BM boundaries

% Based on Bhaskar's original Java code

if nargin < 2
    error('number of arguments must be at least 2!')
elseif nargin < 3
    params.sigma_lat = 16.67;
    params.sigma_ax = 11.6;
    params.distconst = 96.68;
    doplots = false;
    bsc_indep = false;
elseif nargin < 4
    doplots = false;
    bsc_indep = false;
elseif nargin < 5
    bsc_indep = false;
end
sl = 17;
sl = 48;
sl = round(size(img_vol,3)/2);
% sl = 1;

% Maximum distance from ILM to ISOS
maxdist = 386.73; % ~100 pixels in spectralis
% Maximum distance from ISOS to BM
maxdist_bm = 116.02; % ~30 pixels in spectralis
% Minimum distance from ISOS to BM
isosThresh = 20; 
% Median filter outlier threshold distance
if bsc_indep
    dc_thresh = 20;
else
    dc_thresh = 58.01; %~15 pixels in spectralis
end
% Median filter outlier kernel
mf_k = 855.12; % %~7 pixels in spectralis

% Sigma values for smoothing final surfaces
% Through plane direction
sigma_tp_ilm = 91.62; % ~0.75 pixels in spectralis
sigma_tp_isos = 91.62; % ~0.75 pixels in spectralis
sigma_tp_bm = 244.32; % ~2 pixels in spectralis
% Lateral direction
sigma_lat_ilm = 55.56; % ~10 pixels in spectralis
sigma_lat_isos = 55.56; % ~10 pixels in spectralis
sigma_lat_bm = 111.13; % ~20 pixels in spectralis

%% Convert all values from micron to pixel
sigma_lat = params.sigma_lat/(header.ScaleX*1000);
sigma_ax = params.sigma_ax/(header.ScaleZ*1000);
distConst = round(params.distconst/(header.ScaleZ*1000));
maxdist = round(maxdist/(header.ScaleZ*1000));
maxdist_bm = round(maxdist_bm/(header.ScaleZ*1000));
isosThresh = round(isosThresh/(header.ScaleZ*1000));
dc_thresh = round(dc_thresh/(header.ScaleZ*1000));
mf_k = round(mf_k/(header.Distance*1000));
sigma_tp_ilm = sigma_tp_ilm/(header.Distance*1000);
sigma_tp_isos = sigma_tp_isos/(header.Distance*1000);
sigma_tp_bm = sigma_tp_bm/(header.Distance*1000);
sigma_lat_ilm = sigma_lat_ilm/(header.ScaleX*1000);
sigma_lat_isos = sigma_lat_isos/(header.ScaleX*1000);
sigma_lat_bm = sigma_lat_bm/(header.ScaleX*1000);

%% Handle nan values on the borders
% The first few or last few columns might be nan, if they are, replace with
% the first non-nan column
fn = find(all(squeeze(isnan(img_vol(:,1,:))),1));
for i = 1:numel(fn)
    j = 1;
    while all(squeeze(isnan(img_vol(:,j,fn(i)))),1)
        j = j + 1;
    end
    img_vol(:,1:j-1,fn(i)) = repmat(img_vol(:,j,fn(i)),[1 j-1]);
end
fn = find(all(squeeze(isnan(img_vol(:,end,:))),1));
for i = 1:numel(fn)
    j = 0;
    while all(squeeze(isnan(img_vol(:,end-j,fn(i)))),1)
        j = j + 1;
    end
    img_vol(:,end-(j-1):end,fn(i)) = repmat(img_vol(:,end-j,fn(i)),[1 j]);
end

% Change nans at top of image to replicate first non-nan value
fn = squeeze(isnan(img_vol(1,:,:)));
if size(img_vol,3) == 1
    fn = fn';
end
fns = double(fn); j = 2;
while any(fn(:)) && j <= size(img_vol,1)
    fn = squeeze(isnan(img_vol(j,:,:)));
    fns = fns + double(fn);
    j = j + 1;
end
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        img_vol(1:fns(i,j),i,j) = img_vol(fns(i,j)+1,i,j);
    end
end

% Change zeros at top of image to replicate first non-nan value
fn = squeeze(img_vol(1,:,:)==0);
if size(img_vol,3) == 1
    fn = fn';
end
fns = double(fn); j = 2;
while any(all(fn)) && j <= size(img_vol,1)
    fn = squeeze(img_vol(j,:,:)==0);
    fns = fns + double(fn);
    j = j + 1;
end
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        img_vol(1:fns(i,j),i,j) = img_vol(fns(i,j)+1,i,j);
    end
end

% Change zeros at bottom of image to replicate last non-nan value
fn = squeeze(img_vol(end,:,:)==0);
if size(img_vol,3) == 1
    fn = fn';
end
fns = double(fn); j = size(img_vol,1)-1;
while any(all(fn)) && j >= 1
    fn = squeeze(img_vol(j,:,:)==0);
    fns = fns + double(fn);
    j = j - 1;
end
for i = 1:size(img_vol,2)
    for j = 1:size(img_vol,3)
        img_vol(end-fns(i,j)+1:end,i,j) = img_vol(end-fns(i,j),i,j);
    end
end

if doplots
    h1 = figure;
    imagesc(img_vol(:,:,sl)), colormap gray
end

%% Pre-processing

% Gaussian filter data
blurImg = calculateFeatures2D({'anigauss',{sigma_ax,sigma_lat,0,0,[]}},img_vol);

if doplots
    figure,imagesc(blurImg(:,:,sl)), colormap gray
end

% Gradient of the image in the y-direction
grad = -imfilter(blurImg,fspecial('sobel'),'replicate');

if doplots
    figure,imagesc(grad(:,:,sl)), colormap gray
end

%% Find ILM and ISOS boundaries

grad_o = grad;

% Largest gradient, this is either the top layer or the center layer
[max1vals, max1pos] = max(grad,[],1);
max1pos = squeeze(max1pos);
if isvector(max1pos)
    max1pos = max1pos';
    max1vals = max1vals';
end
grad2 = grad;

% Find the second largest gradient to the max gradient at a distance of at
% least distConst away but not more than maxdist away
for i = 1:size(grad,2)
    for j = 1:size(grad,3)
        % min distance
        dc = distConst;
        if (max1pos(i,j)-distConst) < 1
            dc = max1pos(i,j) - 1;
        elseif (max1pos(i,j)+distConst) > size(grad,1)
            dc = size(grad,1) - max1pos(i,j);
        end        
        grad((max1pos(i,j)-dc):(max1pos(i,j)+dc),i,j) = 0;
        
        % max distance
        if (max1pos(i,j)-maxdist) > 0
            grad(1:(max1pos(i,j)-maxdist),i,j) = 0;
        end
        if (max1pos(i,j)+maxdist) <= size(grad,1)
            grad((max1pos(i,j)+maxdist):end,i,j) = 0;
        end
    end
end
[max2vals,max2pos] = max(grad,[],1);
max2pos = squeeze(max2pos);
if isvector(max2pos)
    max2pos = max2pos';
    max2vals = max2vals';
end

ilm = min(max1pos,max2pos);
isos = max(max1pos,max2pos);

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),ilm(:,sl),'r.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),isos(:,sl),'b.','markersize',5)
end

% Gradient at these values
[X Z] = ndgrid(1:size(isos,1),1:size(isos,2));
isos_vals = reshape(grad2(sub2ind(size(grad),isos(:),X(:),Z(:))),size(isos));
ilm_vals = reshape(grad2(sub2ind(size(grad),ilm(:),X(:),Z(:))),size(isos));

%% Find BM boundary

grad = grad_o;

% BM is largest negative gradient below the ISOS
for i = 1:size(grad,2)
    for j = 1:size(grad,3)  
        grad(1:isos(i,j)+isosThresh,i,j) = 0;
        
        if (isos(i,j)+maxdist_bm) <= size(grad,1)
            grad((isos(i,j)+maxdist_bm):end,i,j) = 0;
        end
    end
end

[botvals,bot] = min(grad,[],1); 
bm = squeeze(bot); botvals = squeeze(botvals);
if isvector(bm)
    bm = bm';
end

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),bm(:,sl),'g.','markersize',5)
end

%% Detect outliers

if bsc_indep
    % Detect outliers from total retina thickness
    th = bm-ilm;
    th_med = medfilt2(th,[1 mf_k],'symmetric');
    bpt = abs(th-th_med)>dc_thresh;
else
    % Median filter surfaces to detect outliers
    ilm_med = medfilt2(ilm,[1 mf_k],'symmetric');
    isos_med = medfilt2(isos,[1 mf_k],'symmetric');
    bm_med = medfilt2(bm,[1 mf_k],'symmetric');
    bpt = (abs(ilm-ilm_med) > dc_thresh) | (abs(isos-isos_med) > dc_thresh) ...
          | (abs(bm-bm_med) > dc_thresh);
end

% Fill in outlier points with closest point
ilm(bpt) = nan;
isos(bpt) = nan;
bm(bpt) = nan;
nbpt = 0;
if any(bpt(:))
%     disp([num2str(sum(bpt(:))) ' outliers found in layer fit'])
    nbpt = sum(bpt(:));
    ilm = inpaint_nans(ilm);
    isos = inpaint_nans(isos);
    bm = inpaint_nans(bm);
end

%% Get final boundaries by smoothing

% Finally smooth surfaces
if ~bsc_indep
    ilm = imfilter(ilm,fspecial('gaussian',[1 2*round(3*sigma_tp_ilm)+1],...
                                sigma_tp_ilm),'symmetric');
    isos = imfilter(isos,fspecial('gaussian',[1 2*round(3*sigma_tp_isos)+1],...
                                  sigma_tp_isos),'symmetric');
    bm = imfilter(bm,fspecial('gaussian',[1 2*round(3*sigma_tp_bm)+1],...
                              sigma_tp_bm),'symmetric');
end

ilm = imfilter(ilm,fspecial('gaussian',[2*round(3*sigma_lat_ilm)+1 1],...
                            sigma_lat_ilm),'symmetric');
isos = imfilter(isos,fspecial('gaussian',[2*round(3*sigma_lat_isos)+1 1],...
                              sigma_lat_isos),'symmetric');
bm = imfilter(bm,fspecial('gaussian',[2*round(3*sigma_lat_bm)+1 1],...
                          sigma_lat_bm),'symmetric');     

if doplots
    figure(h1),hold on,plot(1:size(img_vol,2),ilm(:,sl),'m.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),isos(:,sl),'c.','markersize',5)
    figure(h1),hold on,plot(1:size(img_vol,2),bm(:,sl),'y.','markersize',5)
end


% Retina mask label image
retinaMask = zeros(size(img_vol),'uint8');
ilm(ilm<1) = 1;
ilm(ilm>size(img_vol,1)) = size(img_vol,1);
isos(isos<1) = 1;
isos(isos>size(img_vol,1)) = size(img_vol,1);
bm(bm<1) = 1;
bm(bm>size(img_vol,1)) = size(img_vol,1);
for i = 1:size(img_vol,2)
    for j = 1:size(grad,3)
        retinaMask(round(ilm(i,j)):round(isos(i,j)-1),i,j) = 1;
        retinaMask(round(isos(i,j)):round(bm(i,j)),i,j) = 2;
    end
end

boundaries = cat(3,ilm,isos,bm);

% Shifts (we shift so that the bottom layer is in the center of the image)
% shifts = bm' - mean(bm) - (round(size(image,1)/2)-mean(bm));
shifts = bsxfun(@minus,bm,mean(bm,1)+(round(size(img_vol,1)/2)-mean(bm,1)));
