function thicknessMapAMPAO2(boundaryPoints,header,circ_rad,showgrid,handles)

if nargin < 4
    showgrid = false;
end

% Circle radii for rings of ETDRS grid
if nargin < 3
    circ_rad = [0.5 1.5 2.5]*1000;
end

% Compute thickness as distance between boundary points
thicknessVals = squeeze(diff(boundaryPoints,1,3))*header.ScaleZ*1000;

% Rotate and flip to align with fundus
thicknessVals = flipdim(permute(thicknessVals,[2 1 3]),1);

% Flip if left eye for consistency of 'nasal' to the right side
if ~strncmp(header.ScanPosition,'OD',2)
    thicknessVals = flipdim(thicknessVals,2);
end

scaleX = header.ScaleX*1000;
scaleY = header.Distance*1000;

% Resize to 512 x 512 for improved resolution
newsize = [512 512];
scf = [newsize(1)/size(thicknessVals,1) newsize(2)/size(thicknessVals,2)];
scaleX = scaleX/scf(2);
scaleY = scaleY/scf(1);
thicknessVals = imresize(thicknessVals,newsize);

% figure,imagesc([scaleX/2 scaleX*(size(thicknessVals,2)-1)+scaleX/2],...
%                [scaleY/2 scaleY*(size(thicknessVals,1)-1)+scaleY/2],sum(thicknessVals,3))

%% Find fovea center point (looking at minimum total retina thickness)

% Retina thickness
th = sum(thicknessVals,3);

% Only look in the center part of the volume
cr = round([0.35 0.65]*size(thicknessVals,2));
cr2 = round([0.35 0.65]*size(thicknessVals,1));
thi = th(cr2(1):cr2(2),cr(1):cr(2));

% Minimum thickness average in 0.3 mm radius circle
r = 300;
nx = ceil(r/scaleX); ny = ceil(r/scaleY);
py = (-ny:ny)*scaleY; px = (-nx:nx)*scaleX;
[X,Y] = meshgrid(px,py);
R = sqrt(X.^2+Y.^2);
m = R < r;
thi = imfilter(thi,double(m)/sum(m(:)),'replicate');

[~,ptm] = min(thi(:));
[ptm1, ptm2] = ind2sub(size(thi),ptm);

% Fovea center (adding back the crop)
f_cen = [ptm1 + cr2(1)-1, ptm2 + cr(1)-1];

% figure,imagesc(sum(thicknessVals,3)), hold on, plot(f_cen(2),f_cen(1),'g*')

%% Extract average in each ETDRS region

xpts = 1:size(thicknessVals,2);
ypts = 1:size(thicknessVals,1);
[X,Y] = meshgrid(xpts,ypts);

thicknessVals = permute(thicknessVals,[3 1 2]);
thicknessVals_sec = zeros(size(thicknessVals,1),10);

X_ij = (X - squeeze(f_cen(2)))*scaleX;
Y_ij = (Y - squeeze(f_cen(1)))*scaleY;
R_ij = sqrt(X_ij.^2 + Y_ij.^2);

% Center area - radius 0.5
p_center = R_ij < circ_rad(1);
thicknessVals_sec(:,1) = mean(thicknessVals(:,p_center),2);

% Inner ring - between radius 0.5 and 1.5
p_inner = (R_ij < circ_rad(2)) & (R_ij >= circ_rad(1));
p_inner_superior = p_inner & (Y_ij < 0) & (abs(Y_ij) > abs(X_ij));
p_inner_inferior = p_inner & (Y_ij >= 0) & (abs(Y_ij) > abs(X_ij));
p_inner_nasal = p_inner & (X_ij > 0) & (abs(Y_ij) <= abs(X_ij));
p_inner_temporal = p_inner & (X_ij <= 0) & (abs(Y_ij) <= abs(X_ij));

thicknessVals_sec(:,2) = mean(thicknessVals(:,p_inner_superior),2);
thicknessVals_sec(:,3) = mean(thicknessVals(:,p_inner_inferior),2);
thicknessVals_sec(:,4) = mean(thicknessVals(:,p_inner_nasal),2);
thicknessVals_sec(:,5) = mean(thicknessVals(:,p_inner_temporal),2);

% Outer ring - between radius 1.5 and 2.5
p_outer = (R_ij < circ_rad(3)) & (R_ij > circ_rad(2));
p_outer_superior = p_outer & (Y_ij < 0) & (abs(Y_ij) > abs(X_ij));
p_outer_inferior = p_outer & (Y_ij >= 0) & (abs(Y_ij) > abs(X_ij));
p_outer_nasal = p_outer & (X_ij > 0) & (abs(Y_ij) <= abs(X_ij));
p_outer_temporal = p_outer & (X_ij <= 0) & (abs(Y_ij) <= abs(X_ij));

thicknessVals_sec(:,6) = mean(thicknessVals(:,p_outer_superior),2);
thicknessVals_sec(:,7) = mean(thicknessVals(:,p_outer_inferior),2);
thicknessVals_sec(:,8) = mean(thicknessVals(:,p_outer_nasal),2);
thicknessVals_sec(:,9) = mean(thicknessVals(:,p_outer_temporal),2);

p_outer_rim = R_ij >= circ_rad(3);

% Full macula within 5x5 cube centered at fovea
p_rect = abs(X_ij) < 2500 & abs(Y_ij) < 2500;
thicknessVals_sec(:,10) = mean(thicknessVals(:,p_rect),2);

stats = thicknessVals_sec;
% uisave( 'stats','ThiknessMapData')

%% Figure

if showgrid
    axes(handles.axes2);
    %     imagesc(handles.data.slo)
    %     colormap gray
    %     hold on
    
    sector_th_img = nan(size(R_ij));
    sector_th_img(p_center) = sum(thicknessVals_sec(:,1),1);
    sector_th_img(p_inner_superior) = sum(thicknessVals_sec(:,2),1);
    sector_th_img(p_inner_inferior) = sum(thicknessVals_sec(:,3),1);
    sector_th_img(p_inner_nasal) = sum(thicknessVals_sec(:,4),1);
    sector_th_img(p_inner_temporal) = sum(thicknessVals_sec(:,5),1);
    sector_th_img(p_outer_superior) = sum(thicknessVals_sec(:,6),1);
    sector_th_img(p_outer_inferior) = sum(thicknessVals_sec(:,7),1);
    sector_th_img(p_outer_nasal) = sum(thicknessVals_sec(:,8),1);
    sector_th_img(p_outer_temporal) = sum(thicknessVals_sec(:,9),1);
    
    % Overlay white boundary
    imagesc(sector_th_img);
    %     imgg = imagesc(sector_th_img);
    %     imgg.AlphaData = .5;
    colormap winter
    
    t = linspace(0,2*pi,500);
    hold on;
    % Three concentric circles - note color is [1 1 0.99] due to a bug in
    %   print -depsc2
    plot(circ_rad(3)/scaleX*cos(t) + f_cen(2),...
        circ_rad(3)/scaleY*sin(t) + f_cen(1),...
        'Color',[1 1 0.99],'LineWidth',2)
    plot(circ_rad(2)/scaleX*cos(t) + f_cen(2),...
        circ_rad(2)/scaleY*sin(t) + f_cen(1),...
        'Color',[1 1 0.99],'LineWidth',2)
    plot(circ_rad(1)/scaleX*cos(t) + f_cen(2),...
        circ_rad(1)/scaleY*sin(t) + f_cen(1),...
        'Color',[1 1 0.99],'LineWidth',2)
    % Radial lines
    t = [pi/4,3*pi/4,5*pi/4,7*pi/4];
    plot([circ_rad(3)/scaleX*cos(t) + f_cen(2);...
        circ_rad(1)/scaleX*cos(t) + f_cen(2)],...
        [circ_rad(3)/scaleY*sin(t) + f_cen(1);...
        circ_rad(1)/scaleY*sin(t) + f_cen(1)],...
        'Color',[1 1 0.99],'LineWidth',2);
    
    axis square
    %     axis off
    colorbar
    title({'\fontsize{10}Thickness values' ; 'by ETDRS grid sector'})
    %     title(sprintf('Thickness values by ETDRS grid sector'));
    % Flip if left eye for consistency of 'nasal' to the right side

end
