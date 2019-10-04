function [ xdmean,ydmean,phi ] = drlse_centerdisk( rgb_shift,x_d,y_d )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
  I=rgb_shift;
 m= 2*ones(size(I,1),size(I,2));
  m(y_d-25: y_d+20,x_d-25:x_d+40)=-2;
 figure(1), imshow(I), hold on, contour(m, [0,0], 'r','LineWidth',2);
Img=imresize(double(I(:,:,1)),1/4);
%% parameter setting
timestep=1;  % time step
mu=0.2/timestep;  % coefficient of the distance regularization term R(phi)
iter_inner=20;
iter_outer=30;%250
lambda=5; % coefficient of the weighted length term L(phi)
alfa=-3;  % coefficient of the weighted area term A(phi)
epsilon=1.5; % papramater that specifies the width of the DiracDelta function

sigma=.8;    % scale parameter in Gaussian kernel
G=fspecial('gaussian',15,sigma); % Caussian kernel
Img_smooth=conv2(Img,G,'same');  % smooth image by Gaussiin convolution
[Ix,Iy]=gradient(Img_smooth);
f=Ix.^2+Iy.^2;
g=1./(1+f);  % edge indicator function.

% initialize LSF as binary step function
phi=imresize(m,1/4);


% figure(1);
% mesh(-phi);   % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% title('Initial level set function');
% view([-80 35]);

% figure(2);
% imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
% title('Initial zero level contour');
% pause(0.5);

potential=2;  
if potential ==1
    potentialFunction = 'single-well';  % use single well potential p1(s)=0.5*(s-1)^2, which is good for region-based model 
elseif potential == 2
    potentialFunction = 'double-well';  % use double-well potential in Eq. (16), which is good for both edge and region based models
else
    potentialFunction = 'double-well';  % default choice of potential function
end  

% start level set evolution
for n=1:iter_outer
    phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);    
    if mod(n,2)==0
        figure(2);
        imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
    end
end

% refine the zero level contour by further level set evolution with alfa=0
alfa=0;
iter_refine = 10;
phi = drlse_edge(phi, g, lambda, mu, alfa, epsilon, timestep, iter_inner, potentialFunction);

finalLSF=phi;
figure(3);
imagesc(Img,[0, 255]); axis off; axis equal; colormap(gray); hold on;  contour(phi, [0,0], 'r');
hold on;  contour(phi, [0,0], 'r');
str=['Final zero level contour, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
title(str);

% figure;
% mesh(-finalLSF); % for a better view, the LSF is displayed upside down
% hold on;  contour(phi, [0,0], 'r','LineWidth',2);
% view([-80 35]);
% str=['Final level set function, ', num2str(iter_outer*iter_inner+iter_refine), ' iterations'];
% title(str);
% axis on;
% [nrow, ncol]=size(Img);
% axis([1 ncol 1 nrow -5 5]);
% set(gca,'ZTick',[-3:1:3]);
% set(gca,'FontSize',14);
 L=imresize(phi,4);
  [y,x]=find(L<=0);
  xdmean=mean(x);
  ydmean=mean(y);
figure(4),
imshow(I); axis off; axis equal; hold on;  contour(imresize(phi,4), [0,0], 'r');
% hold on
% plot(xdmean,ydmean,'*');
end

