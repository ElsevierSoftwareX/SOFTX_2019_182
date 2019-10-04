%This code is provided for non-comercial use, corrosponding to the following paper
%"Image Restoration Using Gaussian Mixture Models With Spatially Constrained Patch Clustering" 
%Milad Niknejad, Hossein Rabbani, Massoud Babaie-Zadeh, Image Processing, IEEE Transactions on  Vol 24 ,  No. 11, PP. 3624 - 3636 Nov. 2015  

% clc;
% clear all;
% close all;
% %---------------------------------------- Topcon---------------------------
function denoising(handles)

% if FileName(end-2:end) == 'vol'
%     [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (strcat(PathName,FileName));
%     %     OCT3D = BScans;
%     image1 = BScans(:,:,scan_num);
    
 
% image = double(information.oct);
image = double(handles.data.BScans(:,:,handles.ScanNum_denoise));
mn1_MVGauss = image;
% image_number_all=size(mn1_MVGauss,3);
%--------------------------------------------------------------------------
    
ju=1;
refju=10; %jump size for reference patches
mxco=.00005;
N=30; %size of clustering window
s=10; % size of patch: (s-1)x(s-1)
mnc=40;%maximum number of cluster dimention
nmiter=5; %number of iterations of EM-like

% count=1;
% for count=1:image_number_all
    imn=mn1_MVGauss;%(:,:,count); %Gaussianized image
    imn=log(imn+1);
    sigma_n= 3*sqrt(evar(imn)); % estimated variance of noise

    % % % Noise variance estimation using ROI
% %     figure(100),h_im = imagesc(imn);colormap gray
% %     title('Please determine Noise ROI')
% %     e = imellipse;
% %     position = wait(e);
% %     BW_noise = createMask(e,h_im);
% %     sigma_n= 3.5*std (BW_noise(:));
% %     close (figure(100))

%% Gaussianization
% Estimate the parameters of mixture model using EM algorithm
Num_component=5; % number of components in the mixture model
Wsize=3; % size of window for building (Wsize x Wsize) multivariate pdf 
[mu,sigma,coef] = MVEstLapGaussMix_Glob(imn,sigma_n,Num_component,Wsize);

Lcnt=Num_component;

% Gaussianization
% method determines our rule for denoising;
            %if method=0 => without denoising
            %if method=1 => soft threshold            
            % if method=2 => wiener
method=0;
y0= MVGaussianizeDenoise(imn,mu,sigma,coef,sigma_n,Num_component,method,Wsize);
imn2=exp(y0)-1;

%% Denoising with SC-GMM
    %agrregation weights coeficients
    if sigma_n<75
        agcoe=.00006;
    else
        agcoe=.000001;
    end

    iterres=zeros(20,1);
    imn=log(imn2+1);

    imn1=imn;
    sigma1=sigma_n;
    for iter=1:nmiter
    % if mod(iter,3)==0
    %     refju=5;
    % elseif mod(iter,3)==1
    %     refju=6;
    % else
    %     refju=5;
    % end
    [imfi]=main(ju,refju,N,imn1,s,mnc,iter,sigma_n,imn,agcoe);

    imn1=imfi;
    end

    imfi=exp(imfi)-1;
%   figure,imagesc(imfi),axis image; axis off;colormap gray 
    Im_den(:,:,count)=imfi;
% end

% save Im_den