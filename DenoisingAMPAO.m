%This code is provided for non-comercial use, corrosponding to the following paper
%"Image Restoration Using Gaussian Mixture Models With Spatially Constrained Patch Clustering" 
%Milad Niknejad, Hossein Rabbani, Massoud Babaie-Zadeh, Image Processing, IEEE Transactions on  Vol 24 ,  No. 11, PP. 3624 - 3636 Nov. 2015  

% clc;
% clear all;
% close all;
% %---------------------------------------- Topcon---------------------------
function output = DenoisingAMPAO(handles)

% if FileName(end-2:end) == 'vol'
%     [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (strcat(PathName,FileName));
%     %     OCT3D = BScans;
%     image1 = BScans(:,:,scan_num);
    
if strcmp(handles.data.format , 'vol')
    image = double(handles.data.BScan).^.25;
else 
% image = double(information.oct);
% image = double(handles.data.BScan(:,:,handles.ScanNum_denoise));
    image = double(handles.data.BScan);
end

mn1_MVGauss = image;
% image_number_all=size(mn1_MVGauss,3);
%--------------------------------------------------------------------------
    
ju=1;
refju=10; %jump size for reference patches
mxco=.00005;
N=30; %size of clustering window
s=10; % size of patch: (s-1)x(s-1)
mnc=40;%maximum number of cluster dimention
nmiter=20; %number of iterations of EM-like

% scan_num=[1:2];

%------------- Get Slice Numbers
prompt={'Enter value \range of slices to be denoised'};
name = 'Scan Number(s)';
defaultans = {'1'};
options.Interpreter = 'none';
answer = inputdlg(prompt,name,[1 40],defaultans,options);
scan_num=str2num(answer{1,1});
%-------------
uiwait(msgbox('The operation may take a while, Please wait...','Warning','modal'));
 
for count = scan_num
    
    this_scan = count
    imn=mn1_MVGauss(:,:,count); %Gaussianized image
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
    [imfi]=main(ju,refju,N,imn1,s,mnc,iter,sigma_n,imn,agcoe,this_scan);

    imn1=imfi;
    end

    imfi=exp(imfi)-1;
%   figure,imagesc(imfi),axis image; axis off;colormap gray 
    Im_den(:,:,count)=imfi;
end
output.BScan = Im_den;
output.NumBScans = scan_num;
% handles.data.PreProcess = output;
% guidata(hObject,handles);

uisave('output')
% save Im_den