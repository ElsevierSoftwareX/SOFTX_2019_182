function[ output ] = vessel_extraction(rgb_rotate )

SE0=zeros(7);SE1=SE0;SE2=SE0;SE3=SE0;SE4=SE0;SE5=SE0;SE6=SE0;SE7=SE0;SE8=SE0;SE9=SE0;SE10=SE0;SE11=SE0;SE0(4,:)=1;
SE1(4,4)=1;SE1(3,7)=1;SE1(5,1)=1;SE2(4,4)=1;SE2(2,7)=1;SE2(3,6)=1;SE2(6,1)=1;SE2(5,2)=1;
x3=[1 2 3 4 5 6 7]';y3=[7 6 5 4 3 2 1]';xy3=[x3 y3];xy3_ind=sub2ind(size(SE0),xy3(:,1),xy3(:,2));SE3(xy3_ind)=1;
x4=[1 2 4 6 7]';y4=[6 5 4 3 2]';xy4=[x4 y4];xy4_ind=sub2ind(size(SE0),xy4(:,1),xy4(:,2));SE4(xy4_ind)=1;
x5=[1 4 7]';y5=[5 4 3]';xy5=[x5 y5];xy5_ind=sub2ind(size(SE0),xy5(:,1),xy5(:,2));SE5(xy5_ind)=1;
x6=[1 4 7]';y6=[3 4 5]';xy6=[x6 y6];xy6_ind=sub2ind(size(SE0),xy6(:,1),xy6(:,2));SE6(xy6_ind)=1;
x7=[1 2 4 6 7]';y7=[2 3 4 5 6]';xy7=[x7 y7];xy7_ind=sub2ind(size(SE0),xy7(:,1),xy7(:,2));SE7(xy7_ind)=1;
x8=[1 2 3 4 5 6 7]';y8=[1 2 3 4 5 6 7]';xy8=[x8 y8];xy8_ind=sub2ind(size(SE0),xy8(:,1),xy8(:,2));SE8(xy8_ind)=1;
x9=[2 3 4 5 6]';y9=[1 2 4 6 7]';xy9=[x9 y9];xy9_ind=sub2ind(size(SE0),xy9(:,1),xy9(:,2));SE9(xy9_ind)=1;
x10=[3 4 5]';y10=[1 4 7]';xy10=[x10 y10];xy10_ind=sub2ind(size(SE0),xy10(:,1),xy10(:,2));SE10(xy10_ind)=1;
SE11(:,4)=1;
se_siz7(:,:,1)=SE0;se_siz7(:,:,2)=SE1;se_siz7(:,:,3)=SE2;se_siz7(:,:,4)=SE3;se_siz7(:,:,5)=SE4;se_siz7(:,:,6)=SE5;se_siz7(:,:,7)=SE11;
se_siz7(:,:,8)=SE6;se_siz7(:,:,9)=SE7;se_siz7(:,:,10)=SE8;se_siz7(:,:,11)=SE9;se_siz7(:,:,12)=SE10;
% [filename, pathname, filterindex]=uigetfile( ...
%     {'*.jpg','JPEG File (*.jpg)'; ...
%     '*.*','Any Image file (*.*)'}, ...
%     'Pick a microscopic image file');
% bmpin=strcat(pathname,filename);
% fluorangio1 = (load(bmpin));
% %%%%color
 fluorangio=rgb_rotate(:,:,1);
%  fluorangio=fluorangio.projection;
 fluorangio=imresize(fluorangio,[256,256]);
 msk=mask_malihe(fluorangio);
 clear img1
 img1=zeros(size(fluorangio,1)+50,size(fluorangio,2)+50);
img1(26:25+size(fluorangio,1),26:25+size(fluorangio,2))=fluorangio;
 M = medfilt2(img1,[40 40]);
f_norm=zeros(size(img1));f_norm = NORMALIZE_image((double(img1) - double(M)),256);
F_norm=zeros(size(fluorangio));
F_norm=f_norm(26:25+size(fluorangio,1),26:25+size(fluorangio,2));
 
 img0=imcomplement(F_norm);
    max_f=max(max(img0));
min_f=min(min(img0));
for i=1:size(img0,1)
    for j=1:size(img0,2)
   img0(i,j)=256/(max_f-min_f)*(img0(i,j)-min_f);
    end
end
 img = img0;%nonlinear_diffusion( double(img0));
 img=(img./max(img(:)))*256;
 see=strel('disk',1);
img21=imdilate(imerode(img,see),see);
s_c=strel('disk',2);
vec=[1 2 3 4 5 6 7 8];%img=(img/max(img(:)))*256;
clear g1 g2
ii=1;
for i=8%1:size(vec,2)
    s_o=strel('disk',vec(i));
g1=zeros(size(img));g1=imopen(imclose(img,s_c),s_o);
g2(:,:,ii)=img-min(g1,img);ii=ii+1;
end
img2=g2(:,:,1);
[vessel_long,vessel_thin,RESRES,f]= vessel_segmentation_malihe(img,img2);
[result10,result9,ff_2,ff_3,ff_4] = Fillingnew_malihe( bwareaopen(vessel_long,15),f);
f_2_1=zeros(size(ff_2));
f_2_1(:,round(2*(size(ff_2)/3)):end)=ff_2(:,round(2*(size(ff_2)/3)):end);
f_final=f_2_1+ff_3;
f_final=im2bw(f_final,0);
f_final=bwareaopen(f_final,5).*msk;
result10=result10.*msk;
%%
ff_2=im2bw(ff_2,0).*msk;
ff_3=im2bw(ff_3,0).*msk;
ff_4=im2bw(ff_4,0).*msk;
f_final2=f_2_1+ff_4;
f_final2=im2bw(f_final2,0);
f_final2=bwareaopen(f_final2,5).*msk;
%% mahmoodi
ff_44=imresize(ff_4,[401,398]);
ff_440=im2bw(ff_44,.6);

ff_441=bwareaopen(ff_440,9);
result10_1=imresize(result10,[256,256]);
result10_1=im2bw(result10_1,.8);
Vessels_oct=result10_1;
output=imresize(result10_1, [size(rgb_rotate,1),size(rgb_rotate,2)]);

% figure,imshow(result10_1)
% save 'D:\mydata\oct-mokhtary\mokhtari\data_oct\OD-mat\v15340rp.mat'  Vessels_oct
end
% save 'D:\mydata\oct-mokhtary\mokhtari\data base\20\PROJEC VESSEL MACULA R'  Vessels_oct
% imwrite(Vessels_oct,('C:\mydata\oct-mokhtary\mokhtari\DATA_fundus1\685 dr.gif')));