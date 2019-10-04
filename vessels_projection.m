function [ J ] = vessels_projection( projection )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


RGB(:,:,1)=projection;
RGB(:,:,2)=projection;
RGB(:,:,3)=projection;
RGB=double(RGB);
% load('G:\code oct\marzieh\1076dl')
% RGB=imread('2043_r.jpg');
RGB=imresize(RGB,[256,256]);
%%%% tariktarin noghteh mahali
% [y_mm,x_mm]=find(RGB==min(min(RGB(110:160,85:157))));
% [aaa,bbb]=find((x_mm<85|x_mm>157));
% x_mm(aaa,:)=[];
% y_mm(aaa,:)=[];
% [aa,bb]=find((y_mm<110|y_mm>160));
% x_mm(aa,:)=[];
% y_mm(aa,:)=[];
% x_m=mean(x_mm);
% y_m=mean(y_mm);
%  x_m1=128;
%  y_m1=128;
% RGB=imread('tr4.tif');

J= imadjust(RGB(:,:,2),[0.1 0.8],[]);
%   J= imadjust(RGB(:,:,2),[0 .5],[0 .9]);
% J= imadjust(RGB(:,:,2),stretchlim(RGB(:,:,2)),[]);
 RGB(:,:,2)=J;
% figure(1)
% imshow(RGB)
mask=msk(RGB);
% figure
% imshow(RGB)
% threshold = graythresh(IG);
% imagen =~im2bw(IG,.015);
% imagen=~imagen;
% figure
% imshow(imagen)
IG = rgb2gray(RGB);
IG=double(IG);
% figure(2)
% imshow(uint8(IG))
str=std(std(IG));
M=max(max(IG));
m=.4*M;
LUV = colorspace('rgb->luv',RGB+1);
L(:,:)= LUV(:,:,1);
    U (:,:)= LUV(:,:,2);
    V (:,:)= LUV(:,:,3);
    %**************************8
    CL = fdct_usfft(L,1);
    CU = fdct_usfft(U,1);
    CV = fdct_usfft(V,1);
    %***********************
  al=CL{1,1}{1,1};
  au=CU{1,1}{1,1};
  av=CV{1,1}{1,1};
  [al1,au1,av1]= colorcontrast3(al,au,av,m,str);
  CLmod{1,1}{1,1}=al1;
  CUmod{1,1}{1,1}=au1;
  CVmod{1,1}{1,1}=av1;
  for i=1:32
 bl=CL{1,2}{1,i};
 bu=CU{1,2}{1,i};
 bv=CV{1,2}{1,i};
 
[bl1,bu1,bv1]= colorcontrast3(bl,bu,bv,m,str);
  CLmod{1,2}{1,i}=bl1;
  CUmod{1,2}{1,i}=bu1;
  CVmod{1,2}{1,i}=bv1;
  end
  for i=1:32
 cl=CL{1,3}{1,i};
 cu=CU{1,3}{1,i};
 cv=CV{1,3}{1,i};
 
[cl1,cu1,cv1]= colorcontrast3(cl,cu,cv,m,str);
  CLmod{1,3}{1,i}=cl1;
  CUmod{1,3}{1,i}=cu1;
  CVmod{1,3}{1,i}=cv1;
  end
  for i=1:64
 dl=CL{1,4}{1,i};
 du=CU{1,4}{1,i};
 dv=CV{1,4}{1,i};
 
[dl1,du1,dv1]= colorcontrast3(dl,du,dv,m,str);
  CLmod{1,4}{1,i}=dl1;
  CUmod{1,4}{1,i}=du1;
  CVmod{1,4}{1,i}=dv1;
  end
%   for i=1:64
%  dl=CL{1,5}{1,i};
%  du=CU{1,5}{1,i};
%  dv=CV{1,5}{1,i};
%  
%   [dl1,du1,dv1]= colorcontrast2(dl,du,dv,m,str);
%   CLmod{1,5}{1,i}=dl1;
%   CUmod{1,5}{1,i}=du1;
%   CVmod{1,5}{1,i}=dv1;
%   end
  el=CL{1,5}{1,1};
  eu=CU{1,5}{1,1};
  ev=CV{1,5}{1,1};
  [el1,eu1,ev1]= colorcontrast3(el,eu,ev,m,str);
  CLmod{1,5}{1,1}=el1;
  CUmod{1,5}{1,1}=eu1;
  CVmod{1,5}{1,1}=ev1;
  %***************************************
  Lmod = ifdct_usfft(CLmod,1);
  Umod = ifdct_usfft(CUmod,1); 
  Vmod = ifdct_usfft(CVmod,1);
    %************************************************
    LUVmod(:,:,1)=Lmod;
    LUVmod(:,:,2)=Umod;
    LUVmod(:,:,3)=Vmod;
      %********************************************************
RGBmod2 = colorspace('rgb<-luv',LUVmod);
gray=RGBmod2(:,:,2);
se = strel('disk',1);
erodedmask = imerode(mask,se);
% figure(3), imshow(erodedmask)
gray=gray.*erodedmask;
% figure(4)
% imshow(gray)
norm=abs(max(max(gray(:))));
X22=(gray/norm)*256;
rX22=X22;
% figure(5)
% imshow(uint8(rX22))
 sa=.3;
    [M,N] = size(rX22);
x = [-4: 4];
tmp1 = exp(-(x.*x)/(2*sa*sa)); 
tmp1 = max(tmp1)-tmp1; 
ht1 = repmat(tmp1,[9 1]); 
sht1 = sum(ht1(:));
mean1 = sht1/(9*9);
ht1 = ht1 - mean1;
ht1 = ht1/sht1;

h{1} = zeros(11,11);
for i = 1:9
    for j = 1:9
        h{1}(i+3,j+1) = ht1(i,j);
    end
end

for k=1:11
    ag = 15*k;
    h{k+1} = imrotate(h{1},ag,'bicubic');
    h{k+1} = wkeep(h{k+1},size(h{1}));
end

for k=1:12
    R{k} = conv2(rX22, h{k}, 'same');
end

rt = zeros(M,N);
for i=1:M
    for j=1:N
        ER = [R{1}(i,j), R{2}(i,j), R{3}(i,j), R{4}(i,j), R{5}(i,j), R{6}(i,j),... 
                R{7}(i,j), R{8}(i,j), R{9}(i,j), R{10}(i,j), R{11}(i,j), R{12}(i,j)];
        rt(i,j) = max(ER);
    end
end

rmin = abs(min(rt(:)));
for m = 1:M
    for n = 1:N
        rt(m,n) = rt(m,n) + rmin;
    end
end
rmax = max(max(rt));
for m = 1:M
    for n = 1:N
        rt(m,n) = round(rt(m,n)*255/rmax);
    end
end
rt=imresize(rt,[256,256]);
%  figure(6)
%  imshow(uint8(rt))
 %%%%%%%%%%%%%%%blood vessel
 C3 = fdct_usfft(rt,1);
  a=C3{1,1}{1,1};
 D=a;
 a1{1,1}{1,1}=edgenhance(D);
 for i=1:32
 b=C3{1,2}{1,i};
 D=b;
 b1{1,2}{1,i}=edgenhance(D);
 end
  for i=1:32
 c=C3{1,3}{1,i};
  D=c;
 c1{1,3}{1,i}=mim_eslah_zarayeb(D,2);
  end
  for i=1:64
 d=C3{1,4}{1,i};
   D=d;
 d1{1,4}{1,i}=mim_eslah_zarayeb(D,2);
  end
 e=C3{1,5}{1,1};
 D=e;
 e1{1,1}{1,1}=mim_eslah_zarayeb(D,2);
 %***********************
 C4={a1{1,1},b1{1,2},c1{1,3},d1{1,4},e1{1,1}};
 X3 = ifdct_usfft(C4,1);
%  figure(7)
%  imshow(uint8(X3))
  for i=1:256
     for j=1:256
         if X3(i,j)>0
             X4(i,j)=X3(i,j);
         else
             X4(i,j)=0;
         end
     end
 end
NO=max(max(X4(:)));
X4=(X4/NO)*256;
% figure(8)
%  imshow(uint8(X4))

%  NO2=max(max(X4(:)));
%  
%  X5=X4-(1/13)*NO2;
%  figure
%  imshow(uint8(X5))
%  maa=mean2(X5(:));
%  for i=1:256
%      for j=1:256
%          if X5(i,j)>maa
%              X6(i,j)=1;
%          else
%              X6(i,j)=0;
%          end
%      end
%  end
% figure
% imshow(X6)
maa=sum(X4(:));
 ma1=1.5*maa/(256*256);
 for i=1:256
     for j=1:256
         if X4(i,j)>ma1
             X6(i,j)=1;
         else
             X6(i,j)=0;
         end
     end
 end
% figure(9)
% imshow(X6)
X7=X6;
%  X7 = bwmorph(X6,'bridge',1);
% figure
%  imshow(X7)
% se1= strel('disk',1);
% X7 = imdilate(X6,se1);
%  figure
% imshow(X7)
[Label,Num] = bwlabel(X7);
Lmtx = zeros(Num+1,1);
for i=1:256
    for j=1:256
        Lmtx(double(Label(i,j))+1) = Lmtx(double(Label(i,j))+1) + 1;
    end
end
sLmtx = sort(Lmtx);

cp = 100;
for i=1:256
    for j=1:256
        if (Lmtx(double(Label(i,j)+1)) > cp) & (Lmtx(double(Label(i,j)+1)) ~= sLmtx(Num+1,1))
            J(i,j) = 1;
        else
            J(i,j) = 0;
        end
    end
end
% save v682f J
% J=imresize(J,[622,667]);
% J=imresize(J,[622,667]);


%  SE = strel('line',13,120);
% J= imopen(J,SE);


 figure(10); imshow(~J,[]);
end

