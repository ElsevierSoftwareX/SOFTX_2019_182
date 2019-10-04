function [ layer1,layer2,layer3, layer4 ] = find_layer(imagei,image,thereshold )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
X1=image;
level1 =graythresh((X1));
   level=thereshold*level1;
   [Gx,Gy]=imgradient(imagei);
   BW2=im2bw((Gx),level);
    BW = im2bw((X1),level);
X2=abs(1-X1);
    X3=adapthisteq(X2);
    se = strel('disk',5);
        se1 = strel('disk',3);
        afterFilling=imfill(BW,'holes');
        BWareaopen = bwareaopen(afterFilling, 700);
        afterClosing=imclose(BWareaopen,se);
        afterOpening = imopen(afterClosing,se1);
   
  edge=1:768;
%   clear r1;
    for j=1:768
        clear r1
        [r,c]=find(afterOpening(:,j)==1);
        [r2,c2]=find(BW2(:,j)==1);
        a=size(r,1);
        a2=size(r2,1);
%         r1(1)=r(1);
%         for ii=1:a-1
%             [r3,c3]=find(abs(r(ii)-r(ii+1))>=2);
%         end  
%  r1(2:end-1,1)=r(r3,1);
 
%         r1(end+1)=r(a);
%         a2=size(r1,1);
TFS=r-circshift(r,1);

[r4,c4]=find(TFS>=20);
% if size(r4)~=1;
%     r4=max(r4);
% end
if a2~=0
    r1(1)=r2(1);
else
    r1(1)=NaN;
end
if size(r4,1)==1 ;
    r1(1)=r(1);
     r1(2)=r(r4-1);
     r1(3)=r(r4);
    r1(4)=r(end);
else
    r1(1)=NaN;
     r1(2)=NaN;
   r1(3)=NaN;
    r1(4)=NaN;
end
if size(r4,1)==1;
           layer1(j)=r1(1);
            layer2(j)=r1(2);
         layer3(j)=r1(3);
         layer4(j)=r1(4);
       else
           layer1(j)=r1(1);
            layer2(j)=NaN;
            layer3(j)=NaN;
         layer4(j)=r1(4);
end
end

