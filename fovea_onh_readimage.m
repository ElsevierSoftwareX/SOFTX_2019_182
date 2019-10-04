% close all
% clear all
% [FileName,PathName] = uigetfile('*.png','Select xml file');
% vab= imread(FileName);
function[center00,radius0,x0, y0,x00,y00]=fovea_onh_readimage(SLO)
% load SLO_Dr_Alamzadeh_D_157
% SLO=uint8(SLO_Dr_Alamzadeh_D_157);
% SLO = fliplr(SLO);
vab=SLO;
% vab = imadjust(vab);
% vab= vab(80:768-80,1:768-50);
% vab= vab(10:768-100,:);
% vab= vab(50:768-50,79:768);
 ad=vab;
%  ad=rgb2gray(vab);
a=adapthisteq(vab);
side=12;
SE=strel('octagon',side);
b=imclose(a,SE);
c=imopen(b,SE);
A=0.43;
cc=im2bw(c,A);
d=~cc;

[label1 NUM]=bwlabeln(d);
stat1 = regionprops(label1,'all');
D=500;
pixel=((3.14*D*D)/(8*7.33));

for m=1:NUM
  if ((pixel)/4<stat1(m).Area && stat1(m).Area<(pixel*5)/4)
      A(m)=1;
  
  else
      A(m)=0;
  end
end
   
    
for m=1:NUM
 if((stat1(m).Area)/(stat1(m).BoundingBox(3)*stat1(m).BoundingBox(4))<0.5)
     B(m)=0;
 else
     B(m)=1;
     
 end
 
end

for m=1:NUM
    if(stat1(m).BoundingBox(3)*7<stat1(m).BoundingBox(4) ||7*stat1(m).BoundingBox(4)<stat1(m).BoundingBox(3))
       C(m)=0;
    else
     C(m)=1 ;
    end
end

chan=d;
for m=1:NUM
    if(A(m)==1 && B(m)==1 &&C(m)==1)
        chan(stat1(m).PixelIdxList)=1;
    else
      chan(stat1(m).PixelIdxList)=0;
    end
end


[label2 NUM1]= bwlabel(chan);
stat2 = regionprops(label2,'all');



[maxValue,index] = max([stat2.Area]);

x0 = stat2(index).Centroid(1);
y0 = stat2(index).Centroid(2);
centers = [x0 y0];
radius0 = stat2(index).EquivDiameter/2;

 %imshow(vab,[]);
% title('Localised optic nerve head');
% hold on
% viscircles(centers,radius0);
% plot(x0,y0,'+');
% hold on    
%%%%%%%%%%%%%%%%%%%
%finding fovea
if(x0<400 && 360<x0)
figure;
 subplot(1,2,1);

imshow(vab,[]);
title('Original slo image');
 
% subplot(1,2,2);
%  
%  imshow(vab,[]);
% title('Localised optic nerve head');
% hold on
% viscircles(centers,radius0);
% plot(x0,y0,'+');
% hold on      
    
end 


if(x0<360)
  a1 =pi/6+pi/2 ;  
   a2 = -pi/6+pi/2;
   t = linspace(a1,a2);
   %first circle
 x=5*radius0*sin(t)+x0;
 y=5*radius0*cos(t)+y0;
 %secend circle
 x1=10*radius0*sin(t)+x0;
 y1=10*radius0*cos(t)+y0;
 
  figure;
 subplot(2,2,1);

imshow(vab,[]);
title('Original slo image');
 
subplot(2,2,2);
 
 imshow(vab,[]);
title('Localised optic nerve head');
hold on
viscircles(centers,radius0);
plot(x0,y0,'+');
hold on  
subplot(2,2,3);
imshow(vab,[]);
title('Localised macula');
hold on
viscircles(centers,radius0);
plot(x0,y0,'+');
hold on

plot(x,y);
axis equal
 hold on
  plot(x1,y1);
axis equal
hold on

line([x(1) x1(1)],[y(1) y1(1)]);
line([x(100) x1(100)],[y(100) y1(100)]);

 
 q=true(size(vab));
 [x,y]=size(q);
 for i=1:x
    for j=1:y
        if (5*radius0<sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)) && sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0))<10*radius0 && (j-y0)/(i-x0)<0.5 && (j-y0)/(i-x0)>-0.5 )
           q(j,i)=0;
        else 
          q(j,i)=1;  
        end
    end
end
 
 cc1=~q;
%cenrtroid of fovea
[label3 NUM3]=bwlabeln(cc1);
stat3= regionprops(label3,'all');
[maxValue,index] = max([stat3.Area]);
x00 = stat3(index).Centroid(1);
y00 = stat3(index).Centroid(2);
center00=[x00 y00+39];


subplot(2,2,4);
 imshow(vab,[]);
 title('Localized fovea');
 hold on
viscircles(center00,radius0);
plot(x00,y00+39,'+');



 end

if(x0>400)
  
   a1 =5*pi/6+pi/2 ; 
   a2 = 7*pi/6+pi/2;
   t = linspace(a1,a2);
   %first circle
 x=5*radius0*sin(t)+x0;
 y=5*radius0*cos(t)+y0;

 %secend circle
 x1=10*radius0*sin(t)+x0;
 y1=10*radius0*cos(t)+y0;
 
 figure;
 subplot(2,2,1);
imshow(vab,[]);
title('Original slo image');
 
subplot(2,2,2);
imshow(vab,[]);
title('Localized optic nerve head');
hold on
viscircles(centers,radius0);
plot(x0,y0,'+');

 subplot(2,2,3);
imshow(vab,[]);
title('Localized macula');
hold on
viscircles(centers,radius0);
plot(x0,y0,'+');
hold on

 plot(x,y);
axis equal
 hold on
  plot(x1,y1);
axis equal
hold on

line([x(1) x1(1)],[y(1) y1(1)]);
line([x(100) x1(100)],[y(100) y1(100)]);

 q=true(size(ad));
 [x,y]=size(q);
 for i=1:x
    for j=1:y
        if (5*radius0<sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0)) && sqrt((i-x0)*(i-x0)+(j-y0)*(j-y0))<10*radius0 && (j-y0)/(i-x0)<0.5 && (j-y0)/(i-x0)>-0.5 )
           q(j,i)=0;
        else 
          q(j,i)=1;  
        end
    end
end
 
 cc1=~q;

 
%cenrtroid of fovea
[label3 NUM3]=bwlabeln(cc1);
stat3= regionprops(label3,'all');
[maxValue,index] = max([stat3.Area]);
x00 = stat3(index).Centroid(1);
y00 = stat3(index).Centroid(2);
center00=[x00 y00+30];

subplot(2,2,4);
 imshow(vab,[]);
 title('Localised fovea');
 hold on
viscircles(center00,radius0);
plot(x00,y00+30,'+');

FoveaLoc = [x00,y00+39];
ONHLoc = [x0,y0];    
uisave({'FoveaLoc','ONHLoc'})


end
