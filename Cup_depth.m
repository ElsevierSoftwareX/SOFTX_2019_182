function [ILM,cup_culum,Cup_depth ] = Cup_depth( B_scan_image )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
 %=================braye estekhraje laye aval va akhare===============%
for i=B_scan1:B_scan2;
    imagei=d3_registerr(:,:,i);
      X= ridgelet_oct( imagei );
      OCT(:,:,i)=X;
%     imagei=uint8(imagei);
% imagei=filter2(h,imagei);
%  X=imagei;
   %=================braye estekhraje laye aval va akhare===============%
level =graythresh(uint8(imagei));
H=ones(10,3);BW = im2bw(uint8(imagei),level);

BW1=filter2(H,BW);
leve2 = graythresh(uint8(BW1));
BW2 = im2bw(uint8(BW1),leve2);

% ==============braye estekhraje RPE============================;
level3= 1.4 *graythresh(uint8(X));
BW3= im2bw(uint8(X),level3);

BW4=filter2(H,BW3);

level4=graythresh(uint8(BW4));
BW5=im2bw(uint8(BW4),level4);
%===============================
for j=1:512
vectorj=(BW2(:,j));
vectorj=smooth(vectorj,10);
[r1,c1]=find(vectorj>=(.98*max(vectorj)));
rpe= mean(r1);
[r,c]=find(vectorj>=(.73*max(vectorj)));
r=sort(r);
if j==1;
    max3=r(1,1)-100;
end
for h=1:size(r,1)
if r(h,1)<=max3
    r(h,:)=500; 
end
end
r=sort(r);
r1=sort(r1);
edge(j,1)=j;
edge(j,2)=r(1,1);
a=size(r1,1);
edge(j,3)=r1(a,1);
edge(j,5)=rpe;
% find the brightest pixel in image
[r2,c2]=find(uint8(BW5(r(1,1)+75:r1(a,1),j))==1);
sort(r2);
if size(r2,1)~=0
    edge(j,6)=(mean(r2)+r(1,1)+75);
else
    edge(j,6)=NaN;
end
end
h=[1 1 1 1 -1 -1 -1 -1];
edge(:,4)=edge(:,2);

 for j=21:480 
    if abs(edge(j-11,2)-edge(j,2))>=30 && abs(edge(j+10,2)-edge(j,2))>=30
     edge(j,2)=edge(j-1,2);
    end
%     if abs(edge(j-11,6)-edge(j,6))>=20 && abs(edge(j+10,6)-edge(j,6))>=20
%      edge(j,6)=edge(j-1,6);
%     end
%     if edge(j,6)>=459;
% %        if edge(j,6)>=(3/4)*size(imagei,1)
%         edge(j,6)=NaN; 
% 
%     end
 end

 RPE(:,i)=edge(:,6);
smooth1=smooth(edge(:,1),edge(:,2),20);
smooth2=smooth(edge(:,1),edge(:,3),20);
A=num2str(i);
G1=strcat('bscan shomare',A);
ILM(:,i)=smooth1;
end
%=============find the position of vessels for improve the RPE layer=======
%%%%
%============================find rpe break point==========
TF = isnan(RPE);
TS=circshift(TF,1);
TFS=TF-TS;
TFS(1,:)=0;
TFS(1:round(.1*size(imagei,2)),:)=0;
TFS(round(.9*size(imagei,2)):size(imagei,2),:)=0;
% TFS(1:round(0.2*size(imagei,2)),:)=0;
% TFS(round(0.8*size(imagei,2))saveas:size(imagei,2),:)=0;
% 
leftRPE=zeros(1,128);
rightRPE=zeros(1,128);
for i=B_scan1:B_scan2
    A=num2str(i);
G1=strcat('B _ scan shomare',A);
TFSi=TFS(:,i);
[R,C]=find(abs(TFSi)==1);
for ii=1:size(R,1)-1;
    if TFSi(R(ii),1)-TFSi(R(ii+1),1)==-2;
    if (R(ii+1)-R(ii))<=15;
TFSi(R(ii),1)=0;
TFSi(R(ii+1),1)=0;
    end
    end
end
[R,C]=find(abs(TFSi)==1);
 R=sort(R);
if size(R,1)==1 % in if male naboode cup va naboode fasele bine leftrpebreak & rightrpebreak;
     leftrpebreak=R;
rightrpebreak=R;
 else
for ii=1:size(R,1)-1;

    if TFSi(R(ii),1)-TFSi(R(ii+1),1)==2;
       
        RowRpe(ii,1)=R(ii+1)-R(ii);
        RowRpe(ii,2)=R(ii);
        RowRpe(ii,3)=R(ii+1);
    end
   
end
TFS2(:,i)=TFSi;
RowRpe1=sortrows(RowRpe,-1);
leftrpebreak=round(RowRpe1(1,2));
rightrpebreak=round(RowRpe1(1,3));
 end

leftRPE(i)=leftrpebreak;
 rightRPE(i)=rightrpebreak;
%  LeftRPE(B_scan1)=leftRPE(B_scan1);
%   RightRPE(B_scan1)=rightRPE(B_scan1);
%   if i==B_scan1+1;
%   LeftRPE(B_scan1+1)=leftRPE(B_scan1+1);
%   RightRPE(B_scan1+1)=rightRPE(B_scan1+1);
%   end
%  if i>=B_scan1+2
%      if rightRPE(i)-rightRPE(i-2)>=30;
%          RightRPE(i)=rightRPE(i-1);
%          
%      else
%          RightRPE(i)=rightRPE(i);
%      end
% %       if leftRPE(i)-leftRPE(i-2)>=30;
%           LeftRPE(i)=leftRPE(i-1);
          
%       else
%                     LeftRPE(i)=leftRPE(i);
%       end
%   end

% i=87;

end

% ===================post processing=======================================
for i=B_scan1:B_scan2;
     LeftRPE(B_scan1)=leftRPE(B_scan1);
  RightRPE(B_scan1)=rightRPE(B_scan1);
%   if i==B_scan1+1;
  LeftRPE(B_scan1+1)=leftRPE(B_scan1+1);
  RightRPE(B_scan1+1)=rightRPE(B_scan1+1);
%   end
    if abs(leftRPE(i+2)-leftRPE(i))>=15 && abs(leftRPE(i+2)-leftRPE(i+4))>=15;
%         LeftRPE(i+2)=round(mean([leftRPE(i),leftRPE(i+1),leftRPE(i+3),leftRPE(i+4)]));
if leftRPE(i+1)==LeftRPE(i+1)
         LeftRPE(i+2)=round(mean([leftRPE(i+1),leftRPE(i+3)]));
else
    LeftRPE(i+2)=round(mean([leftRPE(i),leftRPE(i+3)]));
end
    else
       LeftRPE(i+2)= leftRPE(i+2);
    end
    if abs(rightRPE(i+2)-rightRPE(i))>=15 && abs(rightRPE(i+2)-rightRPE(i+4))>=15;
%         RightRPE(i+2)=round(mean([rightRPE(i),rightRPE(i+1),rightRPE(i+3),rightRPE(i+4)]));
if RightRPE(i+1)==rightRPE(i+1)
        RightRPE(i+2)=round(mean([rightRPE(i+1),rightRPE(i+3)]));
else
     RightRPE(i+2)=round(mean([rightRPE(i),rightRPE(i+3)]));
end
    else
       RightRPE(i+2)= rightRPE(i+2);
    end
    LeftRPE(B_scan2)=leftRPE(B_scan2);
  RightRPE(B_scan2)=rightRPE(B_scan2);
   LeftRPE(B_scan2-1)=leftRPE(B_scan2-1);
  RightRPE(B_scan2-1)=rightRPE(B_scan2-1);
  
  %============================improve RPE break point with position
  %-vessels around RPE break point===================================
  
%    RPEileft=RPE(LeftRPE-30:LeftRPE+30,i)';
%  Tv = isnan(RPEileft);
%  [r1,c1]=find(Tv==1);
 vesseli=Vessel2(i-1:i+1,LeftRPE(i):LeftRPE(i)+30);
  [r,c]=find(vesseli==1);
  clear c1;
 c1=max(c)+LeftRPE(i);
 if size(c1)~=0 
 LeftRPE1(i)=c1;
  LeftRPE2(i)=c1;
 else
    LeftRPE1(i)=NaN;  
      LeftRPE2(i)= LeftRPE(i);
 end  
 vesselj=Vessel2(i-1:i+1,RightRPE(i)-30:RightRPE(i));
  [r,c]=find(vesselj==1);
  clear c1;
 c1=min(c)+RightRPE(i)-30;
 if size(c1)~=0 
 RightRPE1(i)=c1;
 RightRPE2(i)=c1;
 
 else
    RightRPE1(i)=NaN; 
     RightRPE2(i)= RightRPE(i);
 end  
     row_leftrpebreak=ILM(LeftRPE2(i),i);
  row_rightrpebreak=ILM(RightRPE2(i),i);
  
 x=edge(:,1);
 y=ILM(:,i);
p = polyfit(x,y,40);
F(:,i)=polyval(p,x);
% edge_left=F(10);
% edge_right=F(500);
if LeftRPE2(i)<=RightRPE2(i);
  center=max(F(LeftRPE2(i):RightRPE2(i),i));
else
  center=max(F(RightRPE2(i):LeftRPE2(i),i));
end
cup_row_left=row_leftrpebreak+.5*(center-row_leftrpebreak);
cup_row_right=row_rightrpebreak+.5*(center-row_rightrpebreak);
% cup_row_left=((abs((edge_left-center)/2))+edge_left);
% cup_row_right=((abs(edge_right-center)/2)+edge_right);
 [colum_center,row_center]=find(center==F(:,i));
 [cup_colum_left,y]=find(abs(ILM(LeftRPE2(i):colum_center,i)-cup_row_left)<=10);
if size(cup_colum_left)~=0;
cup_colum_left=mean(cup_colum_left)+LeftRPE2(i)-1;
else
    cup_colum_left=colum_center;
end

[cup_colum_right,y]=find(abs(ILM(colum_center:RightRPE2(i),i)-cup_row_right)<=10);
if size(cup_colum_right)~=0;
cup_colum_right=mean(cup_colum_right)+colum_center-1;
else
    cup_colum_right=colum_center;
end
cup_point(i,:)=[i,cup_colum_left,cup_colum_right,cup_row_left,cup_row_right];
disk_point(i,:)=[i,LeftRPE(i),RightRPE(i),row_leftrpebreak,row_rightrpebreak,LeftRPE2(i),RightRPE2(i)];
center_point(i,:)=[i,center];
%%%%%%%%%%%%%%%

 hold off


end