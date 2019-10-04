clear all,
% close all,
clear;
clc;
BB=num2str(682);
load('683dl.mat');
B_scan1=45;
B_scan2=size(d3,3)-45;
%  d3=d3(:,57:456,:);

[size1,size2,size3]=size(d3);
for i=B_scan1:B_scan2;
% for i=B_scan1-10:B_scan1+28;
    imagei=d3(:,:,i);
% imagei=imagei(:,57:456);
      X= ridgelet_oct( imagei,120,280,300,512);
      OCT(:,:,i)=X;
%     imagei=uint8(imagei);
% imagei=filter2(h,imagei);
%  X=imagei;
   %=================braye estekhraje laye aval va akhare===============%
level =graythresh(uint8(imagei));
H=ones(10,3);BW = im2bw(uint8(imagei),level);

BW1=filter2(H,BW);
% figure(1001), imshow(BW,[]);
% figure(1001), imshow(BW1,[]);

leve2 = graythresh(uint8(BW1));
BW2 = im2bw(uint8(BW1),leve2);

% ==============braye estekhraje RPE============================;
level3=1.5*graythresh(uint8(X));
BW3= im2bw(uint8(X),level3);
% figure, imshow(imagei,[]);
% figure, imshow(BW3,[]);
% figure, imshow(X,[]);
% figure, imshow(ridgelet_oct( imagei,120,280,1,512),[]);
BW4=filter2(H,BW3);
% figure, imshow(BW4,[]);

level4=graythresh(uint8(BW4));
BW5=im2bw(uint8(BW4),level4);
% figure, imshow(BW5,[]);

%===============================
for j=1:size2
vectorj=(BW2(:,j));
% vectorj=smooth(vectorj,10);
[r1,c1]=find(vectorj>=(.98*max(vectorj)));
rpe= mean(r1);
[r,c]=find(vectorj>=(.73*max(vectorj)));
[r4,c4]=find(vectorj>=(.6*max(vectorj)));
% if j==300
% %     figure(1003), plot(vectorj)
% end
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
edge(j,7)=r4(end,1);
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

 for j=21:size2/2 
    if abs(edge(j-11,2)-edge(j,2))>=15 && abs(edge(j+10,2)-edge(j,2))>=15
     edge(j,2)=NaN;
    end
    if abs(edge(j-11,6)-edge(j,6))>=30 
      edge(j,6)=NaN;
    end
 end
 for j=(size2)/2:size2-32  
    if abs(edge(j+11,2)-edge(j,2))>=15 && abs(edge(j+10,2)-edge(j,2))>=15
     edge(j,2)=NaN;
    end
    if  abs(edge(j+20,6)-edge(j,6))>=30
      edge(j,6)=NaN;
    end
 end
 sk=isnan(edge(:,6));
 [ak,bk]=find(sk==0);
 smean=edge(ak,6);
k1=mean((edge(ak,6)));
for j=21:size2 
    if edge(j,6)>=k1+40
        edge(j,6)=NaN;
    end
end

 RPE(:,i)=edge(:,6);
smooth1=smooth(edge(:,1),edge(:,2),20);
smooth2=smooth(edge(:,1),edge(:,3),20);
smooth3=smooth(edge(:,1),edge(:,7),20);
A=num2str(i);
G1=strcat('bscan shomare',A);
ILM(:,i)=smooth1;
ENDlayer(:,i)=smooth3;


end

RPE1=RPE;


maxILM=max(max(ILM));
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
% leftRPE=zeros(1,size3);
% rightRPE=zeros(1,size3);
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
clear R;
clear C;
clear RowRpe;
[R,C]=find(abs(TFSi)==1);
 R=sort(R);
if size(R,1)==0||size(R,1)==1  % in if male naboode cup va naboode fasele bine leftrpebreak & rightrpebreak;
     leftrpebreak=NaN;
rightrpebreak=NaN;
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
if isnan(leftrpebreak)==1
  leftRPE(i)= leftRPE(i-1);
else
leftRPE(i)=leftrpebreak;
end
if isnan(rightrpebreak)==1
  rightRPE(i)= rightRPE(i-1);
else
 rightRPE(i)=rightrpebreak;
end

 LeftRPE(B_scan1)=leftRPE(B_scan1);
  RightRPE(B_scan1)=rightRPE(B_scan1);
  if i==B_scan1+1;
  LeftRPE(B_scan1+1)=leftRPE(B_scan1+1);
  RightRPE(B_scan1+1)=rightRPE(B_scan1+1);
  end
 if i>=B_scan1+2
     if abs(rightRPE(i)-rightRPE(i-2))>=50;
         RightRPE(i)=rightRPE(i-1);
         
     else
         RightRPE(i)=rightRPE(i);
     end
     if abs(leftRPE(i)-leftRPE(i-2))>=50;
          LeftRPE(i)=leftRPE(i-1);
          
      else
                    LeftRPE(i)=leftRPE(i);
      end
  end

% i=87;


end
% LeftRPE=leftRPE;
%  RightRPE=rightRPE;
% ===================post processing=======================================
 for i=B_scan1 :B_scan2;
     imagei=d3(:,:,i);

%%%%%%RightRPE2&LeftRPE2 bara behbode rag hast.
RightRPE2=RightRPE;
RightRPE1=RightRPE;
LeftRPE2=LeftRPE;
LeftRPE1=LeftRPE;
      row_leftrpebreak=ILM(LeftRPE2(i),i);
   row_rightrpebreak=ILM(RightRPE2(i),i);
 x=edge(:,1);
 y=ILM(:,i);
p = polyfit(x,y,40);
F(:,i)=polyval(p,x);
% % edge_left=F(10);
% % edge_right=F(500);
if LeftRPE2(i)<=RightRPE2(i);
  center=max(F(LeftRPE2(i):RightRPE2(i),i));
else
  center=max(F(RightRPE2(i):LeftRPE2(i),i));
end
cup_row_left=row_leftrpebreak+.5*(center-row_leftrpebreak);
cup_row_right=row_rightrpebreak+.5*(center-row_rightrpebreak);
% % cup_row_left=((abs((edge_left-center)/2))+edge_left);
% % cup_row_right=((abs(edge_right-center)/2)+edge_right);
 [colum_center,row_center]=find(center==F(:,i));
 centerpoints(1,i)=colum_center;
 centerpoints(2,i)=row_center;
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
ENDlayer2(1:size2,i)=RPE1(:,i);
for j=LeftRPE2(i):RightRPE2(i)
 if ENDlayer(j,i)>=max(ILM(:,i))
     ENDlayer2(j,i)=ENDlayer(j,i);
 else 
     distance_left=abs(ILM(LeftRPE2(i)-1,i)-RPE1(LeftRPE2(i)-1,i));
     distance_right=abs(ILM(RightRPE2(i),i)-RPE1(RightRPE2(i),i));
     if j<=colum_center
     ENDlayer2(j,i)=ILM(j,i)+ distance_left;
     else
         ENDlayer2(j,i)=ILM(j,i)+ distance_right; 
     end
 end
end

i;
ENDlayer3(:,i)=smooth(ENDlayer2(:,i),50);
% find DISk volume and Cup Volume
% for 
x=[LeftRPE2(i),RightRPE2(i)];
y=[ILM(LeftRPE2(i),i),ILM(RightRPE2(i),i)];
p = polyfit(x,y,1);
upper_line(:,i)=polyval(p,[1:512]);
% find left -section right-section
end






%============================== hazfe bscanha ba cup kochik================
      DISK=disk_point;
    CUP=cup_point;
    clear RadialdiskL,clear RadialdiskR,clear Radialdisk
for  i=B_scan1:B_scan2;
    imagei=d3(:,:,i);
     
%         A=num2str(i);
%         AC=strcat('F:\mydata\oct-mokhtary\mokhtari\data_oct\result\',BB);
%         mkdir(AC);
%         AD=strcat(AC,'\b-scan',A);
% G1=strcat(A);
% G2=strcat(' rigelete B _ scan shomare',A);

f1=figure(2),
 imshow(imagei,[]),title(G1);
 hold on 
 plot(edge(:,1),ILM(:,i),'r');
%  plot(edge(:,1),RPE(:,i),'r');
%  plot(edge(:,1),ENDlayer2(:,i),'r');
  plot(edge(:,1),ENDlayer3(:,i),'g');
  plot(leftRPE(i),1:650);
 plot(rightRPE(i),1:650);
plot(LeftRPE(i),1:650);
 plot(RightRPE(i),1:650);

plot(LeftRPE(i),ILM(LeftRPE(i),i),'*b');
plot(RightRPE(i),ILM(RightRPE(i),i),'*b');
% plot(LeftRPE1(i),ILM(LeftRPE(i),i),'*g');
% plot(RightRPE1(i),ILM(RightRPE(i),i):ENDlayer3(RightRPE(i),i),'.g');
% plot(RightRPE1(i),ILM(RightRPE(i),i),'*g');
% plot(LeftRPE1(i),ILM(LeftRPE(i),i):ENDlayer3(LeftRPE(i),i),'.g');
% plot(size2/2,ILM(size2/2,i),'*g');
% plot(LeftRPE1(i):RightRPE1(i),upper_line(LeftRPE1(i):RightRPE1(i),i),'.g');
% plot(LeftRPE2(i),ENDlayer3(LeftRPE2(i),i),'.r');

hold off

 pause(.5);

    end




