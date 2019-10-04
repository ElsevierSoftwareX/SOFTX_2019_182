function output = ONHSegAMPAO(handles)

if strcmp(handles.data.format , 'vol')
    data = double(handles.data.BScan).^.25;
else
    data = double(handles.data.BScan);
end

% BB=num2str(6330);
d3_registerr=data;
projection=fliplr(projection_oct(d3_registerr));
B_scan1=floor(size(d3_registerr,3)/10);
if B_scan1<10
    B_scan1=4;
end
B_scan2=size(d3_registerr,3)-B_scan1;
if B_scan2>100
    B_scan2=100;
end

for i=B_scan1:B_scan2;
    Scan_Number = i
% for i=B_scan1+17;
    imagei=d3_registerr(:,:,i);
      X= ridgelet_oct( imagei,120,280,300,512 );
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
level3= 1.4 *graythresh(uint8(X));
BW3= im2bw(uint8(X),level3);

BW4=filter2(H,BW3);

level4=graythresh(uint8(BW4));
BW5=im2bw(uint8(BW4),level4);
%===============================
for j=1:512
vectorj=(BW2(:,j));
% vectorj=smooth(vectorj,10);
[r1,c1]=find(vectorj>=(.98*max(vectorj)));
rpe= mean(r1);
[r,c]=find(vectorj>=(.73*max(vectorj)));
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
% histogram equlaizer And vessel Extraction
% [r,c]=find(projection<=1);
% projection(r,c)=mean(mean(projection));
[ projection1 ] = histogramequlaizer( uint8(projection) );
vessels_projection= vessel_extraction( projection1 );
 vessels_projection=imresize(vessels_projection,[128,512]);
 vessels_projection1=imfill(vessels_projection,'holes');
% find the boundary of disk with Localized_seg 
I=projection;
m = false(size(I,1),size(I,2));   %-- create initial mask
m(54:74,192:380) = true;
[seg,phi] = localized_seg(I, m, 500);  %-- run segmentation
% S=edge(seg,'Canny');
% [R,C]=find(S==1);
Vessel1=vessels_projection1;
[r,c]=find(seg==1);
Vessel2=Vessel1;
Vessel2(r,c)=0;
% Vessels_oct=double(Vessels_oct);
%  Vessels_oct=imresize(Vessels_oct,[128,512]);
% imagesc(Vessels_oct);
% imshow(Vessels_oct,[]);
RPE1=RPE;
%  for i=B_scan1:B_scan2;
%     

% [r,c]=find(Vesseli==1);
% RPEi=RPE(:,i)';

% Tv = isnan(RPE(:,i))';
% [r1,c1]=find(Tv==1);
% RPEi2=RPEi;
% RPEi2(c1)=[];
% meanR=mean(RPEi2);
% for ii=1:size(c,2);
% if Tv(1,c(ii))==1
%     RPEi(1,c(ii))=meanR;
% end
%     
% end
% RPE(:,i)=RPEi';
%  end
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


%============================== hazfe bscanha ba cup kochik================
      DISK=disk_point;
    CUP=cup_point;
for i=B_scan1:B_scan2;
     if cup_point(i,3)==0;
        cup_point(i,3)=1;
    end
    if cup_point(i,2)==0;
        cup_point(i,2)=1;
    end
    if abs(disk_point(i,4)-center_point(i,2))<=30 
        cup_point(i,:)=NaN ;
        disk_point(i,:)=NaN ;
    else
        A=num2str(i);
%         AC=strcat('D:\kafieh\students\marzieh mokhtari-anpao\marzieh\result\',BB);
%         mkdir(AC);
%         AD=strcat(AC,'\b-scan',A);
G1=strcat(A);
G2=strcat(' rigelete B _ scan shomare',A);
%         f3=figure(3),title(G2);
 imagei=d3_registerr(:,:,i);
%  subplot(2,1,1)
% blue(1:40,:)=[];
% red(1:40,:)=[];
%  imshow(fliplr(projection),[50,70]);
% hold on
%  subplot('position',[0,.49,.5,.49])
% figure(99),
% imshow(OCT(:,:,i),[]),title(G2);
% hold on 
%   plot(edge(:,1),ILM(:,i),'r');
%  plot(edge(:,1),RPE(:,i),'g');
%   plot(leftRPE(i),1:650);
%  plot(rightRPE(i),1:650);
% plot(LeftRPE(i),1:650);
%  plot(RightRPE(i),1:650);
% 
% plot(LeftRPE(i),ILM(LeftRPE(i),i),'*b');
% plot(RightRPE(i),ILM(RightRPE(i),i),'*b');
% plot(LeftRPE1(i),ILM(LeftRPE(i),i),'*g');
% plot(LeftRPE1(i),1:650,'g');
% plot(RightRPE1(i),ILM(RightRPE(i),i),'*g');
% plot(RightRPE1(i),1:650,'g');
% 
% 
% %  plot(edge(:,1),200*TFS2(:,i)+400);
% % plot(x,F(:,i),'r');
% % plot(x,RPE1(:,i),'r');
% % plot(1:colum_center,cup_row_left); plot(colum_center:512,cup_row_right);
% plot(cup_point(i,3),ILM(round(cup_point(i,3)),i),'*r');
% %  plot(cup_colum_right1(i),ILM(round(cup_colum_right1(i)),i),'*r');
%  plot(cup_point(i,2),ILM(round(cup_point(i,2)),i),'*r');
% %  subplot(2,1,2);
% hold off
figure(100),
%  subplot('position',[.5,.49,.5,.49])
 
 imshow(imagei,[]),title(G1);
 hold on 
 plot(edge(:,1),ILM(:,i),'r');
 plot(edge(:,1),RPE(:,i),'g');
  plot(leftRPE(i),1:650);
 plot(rightRPE(i),1:650);
plot(LeftRPE(i),1:650);
 plot(RightRPE(i),1:650);

plot(LeftRPE(i),ILM(LeftRPE(i),i),'*b');
plot(RightRPE(i),ILM(RightRPE(i),i),'*b');
plot(LeftRPE1(i),ILM(LeftRPE(i),i),'*g');
plot(LeftRPE1(i),1:650,'g');
plot(RightRPE1(i),ILM(RightRPE(i),i),'*g');
plot(RightRPE1(i),1:650,'g');
%  plot(edge(:,1),200*TFS2(:,i)+400);
% plot(x,F(:,i),'r');
% plot(1:colum_center,cup_row_left); plot(colum_center:512,cup_row_right);
plot(cup_point(i,3),ILM(round(cup_point(i,3)),i),'*r');
%  plot(cup_colum_right1(i),ILM(round(cup_colum_right1(i)),i),'*r');
 plot(cup_point(i,2),ILM(round(cup_point(i,2)),i),'*r');
%  subplot(2,1,2);
hold off
% subplot('position',[0,0,.5,.49]) 
% figure,
% imshow((imresize(Vessel2,[650,512])),[]);
% hold on
% plot(disk_point(i,3),round(disk_point(i,1)*650/128),'*b',disk_point(i,2),round(disk_point(i,1)*650/128),'*b')
%  plot(cup_point(i,3),round(cup_point(i,1)*650/128),'*r',cup_point(i,2),round(cup_point(i,1)*650/128),'*r')
% plot(LeftRPE1(1,i),round(disk_point(i,1)*650/128),'*g',RightRPE1(1,i),round(disk_point(i,1)*650/128),'*g')
%  
%  
%  hold off
%   subplot('position',[.5,0,.5,.49]) 
% figure(101),
% imshow((imresize(projection,[650,512])),[50,70]);
% hold on
% plot(disk_point(i,3),round(disk_point(i,1)*650/128),'*b',disk_point(i,2),round(disk_point(i,1)*650/128),'*b')
%   plot(cup_point(i,3),round(cup_point(i,1)*650/128),'*r',cup_point(i,2),round(cup_point(i,1)*650/128),'*r')
% plot(LeftRPE1(1,i),round(disk_point(i,1)*650/128),'*g',RightRPE1(1,i),round(disk_point(i,1)*650/128),'*g')
%   hold off
%   saveas(f3,G1,'jpg');
%     saveas(f3,G1,'fig');

  pause(2);
    end
end
% figure(5);
% bar(B_scan1:B_scan2,LeftRPE(B_scan1:B_scan2));
% figure(6)
% bar(B_scan1:B_scan2,RightRPE(B_scan1:B_scan2));

% blue(1:40,:)=[];
% red(1:40,:)=[];
figure(102)
 imshow((projection),[50,70]);
 hold on
 plot(disk_point(:,3),disk_point(:,1),'*b',disk_point(:,2),disk_point(:,1),'*b')
  plot(cup_point(:,3),cup_point(:,1),'*r',cup_point(:,2),cup_point(:,1),'*r')
  plot(LeftRPE1,disk_point(:,1),'*g',RightRPE1,disk_point(:,1),'*g')
 
  cup_wide=sqrt((cup_point(:,2)-cup_point(:,3)).^2+(cup_point(:,4)-cup_point(:,5)).^2);
disk_wide=sqrt((disk_point(:,2)-disk_point(:,3)).^2+(disk_point(:,4)-disk_point(:,5)).^2);
CDR=cup_wide./disk_wide;
tf=isnan(CDR);
[a,b]=find(tf==1);
CDR1=CDR;
CDR(a)=[];
stdCDR=std(CDR);
meanCDR=mean(CDR);
A=num2str(stdCDR);

A1=num2str(meanCDR);
% G=strcat(BB , 'mean=',A1,'  std=',A);
% figure(4)
%  bar(B_scan1:B_scan2,cup_wide(B_scan1:B_scan2),.8,'r');
%  hold on
% bar(B_scan1:B_scan2,disk_wide(B_scan1:B_scan2),.5,'b');
% hold on
% bar(B_scan1:B_scan2,CDR1(B_scan1:B_scan2).*100,.2,'g'),legend('cup-wide','disk-wide','CDR%'),title(G);

% %======================fundus show=====================================%
% 
%  
% %======================================================================
% disk_position=zeros(128,512);
% tf= isnan(disk_point);
% [r,c]=find(double(tf)==1);
%  disk_point(r,1:7)=0;
% for ii=1:size(disk_point,1);
% if disk_point(ii,1)~=0;
% disk_position(disk_point(ii,1),round(disk_point(ii,6)))=1;
% disk_position(disk_point(ii,1),round(disk_point(ii,7)))=1;
% end
% end
% figure(17)
% disk_position= imresize(disk_position,[257,229]);
% %   disk_position=fliplr(disk_position);
% 
% imshow(disk_position,[]);
% 
% cup_position=zeros(128,512);
% tf= isnan(cup_point);
% [r,c]=find(double(tf)==1);
%  cup_point(r,1:5)=0;
% for ii=1:size(cup_point,1);
% if cup_point(ii,1)~=0;
% cup_position(cup_point(ii,1),round(cup_point(ii,2)))=1;
% cup_position(cup_point(ii,1),round(cup_point(ii,3)))=1;
% end
% end
% cup_position= imresize(cup_position,[257,229]);
% %  cup_position= fliplr(cup_position);
% figure(18)
% imshow(cup_position,[]);
% 
% %======================= find fundus alignment ============================%
%    %%%%%%%%%% disk and cup point match to fundus image;
%    cup=zeros(size(rgb_rotate));
% disk=zeros(size(rgb_rotate));
% disk(ydmean1-128:ydmean1+128,xdmean1-128:xdmean1+100)=disk_position;
% cup(ydmean1-128:ydmean1+128,xdmean1-128:xdmean1+100)=cup_position;
% [ycup,xcup]=find(cup>=.1);
% [ydisk,xdisk]=find(disk>=.1);
% figure(20),imshow(rgb_rotate);
% hold on
%  plot(xcup,ycup,'r.'); plot(xdisk,ydisk,'b.');
% % Par = CircleFitByKasa([xcup,ycup]);
% % 
% % center=[Par(1),Par(2)];
% % viscircles(center,Par(3));
% % Par1 = CircleFitByKasa([xdisk,ydisk]);
% % 
% % center1=[Par1(1),Par1(2)];
% % viscircles(center1,Par1(3));
% %%%%%%%%%%%%%%%save data%%%%%%%%%%%%%%%%%
%  AD=strcat('D:\kafieh\students\marzieh mokhtari-anpao\marzieh\result\',BB);
%  mkdir(AD)
%  AD1=strcat(AD,'\cup_point.mat');
%  save(AD1,'cup_point');
% AD2=strcat(AD,'\disk_point.mat');
% save(AD2,'disk_point');
% AD3=strcat(AD,'\ILM.mat');
% save(AD3,'ILM');
% AD4=strcat(AD,'\RPE1.mat');
% save(AD4,'RPE');
% AD4=strcat(AD,'\LeftRPE.mat');
% save(AD4,'LeftRPE');
% AD40=strcat(AD,'\LeftRPE2.mat');
% save(AD40,'LeftRPE2');
% AD4=strcat(AD,'\RightRPE.mat');
% save(AD4,'RightRPE');
% AD40=strcat(AD,'\RightRPE2.mat');
% save(AD40,'RightRPE2');
% AD5=strcat(AD,'\Rigelet_OCT.mat');
% save(AD5,'OCT');
% AD6=strcat(AD,'\RGB');
% print(AD6,'-dtiff','-r500',figure(20));
% AD7=strcat(AD,'\XYDisk.mat');
% XYdisk=[xdisk,ydisk];
% save(AD7,'XYdisk');
% AD8=strcat(AD,'\XYCup.mat');
% XYCup=[xcup,ycup];
% save(AD8,'XYCup');
% AD8=strcat(AD,'\RPE.mat');
% save(AD8,'RPE');
% AD8=strcat(AD,'\center_point.mat');
% save(AD8,'center_point');
% AD9=strcat(AD,'\projection');
% print(AD9,'-dtiff','-r500',figure(16));


% %%%%%%%%%%%%%%%save data%%%%%%%%%%%%%%%%%
output.cup_point = cup_point;
output.disk_point = disk_point;
output.ILM = ILM;
output.RPE = RPE;

uisave('output')
