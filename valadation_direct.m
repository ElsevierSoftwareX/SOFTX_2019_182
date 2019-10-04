function [vessel_long,vessel_thin,RESfinal5] = valadation_direct( AS1,f,QP,img_direct1,se_siz7)

%% ASLI
AS2=[];AS2=f;clear QQP
max_f=max(max(f));
min_f=min(min(f));
for i=1:size(f,1)
    for j=1:size(f,2)
   f(i,j)=256/(max_f-min_f)*(f(i,j)-min_f);
    end
end
% for i=1:8
%     UU=[];UU=[QP(:,:,i)~=0]-bwareaopen(QP(:,:,i),9);UU1=bwareaopen(QP(:,:,i),9);
%     UU2=[UU~=0]-[UU.*f>max(f(:))*.6]+[UU1~=0];
%     QP(:,:,i)=UU2.*QP(:,:,i);
% end
    
ii1=0;
for ll=2:2:8
    ss1=ll+1;ii1=ii1+1;
    teta_new(ii1)=((2*pi/8-(ll-1)*pi/8)+(2*pi/8-ll*pi/8))/2;
    if teta_new(ii1)*180/pi<0
    teta_new(ii1)=teta_new(ii1)+pi;
else
 teta_new(ii1)=teta_new(ii1);end
if teta_new(ii1)*180/pi<22.5 
    tetahat_new(ii1)=0+90;
elseif teta_new(ii1)*180/pi>=22.5 && teta_new(ii1)*180/pi<67.5
    tetahat_new(ii1)=45+90;
elseif teta_new(ii1)*180/pi>=67.5 && teta_new(ii1)*180/pi<112.5
    tetahat_new(ii1)=0;
elseif teta_new(ii1)*180/pi>=112.5 && teta_new(ii1)*180/pi<157.5
    tetahat_new(ii1)=45;
elseif teta_new(ii1)*180/pi>=157.5 && teta_new(ii1)*180/pi<180
    tetahat_new(ii1)=0+90;
end
    if ss1>8
        ss1=ss1-8;end
   
    
    for i=1:size(f,1)
        for j=1:size(f,2)
    QQP(i,j,ii1)=max(QP(i,j,ll),QP(i,j,ss1));
        end
    end
   
   QQP(:,:,ii1) = bwmorph([QQP(:,:,ii1)~=0],'bridge').*QQP(:,:,ii1);
QQP(:,:,ii1)=bwmorph(QQP(:,:,ii1)~=0,'skel').*QQP(:,:,ii1);
%     QQP(:,:,ii1)=QQP(:,:,ii1)./(max(max(QQP(:,:,ii1))));
    for i=1:size(f,1)
        for j=1:size(f,2)
    img_direct2(i,j,ii1)=max([img_direct1(i,j,ll)>0].*img_direct1(i,j,ll),[img_direct1(i,j,ss1)>0].*img_direct1(i,j,ss1));
        end
     end

end


for ii=1:4
    we=zeros(size(f));we1=we;we2=we;we3=we;
    %remain_zah(:,:,ii)=zeros(size(f));
    we = QQP(:,:,ii);%bwmorph([QQP(:,:,ii)~=0],'bridge');
    intensity1=zeros(size(f));intensity3=intensity1;intensity5=intensity1;
[L,num]=bwlabel(we);
center_img(:,:,ii)=zeros(size(f));valiation_img(:,:,ii)=zeros(size(f));
count_intensity=zeros(1,num);
 %bwmorph(remain_zah(:,:,ii)~=0,'skel').*intensity1;
ggg=we;
for i=1:num
    clear winvec xx yy xy xy_ind
   [xx yy]= find(L==i);xy=[xx yy];
   winvec=zeros(size(xx,2),1);
   xy_ind=sub2ind(size(f),xy(:,1),xy(:,2));
   winvec=ggg(xy_ind);%%%%%%%%%%f or ggg???????????
   intensity_seg=geomean([mean(winvec)  max(winvec)]);
   intensity1(xy_ind)=intensity_seg;
   count_intensity(1,i)=intensity_seg;
end
clear mm threshold count_intensity1
[mm]=find(count_intensity<graythresh(count_intensity));count_intensity1=zeros(1,size(mm,2));
count_intensity1=count_intensity(mm);
threshold=geomean([std(count_intensity1)  min(count_intensity1)]);
threshold1=geomean([std(count_intensity1)  max(count_intensity1)]);
for i=1:num
    clear x y xy xyind
    [x y]=find(L==i);xy=[x y];xyind=sub2ind(size(L),xy(:,1),xy(:,2));
    if count_intensity(i)>threshold1%graythresh(count_intensity)
        intensity3(xyind)=intensity1(xyind);
    elseif count_intensity(i)<graythresh(count_intensity) & count_intensity(i)>threshold
        intensity5(xyind)=intensity1(xyind);end
end
intensity5=bwareaopen(intensity5~=0,2).*intensity5;
valiation_img(:,:,ii)=intensity5;
center_img(:,:,ii)=intensity3;
end
for ii=1:4
    clear we L num data A mat max_mat len th_len mm th_len1
    we=zeros(size(f));we=valiation_img(:,:,ii);
    [L num]=bwlabel(we);
    data=regionprops(L,'BoundingBox') ; T1=zeros(1,num);
for i=1:num
     A=data(i) ; 
       mat=A.BoundingBox ; 
       max_mat=max(mat(3),mat(4)) ;len(1,i)=max_mat;
       clear x y xyind xy
       [x y]=find(L==i);xy=[x y];xyind=sub2ind(size(f),xy(:,1),xy(:,2));
       T1(i)=we(x(1),y(1));
end
clear mm
th_len=3;%round(mean(len)-.4*std(len));
th_len1(ii)=th_len;
mm=find(len<th_len);T2=zeros(1,size(mm,2));
for jj=1:size(mm,2)
    T2(jj)=T1(mm(jj));end
newimg=zeros(size(f));
for i=1:num
    if len(1,i)>th_len & T1(i)>geomean([std(T2) min(T2)])
    clear x y xy xyind
    [x y]=find(L==i);xy=[x y];xyind=sub2ind(size(f),xy(:,1),xy(:,2));
    newimg(xyind)=we(xyind);end
end
valiation_img2(:,:,ii)=zeros(size(f));valiation_img2(:,:,ii)=newimg;
end
RES1=zeros(size(f));
for i=1:4
valiation_img000(:,:,i)=bwareaopen(valiation_img(:,:,i),4);end
for i=1:256
    for j=1:256
        RES1(i,j)=max(valiation_img000(i,j,:));
    end
end
RES2=zeros(size(f));
for i=1:256
    for j=1:256
        RES2(i,j)=max(center_img(i,j,:));
    end
end

for i=2:256-1
    for j=2:256-1
        if RES1(i,j)==0 & RES1(i-1,j)*RES1(i+1,j)*RES1(i,j-1)*RES1(i,j+1)~=0
            RES1(i-1,j)=0;RES1(i+1,j)=0;RES1(i,j-1)=0;RES1(i,j+1)=0;end
    end
end
se=strel('square',2);RES1=[RES1~=0];
clear x y xy xyind THH

[x y]=find(RES1~=0);xy=[x y];xyind=sub2ind(size(f),xy(:,1),xy(:,2));
THH=f(xyind);
%RES1=[f>(mean(THH)-1.3*std(THH))].*RES1;%%%%%%%%%%%%%%%%%%%%%%%%%&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&mitone hazf she
RES01=RES1-imopen(RES1,se);RES02=RES01;%bwareaopen(RES01,3);
RES03=RES02-bwareaopen(RES02,15);RES04=bwareaopen(RES02,15);
clear L num
[L num]=bwlabel(RES03);
RES05=RES03;
for i=1:num
    clear x y xy xyind win
    [x y]=find(L==i);xy=[x y];xyind=sub2ind(size(f),xy(:,1),xy(:,2));
    win=f(xyind);win1=RES03(xyind);l=0;
    %while stop==0
    while max(win)-min(win)>50
       l=l+1;
        win=([1-(win==max(win))].*win);
    end
    if l~=0
         RES05(xyind)=(win~=0).*win1;end
   % end
end
%RES05=bwareaopen(RES05,3);
RES06=RES05+RES04;

RES08=RES06-bwareaopen(RES06,15);RES08=bwmorph(RES08,'skel');RES001=RES08;RES09=bwareaopen(RES06,15);
for i=2:256-1
    for j=2:256-1
        clear win
        win=RES001(i-1:i+1,j-1:j+1);
        if numel(find(win~=0))>=4
           RES001(i,j)=0;end
    end
end
%RESRES01=bwareaopen(RESRES01,3);
RES010=RES09+RES001;


RES2=[RES2~=0];%RES2=bwareaopen(RES2,6);
RESRES=RES2+RES010;RESRES=(RESRES~=0);
for i=2:256-1
    for j=2:256-1
        if RESRES(i,j)==0 & RESRES(i-1,j)*RESRES(i+1,j)*RESRES(i,j-1)*RESRES(i,j+1)~=0
            RESRES(i-1,j)=0;RESRES(i+1,j)=0;RESRES(i,j-1)=0;RESRES(i,j+1)=0;end
    end
end
RESRES1=RESRES-bwareaopen(RESRES,15);RESRES1=bwmorph(RESRES1,'skel');RESRES01=RESRES1;RESRES10=bwareaopen(RESRES,15);
for i=2:256-1
    for j=2:256-1
        clear win
        win=RESRES1(i-1:i+1,j-1:j+1);
        if numel(find(win~=0))>=4
           RESRES01(i,j)=0;end
    end
end
%RESRES01=bwareaopen(RESRES01,3);
RESRES11=RESRES10+RESRES01;
RESRES12=bwareaopen(RESRES11,6);
RESfinal5=RESRES12;
vessel_long=RES2+RES04;
vessel_thin=RES001;

