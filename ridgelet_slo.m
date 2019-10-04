function [ X,RT1 ] = ridgelet_slo( im )
[m,n]=size(im);
im1 = double(im);
im1= imresize(im1,[256,256]);
RT = ridgelet( im1,4,1);
RT1 = ridgelet( im1,4,0);
level =graythresh(uint8(RT1));
    BW = im2bw(uint8(RT1),.01*level);
%     RT1=BW+RT1;
% RT2=RT5;
RT4=zeros(size(RT1));
RT2=zeros(size(RT1));
 RT2(:,120:280)=RT1(:,120:280);
RT2=RT1;
RT3=RT2;
 RT3(1:300,:)=0;
 RT4(300:512,:)=RT2(300:512,:);
%  RT4=RT2;
 X=iridgelet(RT4,4,0);
 X=imresize(X,[m,n]);
 
end

