function [ X,RT1 ] = ridgelet_oct( im,freq1,freq2,degree1,degree2 )
[m,n]=size(im);
im1 = double(im);
im1= imresize(im1,[256,256]);
RT = ridgelet( im1,4,1);
RT1 = ridgelet( im1,4,0);
 RT5=RT1;
RT2=RT5;
RT4=zeros(size(RT1));
RT2=zeros(size(RT1));
RT2(:,freq1:freq2)=RT5(:,freq1:freq2);
% fereq1=120;
% fereq2=280;
% degre1=300;
% degre2=512;
RT3=RT2;
 RT3(1:300,:)=0;
%  RT4(1:512,:)=RT2(1:512,:);
  RT4(degree1:degree2,:)=RT2(degree1:degree2,:);

 X=iridgelet(RT4,4,0);
 X=imresize(X,[m,n]);
 
end

