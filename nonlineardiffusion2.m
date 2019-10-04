% clear all
%close all
tic;
%  I=imread('a1.jpg');
%  I=rgb2gray(I);
%  I=im2double(I);
%  I1=I;

I2=p;
I1=I2;
 %st=strel('disk',5);
%I1=imtophat(I1,st);
%II=I1; 
delt=.125;
K=.151
alfa=2;
[M,N]=size(I1);
figure,imshow(I1,[])
COEF=1/sqrt(2);
h=1;
for repit=1:30
       for m=h+1:M-h
          for n=h+1:N-h
    dI1=I1(m-h,n-h)-I1(m,n);
    dI2=I1(m,n-h)-I1(m,n);
    dI3=I1(m+h,n-h)-I1(m,n);
     dI4=I1(m-h,n)-I1(m,n);
     dI5=I1(m+h,n)-I1(m,n);
     dI6=I1(m-h,n+h)-I1(m,n);
     dI7=I1(m,n+1)-I1(m,n);
     dI8=I1(m+h,n+h)-I1(m,n);
     D(1)=1/(1+abs(dI1/K)^alfa);
     D(2)=1/(1+abs(dI2/K)^alfa);
     D(3)=1/(1+abs(dI3/K)^alfa);
     D(4)=1/(1+abs(dI4/K)^alfa);
     D(5)=1/(1+abs(dI5/K)^alfa);
     D(6)=1/(1+abs (dI6/K)^alfa);
     D(7)=1/(1+abs(dI7/K)^alfa);
     D(8)=1/(1+abs(dI8/K)^alfa);
     I1(m,n)=I1(m,n)+delt*D*[COEF*dI1,dI2,COEF*dI3,dI4,dI5,COEF*dI6,dI7,COEF*dI8]';
     end
end
end

figure,imshow(I1,[]);     
imwrite (I1,'I.jpg');
toc
BW = edge(I1,'canny');
figure, imshow(BW)