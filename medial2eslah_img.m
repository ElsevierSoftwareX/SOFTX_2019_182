function [SEEDres1,SEEDres22,SEEDres222,img_direct0,H22] = medial2eslah_img( x0,sigma,teta1,f,iicoef )
global Gxx Gxy Gyy Gx Gy G beta gama
% f1=zeros(256);f1=double(f1);
% max_f=max(max(f));
% min_f=min(min(f));
% for i=1:size(f,1)
%     for j=1:size(f,2)
%    f1(i,j)=256/(max_f-min_f)*(f(i,j)-min_f);
%     end
% end
% f11=f;
% for i=2:256-1
% for j=2:256-1
% if f(i,j)<0
%     esf=f(i-1:i+1,j-1:j+1);
% f11(i,j)=mean(esf(:));end
% end
% end
% back1=zeros(256);
% max_back=max(max(back));
% min_back=min(min(back));
% for i=1:size(back,1)
%     for j=1:size(back,2)
%    back1(i,j)=255/(max_back-min_back)*(back(i,j)-min_back);
%     end
% end
% max_estimation_back=max(max(estimation_back));
% min_estimation_back=min(min(estimation_back));
% for i=1:size(back,1)
%     for j=1:size(back,2)
%    estimation_back2(i,j)=255/(max_estimation_back-min_estimation_back)*(estimation_back(i,j)-min_estimation_back);
%     end
% end
% back3 = back/max(max(abs(back)));f3 = f/max(max(abs(f)));estimation_back1 = estimation_back/max(max(abs(estimation_back)));
if teta1*180/pi<0
    teta=teta1+pi;
else
 teta=teta1;end 
if teta*180/pi<22.5 
    tetahat=0+90;
elseif teta*180/pi>=22.5 & teta*180/pi<67.5
    tetahat=45+90;
elseif teta*180/pi>=67.5 & teta*180/pi<112.5
    tetahat=0;
elseif teta*180/pi>=112.5 & teta*180/pi<157.5
    tetahat=45;
elseif teta*180/pi>=157.5 & teta*180/pi<180
    tetahat=0+90;
end
%%%edge avaliye
% TMAX=0.1;
% [Img_filtered, nIter, dTT] = twodncdf(x0, TMAX,sigma);
%x00=Img_filtered;
im=zeros(size(x0));im=1*x0;%imcomplement(x0);
img_direct=(im)/(max(max(abs(x0))));
img_direct00=img_direct;%conv2(img_direct,G,'same');
 %%%%
 %x0=x0.^5;
%  Ibot = imbothat(x0, se00);
%  Itop = imtophat(x0, se00);
%  x = imsubtract(imadd(Itop, x0), Ibot);



 %%%%
 %teta=0;%teta1;
 %f11=zeros(size(img_direct));
 %f11=conv2(f1,G,'same');DEL=del2(f11);
 Ixx=conv2(img_direct,Gxx,'same');fxx=conv2(f,Gxx,'same');
Ixy=conv2(img_direct,Gxy,'same');fxy=conv2(f,Gxy,'same');
Iyy=conv2(img_direct,Gyy,'same');fyy=conv2(f,Gyy,'same');
Ix=conv2(img_direct,Gx,'same');Iy=conv2(img_direct,Gy,'same');fx=conv2(f,Gx,'same');fy=conv2(f,Gy,'same');
% r=floor(4*sigma);
% Ixx=Ixx(r:size(x,1)+r-1,r:size(x,2)+r-1);fxx=fxx(r:size(f,1)+r-1,r:size(f,2)+r-1);
% Ixy=Ixy(r:size(x,1)+r-1,r:size(x,2)+r-1);fxy=fxy(r:size(f,1)+r-1,r:size(f,2)+r-1);
% Iyy=Iyy(r:size(x,1)+r-1,r:size(x,2)+r-1);fyy=fyy(r:size(f,1)+r-1,r:size(f,2)+r-1);
% Ix=Ix(r:size(x,1)+r-1,r:size(x,2)+r-1);Iy=Iy(r:size(x,1)+r-1,r:size(x,2)+r-1);fx=fx(r:size(f,1)+r-1,r:size(f,2)+r-1);fy=fy(r:size(f,1)+r-1,r:size(f,2)+r-1);

Ixx_prim=Ixx*(cos(teta))^2 + Ixy*(sin(2*teta)) + Iyy*(sin(teta))^2;
Iyy_prim=Ixx*(sin(teta))^2 - Ixy*(sin(2*teta)) + Iyy*(cos(teta))^2;
Ixy_prim=-.5*Ixx*(sin(2*teta)) + Ixy*(cos(2*teta)) + .5*Iyy*(sin(2*teta));
Ix_prim=Ix*(cos(teta)) + Iy*(sin(teta));Iy_prim=Iy*(cos(teta)) + Ix*(-sin(teta));
fxx_prim=fxx*(cos(teta))^2 + fxy*(sin(2*teta)) + fyy*(sin(teta))^2;
fyy_prim=fxx*(sin(teta))^2 - fxy*(sin(2*teta)) + fyy*(cos(teta))^2;
fxy_prim=-.5*fxx*(sin(2*teta)) + fxy*(cos(2*teta)) + .5*fyy*(sin(2*teta));
fx_prim=fx*(cos(teta)) + fy*(sin(teta));fy_prim=fy*(cos(teta)) + fx*(-sin(teta));
grad_imgdirect=sqrt(Ix_prim.^2+Iy_prim.^2);
for m=1:size(img_direct,1)
    for n=1:size(img_direct,2)
        [eigvector,eigvalue] = eig((sigma^2)*[Ixx_prim(m,n) Ixy_prim(m,n);Ixy_prim(m,n) Iyy_prim(m,n)]);
        %[eigvectorf,eigvaluef] = eig([fxx_prim(m,n) fxy_prim(m,n);fxy_prim(m,n) fyy_prim(m,n)]);
        landa1=eigvalue(2,2);eigvector1=eigvector(:,2);%landaf1=eigvaluef(2,2);eigvectorf1=eigvectorf(:,2);
        landa2=eigvalue(1,1);eigvector2=eigvector(:,1);%landaf2=eigvaluef(1,1);eigvectorf2=eigvectorf(:,1);
        w=max(abs(landa1),abs(landa2));
        if abs(landa1)==w
            h22(m,n)=landa1;%/round(2*sigma);%/sqrt(landa1^2+landa2^2);
            %teta22(m,n)=atan((eigvector1(1,1)/eigvector1(2,1)));%tetaf22(m,n)=atan((eigvectorf1(2,1)/eigvectorf1(1,1)));
            h11(m,n)=(landa2);%/round(2*sigma);%/sqrt(landa1^2+landa2^2);
%             teta11(m,n)=atan((eigvector2(1,1)/eigvector2(2,1)));%hf22(m,n)=landaf1;hf11(m,n)=landaf2;
%             v0x(m,n)=eigvector(1,2);v0y(m,n)=eigvector(2,2);
%             vx(m,n)=v0y(m,n)*cosd(teta*180/pi)-v0x(m,n)*sind(teta*180/pi);
%             vy(m,n)=v0y(m,n)*sind(teta*180/pi)+v0x(m,n)*cosd(teta*180/pi);vy(m,n)=-vy(m,n);
             %vv=vx(m,n);vx(m,n)=vy(m,n);vy(m,n)=vv;%bayad dar zah avaz shavad na inja
        elseif abs(landa2)==w
            h22(m,n)=landa2;%/round(2*sigma);%/sqrt(landa1^2+landa2^2);
            %teta22(m,n)=atan((eigvector2(1,1)/eigvector2(2,1)));%tetaf22(m,n)=atan((eigvectorf2(2,1)/eigvectorf2(1,1)));
            h11(m,n)=(landa1);%/round(2*sigma);%/sqrt(landa1^2+landa2^2);
%             teta11(m,n)=atan((eigvector1(1,1)/eigvector1(2,1)));%hf22(m,n)=landaf2;hf11(m,n)=landaf1;
%              v0x(m,n)=eigvector(1,1);v0y(m,n)=eigvector(2,1);
%               vx(m,n)=v0y(m,n)*cosd(teta*180/pi)-v0x(m,n)*sind(teta*180/pi);
%             vy(m,n)=v0y(m,n)*sind(teta*180/pi)+v0x(m,n)*cosd(teta*180/pi);vy(m,n)=-vy(m,n);
%            %vv=vx(m,n);vx(m,n)=vy(m,n);vy(m,n)=vv;%bayad dar zah avaz shavad na inja
        end

    end
end
h222=[h22<=0].*h22;h222=abs(h222);h222=h222/(max(h222(:)));
SS=sqrt(h22.^2+h11.^2);SS=(SS/max(SS(:))).*img_direct;
%%%%%%%mal f
% for m=1:size(img_direct,1)
%     for n=1:size(img_direct,2)
%         %[eigvector,eigvalue] = eig([Ixx_prim(m,n) Ixy_prim(m,n);Ixy_prim(m,n) Iyy_prim(m,n)]);
%         [eigvectorf,eigvaluef] = eig([fxx_prim(m,n) fxy_prim(m,n);fxy_prim(m,n) fyy_prim(m,n)]);
%         landaf1=eigvaluef(2,2)/sqrt(eigvaluef(2,2)^2+eigvaluef(1,1)^2);eigvectorf1=eigvectorf(:,2);
%         landaf2=eigvaluef(1,1)/sqrt(eigvaluef(2,2)^2+eigvaluef(1,1)^2);eigvectorf2=eigvectorf(:,1);
%         w=max(abs(landaf1),abs(landaf2));
%         if abs(landaf1)==w
%             hf22(m,n)=landaf1;tetaf22(m,n)=atan((eigvectorf1(2,1)/eigvectorf1(1,1)));%tetaf22(m,n)=atan((eigvectorf1(2,1)/eigvectorf1(1,1)));
%             hf11(m,n)=(landaf2);tetaf11(m,n)=atan((eigvectorf2(2,1)/eigvectorf2(1,1)));%hf22(m,n)=landaf1;hf11(m,n)=landaf2;
%             v0xf(m,n)=eigvectorf(1,1);v0yf(m,n)=eigvectorf(2,1);
%             vxf(m,n)=v0xf(m,n)*cosd(teta*180/pi)-v0yf(m,n)*sind(teta*180/pi);
%             vyf(m,n)=v0xf(m,n)*sind(teta*180/pi)+v0yf(m,n)*cosd(teta*180/pi);vyf(m,n)=-vyf(m,n);
%             
%         elseif abs(landaf2)==w
%             hf22(m,n)=landaf2;tetaf22(m,n)=atan((eigvectorf2(2,1)/eigvectorf2(1,1)));%tetaf22(m,n)=atan((eigvectorf2(2,1)/eigvectorf2(1,1)));
%             hf11(m,n)=(landaf1);tetaf11(m,n)=atan((eigvectorf1(2,1)/eigvectorf1(1,1)));%hf22(m,n)=landaf2;hf11(m,n)=landaf1;
%             v0xf(m,n)=eigvectorf(1,2);v0yf(m,n)=eigvectorf(2,2);
%               vxf(m,n)=v0xf(m,n)*cosd(teta*180/pi)-v0yf(m,n)*sind(teta*180/pi);
%             vyf(m,n)=v0xf(m,n)*sind(teta*180/pi)+v0yf(m,n)*cosd(teta*180/pi);vyf(m,n)=-vyf(m,n);
%         end
%     end
% end
%%%%end
%  h222=zeros(256);R=h222; vesselnes=R;
% vesselnes=zeros(size(img_direct));  
% for i=1:256
%     for j=1:256
%         if h22(i,j)<0 && msk(i,j)~=0 & img_direct(i,j)>0
%             R(i,j)=abs(h11(i,j)/h22(i,j));
%             h222(i,j)=abs(h22(i,j));l=sqrt(h11(i,j)^2+h22(i,j)^2);
%             vesselnes(i,j)=exp(-R(i,j)^2/(2*(beta^2))).*(1-exp(-l^2/(2*(gama^2))));
%         end
%     end
% end
% %%%%%%grad
% 
% grad=sqrt(fx_prim.^2+fy_prim.^2);cc=sqrt(h22.^2+h11.^2);
% gradx=sqrt(Ix_prim.^2+Iy_prim.^2);
% sobel_x=[1 2 1;0 0 0;-1 -2 -1];
% sobel_y=[1 0 -1;2 0 -2;1 0 -1];
% grad_test=sqrt(conv2(f11,sobel_x,'same').^2+conv2(f11,sobel_y,'same').^2);
% 
% %mur=(abs(h22)-abs(h11))./(abs(h22)+abs(h11));
% %%%%%%%%%%%%%%%med
%  med=zeros(size(img_direct));b1=med;
%  b2=b1;%h22f=b1;h11f=b1;
% for m=9:256-9
%     for n=9:256-9
%     b2(m,n)= grad_test( m+round(2*sigma)* round(vy(i,j)),n+round(2*sigma)* round(vx(i,j)));
%      b1(m,n)= grad_test( m+round(2*sigma)* round(-vy(i,j)),n+round(2*sigma)* round(-vx(i,j)));
%         %rr=(2*sigma*round(vect));
%         
%            p1=exp(-(1-(b1(m,n)/((b1(m,n)+b2(m,n))/2)))/(2*sigma^2));p2=exp(-(1-(b2(m,n)/((b1(m,n)+b2(m,n))/2)))/(2*sigma^2));
%            med(m,n)=1/2*(p1*b1(m,n)+p2*b2(m,n));
%     end
% end
%   med1=zeros(size(img_direct));med1=med-grad_test;med2=zeros(size(img_direct));
% for i=1:size(img_direct,1)
%     for j=1:size(img_direct,2)
%         if  abs(med1(i,j))>grad_test(i,j) %&& x(i,j)>0
%             med2(i,j)=img_direct(i,j);
%         end
%     end
% end
% %%%%%%%%end med

    res5=zeros(size(img_direct));
for i=2:size(img_direct,1)-1
    for j=2:size(img_direct,2)-1
        if      img_direct00(i,j)>0 && Ixx_prim(i,j)<0 && h22(i,j)<0 %& h11(i,j)<0 %& abs(h11(i,j))>b1%abs(h22(i,j))/2%& Ix_prim(i+round(-vy(i,j)),j+round(-vx(i,j)))*Ix_prim(i+round(vy(i,j)),j+round(vx(i,j)))<0 % & teta11(i,j)>teta-11.25*pi/180 & teta11(i,j)<teta+11.25*pi/180 % && med2(i,j)~=0%&& h22f(i,j)<0 %&& (teta22(i,j))<=pi/16 && (teta22(i,j))>=-pi/16
             res5(i,j)=img_direct(i,j);
       
        end
    end
end
%%%%ezafe


zah=zeros(size(res5));
clear vect;
% m_hat=ones(1,floor(2*sigma));
% m_hat(1,floor(2*sigma)/2+1:floor(2*sigma))=-1;
gradfinal=res5;%gradcalculatefinal(tetahat,f,res5,sigma);
for i=1+round(2*sigma):256-round(2*sigma)%9:256-9
    for j=1+round(2*sigma):256-round(2*sigma)%9:256-9
  if  tetahat==45  && res5(i,j)~=0 && gradfinal(i,j)~=0 && img_direct(i,j)>0  && Ixx_prim(i,j)<0 && Ix_prim(i+round(1),j+round(-1))*Ix_prim(i+round(-1),j+round(1))<0    && f(i,j)>f(i+round(-1),j+round(1))  && f(i,j)>f(i+round(1),j+round(-1))
      zah(i,j)=img_direct(i,j);

  elseif  tetahat==135  && res5(i,j)~=0 && gradfinal(i,j)~=0 && img_direct(i,j)>0  && Ixx_prim(i,j)<0 && Ix_prim(i+round(1),j+round(1))*Ix_prim(i+round(-1),j+round(-1))<0   &&  f(i,j)>f(i+round(-1),j+round(-1))  && f(i,j)>f(i+round(1),j+round(1)) 
          zah(i,j)=img_direct(i,j);

        

      
  elseif    tetahat==90 && res5(i,j)~=0 && gradfinal(i,j)~=0  && img_direct(i,j)>0 && Ixx_prim(i,j)<0 && Ix_prim(i+round(-1),j+round(0))*Ix_prim(i+round(1),j+round(0))<0   && f(i,j)>f(i+round(1),j+round(0))  && f(i,j)>f(i+round(-1),j+round(0)) %&& Ixx_prim(i,j)<0 %& Ixx_prim(i+round(vy(i,j)),j+round(vx(i,j)))>0 & Ixx_prim(i+round(-vy(i,j)),j+round(-vx(i,j)))>0   %&& vesselnes(i,j)>vesselnes(i-round(vy(i,j)),j-round(vx(i,j))) & vesselnes(i,j)>vesselnes(i+round(vy(i,j)),j+round(vx(i,j)))%& abs(Ixx_prim(i,j))>abs(Ixx_prim(i-round(vy(i,j)),j-round(vx(i,j)))) & abs(Ixx_prim(i,j))>abs(Ixx_prim(i+round(vy(i,j)),j+round(vx(i,j))))  %& grad(i,j)<grad(xy_ind)  & grad(i,j)<grad(xy_ind) %&& abs(Ix_prim(i,j))<abs(Ix_prim(i+round(vx(i,j)),j+round(vy(i,j))))  && abs(Ix_prim(i,j))<abs(Ix_prim(i-round(vx(i,j)),j-round(vx(i,j))))
          zah(i,j)=img_direct(i,j);

        
  elseif    tetahat==0 && res5(i,j)~=0 && gradfinal(i,j)~=0 &&  img_direct(i,j)>0  && Ixx_prim(i,j)<0 && Ix_prim(i+round(0),j+round(-1))*Ix_prim(i+round(0),j+round(1))<0   && f(i,j)>f(i+round(0),j+round(1))  && f(i,j)>f(i+round(0),j+round(-1)) %&& Ixx_prim(i,j)<0 %& Ixx_prim(i+round(vy(i,j)),j+round(vx(i,j)))>0 & Ixx_prim(i+round(-vy(i,j)),j+round(-vx(i,j)))>0   %&& vesselnes(i,j)>vesselnes(i-round(vy(i,j)),j-round(vx(i,j))) & vesselnes(i,j)>vesselnes(i+round(vy(i,j)),j+round(vx(i,j)))%& abs(Ixx_prim(i,j))>abs(Ixx_prim(i-round(vy(i,j)),j-round(vx(i,j)))) & abs(Ixx_prim(i,j))>abs(Ixx_prim(i+round(vy(i,j)),j+round(vx(i,j))))%& grad(i,j)<grad(xy_ind)  & grad(i,j)<grad(xy_ind) %&& abs(Ix_prim(i,j))<abs(Ix_prim(i+round(vx(i,j)),j+round(vy(i,j))))  && abs(Ix_prim(i,j))<abs(Ix_prim(i-round(vx(i,j)),j-round(vx(i,j))))
          zah(i,j)=img_direct(i,j);
 end 
         
    end
end
zah=[zah~=0].*res5;
%% TAVAJO TAVAJO

% [L num]=bwlabel(zah~=0);len=[];
% data=regionprops(L,'BoundingBox') ; 
% for i=1:num
%      A=data(i) ; 
%        mat=A.BoundingBox ; 
%        max_mat=max(mat(3),mat(4)) ;len(1,i)=max_mat;
%     clear winvec
%     [xx yy]= find(L==i);xy=[xx yy];
%    winvec=zeros(size(xx,2),1);
%    xy_ind=sub2ind(size(f),xy(:,1),xy(:,2));
%    winvec=zah(xy_ind);
%    intensity_seg=geomean([mean(winvec)  max(winvec)]);
%    intensity(1,i)=intensity_seg; 
% end
% 
% h=zeros(1,max(len));
% for i=1:length(len)
%     r=len(i);
%     h(r)=h(r)+1;end
% zahtest=zah;
%    for i=1:num
%        if intensity(1,i)<mean(intensity)+std(intensity) & len(1,i)< 3 %& intensity(1,i)<mean(intensity)
%        clear x y
%        [x y]=(find(L==i));xy=[x y];
%        xy_ind=sub2ind(size(zah),xy(:,1),xy(:,2));
%        zahtest(xy_ind)=0;
%        end
%    end
%% 
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    zahtest1=zah;countx=-1;
%    for i2=1:64:256
%        countx=countx+1;county=0;
%        for j2=1:64:256
%            
%            window=zah(i2:i2+63,j2:j2+63);
%            [L num]=bwlabel(window~=0);
%            data=regionprops(L,'BoundingBox') ; 
% for i=1:num
%      A=data(i) ; 
%        mat=A.BoundingBox ; 
%        max_mat=max(mat(3),mat(4)) ;len(1,i)=max_mat;
%     clear winvec
%     [xx yy]= find(L==i);xy=[xx+countx*64 yy+county*64];
%    winvec=zeros(size(xx,2),1);
%    xy_ind=sub2ind(size(zah),xy(:,1),xy(:,2));
%    winvec=zah(xy_ind);
%    intensity_seg=geomean([mean(winvec)  max(winvec)]);
%    intensity(1,i)=intensity_seg; 
% end
% for i=1:num
%        if len(1,i)<round(mean(len)) && intensity(1,i)<mean(intensity)+std(intensity)
%        clear x y
%        [x y]=(find(L==i));xy=[x+countx*64 y+county*64];
%        xy_ind=sub2ind(size(zah),xy(:,1),xy(:,2));
%        zahtest1(xy_ind)=0;
%        end
% end
% county=county+1;
%        end
%    end
 %%%%%%%end zahtest1          
    
      
% zah_0_0=zah;
% for i=2:256-1
%     for j=2:256-1
%         if zah(i,j)~=0
%         win=grad_test(i-1:i+1,j-1:j+1);
%         win1=zah(i-1:i+1,j-1:j+1);[xx yy]=find(win1~=0);
%         if win(2,2)==min(win(:)) & length(xx)==1
%             zah_0_0(i,j)=0;
%         end
%         end
%     end
% end
%         
% for i=9:256-9
%     for j=9:256-9
%       
%         if res5(i,j)~=0  && (x(i,j))>(x(i-round(vx(i,j)),j-round(vy(i,j)))) && (x(i,j))>(x(i+round(vx(i,j)),j+round(vy(i,j)))) && x(i,j)>0 && Ix_prim(i-round(vx(i,j)),j-round(vy(i,j)))*Ix_prim(i+round(vx(i,j)),j+round(vy(i,j)))<0   && DEL(i,j)<0 && vesselnes(i,j)>vesselnes(i-round(vx(i,j)),j-round(vy(i,j))) && vesselnes(i,j)>vesselnes(i+round(vx(i,j)),j+round(vy(i,j))) && grad(i,j)<grad(i+round(sigma)*round(vx(i,j)),j+round(sigma)*round(vy(i,j)))  && grad(i,j)<grad(i-round(sigma)*round(vx(i,j)),j-round(sigma)*round(vy(i,j))) %&& abs(Ix_prim(i,j))<abs(Ix_prim(i+round(vx(i,j)),j+round(vy(i,j))))  && abs(Ix_prim(i,j))<abs(Ix_prim(i-round(vx(i,j)),j-round(vx(i,j))))
%           zah(i,j)=vesselnes(i,j);  
%         end
%     end
% end
% for i=2:256-1
%     for j=2:256-1
%         if res5(i,j)~=0   && tetahat==90  && abs(Ix_prim(i,j))<abs(Ix_prim(i-1,j))  && abs(Ix_prim(i,j))<abs(Ix_prim(i+1,j)) && x(i,j)>0  && x(i,j)>x(i-1,j) && x(i,j)>x(i+1,j) && vesselnes(i,j)>vesselnes(i-1,j) && vesselnes(i,j)>vesselnes(i+1,j) &&  Ix_prim(i-1,j)*Ix_prim(i+1,j)<0 && grad(i,j)<grad(i-1,j)  && grad(i,j)<grad(i+1,j) && DEL(i,j)<0  %&& abs(DEL(i,j))>abs(DEL(i-1,j)) && abs(DEL(i,j))>abs(DEL(i+1,j))%&& DEL(i+1,j)*DEL(i-1,j)>0&& f(i,j)>f(i-1,j) && f(i,j)>f(i+1,j) %%%
%             zah(i,j)=vesselnes(i,j);
%         elseif res5(i,j)~=0   && tetahat==0  && abs(Ix_prim(i,j))<abs(Ix_prim(i,j-1))  && abs(Ix_prim(i,j))<abs(Ix_prim(i,j+1)) && x(i,j)>0  && x(i,j)>x(i,j-1) && x(i,j)>x(i,j+1) && vesselnes(i,j)>vesselnes(i,j-1) && vesselnes(i,j)>vesselnes(i,j+1) &&  Ix_prim(i,j-1)*Ix_prim(i,j+1)<0 && f(i,j)>f(i,j-1) && grad(i,j)<grad(i,j-1)  && grad(i,j)<grad(i,j+1) && DEL(i,j)<0  %&& abs(DEL(i,j))>abs(DEL(i,j+1)) && abs(DEL(i,j))>abs(DEL(i,j-1))%&& DEL(i,j-1)*DEL(i,j+1)>0&& f(i,j)>f(i,j+1) %&& grad(i,j)<grad(i,j-1)  && grad(i,j)<grad(i,j+1)%%
%             zah(i,j)=vesselnes(i,j);
%         elseif res5(i,j)~=0   && tetahat==45  && abs(Ix_prim(i,j))<abs(Ix_prim(i-1,j+1))  && abs(Ix_prim(i,j))<abs(Ix_prim(i+1,j-1)) && x(i,j)>0 && x(i,j)>x(i-1,j+1) && x(i,j)>x(i+1,j-1) && vesselnes(i,j)>vesselnes(i-1,j+1) && vesselnes(i,j)>vesselnes(i+1,j-1) &&  Ix_prim(i+1,j-1)*Ix_prim(i-1,j+1)<0 && grad(i,j)<grad(i-1,j+1)  && grad(i,j)<grad(i+1,j-1) && DEL(i,j)<0  %&& abs(DEL(i,j))>abs(DEL(i-1,j+1)) && abs(DEL(i,j))>abs(DEL(i+1,j-1))%&& DEL(i+1,j-1)*DEL(i-1,j+1)>0 && f(i,j)>f(i-1,j+1) && f(i,j)>f(i+1,j-1)%%%
%         zah(i,j)=vesselnes(i,j);
%         elseif res5(i,j)~=0   && tetahat==135  && abs(Ix_prim(i,j))<abs(Ix_prim(i-1,j-1)) && vesselnes(i,j)>vesselnes(i-1,j-1) && vesselnes(i,j)>vesselnes(i+1,j+1) && abs(Ix_prim(i,j))<abs(Ix_prim(i+1,j+1)) && x(i,j)>0 && x(i,j)>x(i-1,j-1) && x(i,j)>x(i+1,j+1)  &&  Ix_prim(i+1,j+1)*Ix_prim(i-1,j-1)<0  && grad(i,j)<grad(i-1,j-1)  && grad(i,j)<grad(i+1,j+1) && DEL(i,j)<0  %&& abs(DEL(i,j))>abs(DEL(i-1,j-1)) && abs(DEL(i,j))>abs(DEL(i+1,j+1))%&& DEL(i-1,j-1)*DEL(i+1,j+1)>0 && f(i,j)>f(i-1,j-1) && f(i,j)>f(i+1,j+1)%
%             zah(i,j)=vesselnes(i,j);
%         end
%     end
% end

      
% zah1=zah;
% max_zah=max(max(zah));
% min_zah=min(min(zah));
% for i=1:size(zah,1)
%     for j=1:size(zah,2)
%    zah_new(i,j)=1/(max_zah-min_zah)*(zah(i,j)-min_zah);
%     end
% end
%zah.*CC.*h222;
zahtest=zah;
zah01=(zahtest~=0).*SS;%(abs(h22)./max(max(abs(h22))));%test;%[zahtest~=0].*img_direct;
%zah01=[zah1~=0].*[h222/max(h222(:))];
%zah11=zah01.*cc.*h222;
clear xx yy xy xyind thresh
[xx yy]=find(zah01~=0);xy=[xx yy];xyind=sub2ind(size(zah01),xy(:,1),xy(:,2));
thres=zah01(xyind);

% zah2=[zah01>TH ];
[counts,bins] = imhist(thres);
property=counts/sum(counts);
clear xx yy
[xx yy]=find(counts~=0 );property1=zeros(length(xx),1);bins1=property1;
for i=1:length(xx)
property1(i,1)=property(xx(i),yy(i));bins1(i,1)=bins(xx(i),yy(i));len_bins1(1,i)=xx(i);
end
zah_001=zeros(size(zah));
if length(xx)~=0 & graythresh(bins1)~=0
[xx001 yy001]=find(bins>=graythresh(bins1));
for i=1:length(xx001)
    [mm20 nn20]=find(zah01>(xx001(i)-1.5)/255 & zah01<(xx001(i)-0.5)/255);
    mn20=[mm20 nn20];mn20_ind=sub2ind(size(zah),mn20(:,1),mn20(:,2));
    zah_001(mn20_ind)=zahtest(mn20_ind);end
[xx002 yy002]=find(bins1<graythresh(bins1));proth=[];proth=property1(xx002);[xx007 yy007]=find(proth==max(proth));
proth0=proth(xx007:xx002(end));bins5=bins(xx007:xx002(end));
th002=mean(proth0)+1*std(proth0);%%%%%THRESHOLD*********************************
proth1=abs(proth0-th002); [xx003 yy003]=find(proth1==min(proth1)); 
[xx004 yy004]=find(bins==bins5(xx003(end)));[xx005 yy005]=find(bins<graythresh(bins1));
for i=xx004:xx005(end)
    [mm30 nn30]=find(zah01>(i-1.5)/255 & zah01<(i-0.5)/255);
    mn30=[mm30 nn30];mn30_ind=sub2ind(size(zah),mn30(:,1),mn30(:,2));
    zah_001(mn30_ind)=zahtest(mn30_ind);end
end

SEED30=(zah_001~=0).*iicoef;

%SEED40 = regionmedial2eslah_new(img_direct,SEED30,h22,iicoef,tetahat ,res5,Ix_prim,f,grad,vx,vy);
SEEDres = SEED30.*res5;%bwmorph([SEED30~=0],'bridge');% bwareaopen(SEED40~=0,6);
%bwmorph(SEEDres~=0,'skel').*[res5];
SEEDres2=bwareaopen((SEEDres~=0),2);
SEEDres22=(SEEDres2~=0).*(res5);
%SEEDres22=bwmorph(SEEDres2~=0,'skel').*[res5];
H22=SS;%(abs(h22)./max(max(abs(h22))));%h11.^2+h22.^2;
SEEDres222=(SEEDres2~=0).*SS;%(abs(h22)./max(max(abs(h22))));%(h11.^2+h22.^2);
clear SEEDres1
SEEDres1=(SEED30~=0).*img_direct;
img_direct0=[img_direct00>0].*img_direct00;
