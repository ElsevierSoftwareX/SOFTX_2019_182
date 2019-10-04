function [vessel_long,vessel_thin,RESfinal5,f]= vessel_segmentation(img0,img)
img0=double(img0);
SE0=zeros(7);SE1=SE0;SE2=SE0;SE3=SE0;SE4=SE0;SE5=SE0;SE6=SE0;SE7=SE0;SE8=SE0;SE9=SE0;SE10=SE0;SE11=SE0;SE0(4,:)=1;
SE1(4,4)=1;SE1(3,7)=1;SE1(5,1)=1;SE2(4,4)=1;SE2(2,7)=1;SE2(3,6)=1;SE2(6,1)=1;SE2(5,2)=1;
x3=[1 2 3 4 5 6 7]';y3=[7 6 5 4 3 2 1]';xy3=[x3 y3];xy3_ind=sub2ind(size(SE0),xy3(:,1),xy3(:,2));SE3(xy3_ind)=1;
x4=[1 2 4 6 7]';y4=[6 5 4 3 2]';xy4=[x4 y4];xy4_ind=sub2ind(size(SE0),xy4(:,1),xy4(:,2));SE4(xy4_ind)=1;
x5=[1 4 7]';y5=[5 4 3]';xy5=[x5 y5];xy5_ind=sub2ind(size(SE0),xy5(:,1),xy5(:,2));SE5(xy5_ind)=1;
x6=[1 4 7]';y6=[3 4 5]';xy6=[x6 y6];xy6_ind=sub2ind(size(SE0),xy6(:,1),xy6(:,2));SE6(xy6_ind)=1;
x7=[1 2 4 6 7]';y7=[2 3 4 5 6]';xy7=[x7 y7];xy7_ind=sub2ind(size(SE0),xy7(:,1),xy7(:,2));SE7(xy7_ind)=1;
x8=[1 2 3 4 5 6 7]';y8=[1 2 3 4 5 6 7]';xy8=[x8 y8];xy8_ind=sub2ind(size(SE0),xy8(:,1),xy8(:,2));SE8(xy8_ind)=1;
x9=[2 3 4 5 6]';y9=[1 2 4 6 7]';xy9=[x9 y9];xy9_ind=sub2ind(size(SE0),xy9(:,1),xy9(:,2));SE9(xy9_ind)=1;
x10=[3 4 5]';y10=[1 4 7]';xy10=[x10 y10];xy10_ind=sub2ind(size(SE0),xy10(:,1),xy10(:,2));SE10(xy10_ind)=1;
SE11(:,4)=1;
se_siz7(:,:,1)=SE0;se_siz7(:,:,2)=SE1;se_siz7(:,:,3)=SE2;se_siz7(:,:,4)=SE3;se_siz7(:,:,5)=SE4;se_siz7(:,:,6)=SE5;se_siz7(:,:,7)=SE11;
se_siz7(:,:,8)=SE6;se_siz7(:,:,9)=SE7;se_siz7(:,:,10)=SE8;se_siz7(:,:,11)=SE9;se_siz7(:,:,12)=SE10;
img=(img/max(img(:)))*255;
% clear I_th0 I_th
% for i=1: size(se_siz7,3)  
%     I_th0(:,:,i)=img-min(imopen(imclose(img,se_siz7(:,:,i)),se_siz7(:,:,i)),img);
% %     I_th2(:,:,i)=non_diff_MA(I_th(:,:,1));
% %     U=[];U=I_th2(:,:,i);
% %     x=[];y=[];[x y]=find(U~=0);xy=[x y];xy_ind=[];xy_ind=sub2ind(size(U),xy(:,1),xy(:,2));
% %     th=[];th=U(xy_ind);
% %     U=[U>mean(th)*.6].*U;U=bwareaopen(U,5).*U;
%     I_th(:,:,i)=I_th0(:,:,i);%U;
% 
% end
% clear w
% for k=1:size(se_siz7,3)
% for i=1:size(img,1)
%     for j=1:size(img,2)
%         SUM=sum(I_th(i,j,:));
%         if SUM~=0
%             w(i,j,k)=I_th(i,j,k)./SUM;
%         elseif SUM==0
%             w(i,j,k)=0;%I_th(i,j,k)./10^-9;
%         end
%     end
% end
% end
% II=zeros(size(img));
% for i=1:size(img,1)
%         for j=1:size(img,2)
%             II(i,j)=sum(w(i,j,:).*(I_th(i,j,:)));
%         end
% end
II=img;
COEF0 = fdct_wrapping(img,1,1);%fdct_wrapping(img11,1,1);
%%%%%%%------------------showing different bands------------------%%%%%%%%%
E = cell(size(COEF0));
for s=1:length(COEF0)
    E{s} = cell(size(COEF0{s}));
    for w=1:length(COEF0{s})
        A = COEF0{s}{w};
        for i=1:size(A,1)
            for j=1:size(A,2)
                E{s}{w}(i,j) = 0;
            end
        end
    end
end
 f0=[];f0=II;%AS1;
f=zeros(size(img));f=f0;%bwareaopen([f0~=0],20).*f0;
scale_min=0.4;scale_max=2;
Ns=5;
for i=1:Ns
sig(i)=exp(log(scale_min)+((i-1)*(log(scale_max)-log(scale_min)))/(Ns-1));
end
clear Gxx Gxy Gyy Gx Gy beta gama G
global Gxx Gxy Gyy Gx Gy beta gama G
 beta=0.8;
 gama=1.2;
 
 I_0=[];
 for ii=1:size(COEF0{2},2)/2
    count=1;E0=E;
    for jj=2:size(COEF0,2)
       a1= mod(jj-1,2);
       if a1==0 && jj~=2
           count=2*count;end
       ll=ii*count;
       for k=1:count
           E0{jj}{ll}=COEF0{jj}{ll};
           ll=ll-1;
       end 
    end
   I_0(:,:,ii) = ifdct_wrapping(E0,1,size(img,1),size(img,2));%figure,imagesc(I_0) 

   teta(ii)=((2*pi/8-(ii-1)*pi/8)+(2*pi/8-ii*pi/8))/2;end

for kk=3:size(sig,2)
    
    sigma=sig(kk);ll=1;
 Gxx=zeros(floor(6*sigma)+1,floor(6*sigma)+1);Gxy=Gxx;Gyy=Gxx;G=Gxx;Gx=Gxx;Gy=Gxx;
    k_0=0;
    for i=-floor(3*sigma):floor(3*sigma)
        k_0=k_0+1;l_0=0;
        for j=-floor(3*sigma):floor(3*sigma)
            l_0=l_0+1;
            G(k_0,l_0)=(1/(2*pi*sigma^2))*exp(-(i^2 + j^2)/(2*sigma^2));
            Gx(k_0,l_0)=(-1*(i)/(sigma^2))*G(k_0,l_0);
            Gy(k_0,l_0)=(-1*(j)/(sigma^2))*G(k_0,l_0);
            Gxx(k_0,l_0)=(-1+((i)^2)/(sigma^2))*(1/sigma^2)*G(k_0,l_0);
            Gyy(k_0,l_0)=(-1+((j)^2)/(sigma^2))*(1/sigma^2)*G(k_0,l_0);
            Gxy(k_0,l_0)=(((i)*(j))/(sigma^4))*G(k_0,l_0);
        end
%          G(k_0,:) = max(G(k_0,:))-G(k_0,:);
%          Gx(k_0,:) = max(Gx(k_0,:))-Gx(k_0,:);Gy(k_0,:) = max(Gy(k_0,:))-Gy(k_0,:);
%          Gxx(k_0,:) = max(Gxx(k_0,:))-Gxx(k_0,:);Gyy(k_0,:) = max(Gyy(k_0,:))-Gyy(k_0,:);
%          Gxy(k_0,:) = max(Gxy(k_0,:))-Gxy(k_0,:);
    end
%     sht1 = sum((G(:))); mean = sht1/(size(G,1)*size(G,2));G = G - mean;G = G/sht1;
%     sht1 = sum(Gx(:)); mean = sht1/(size(Gx,1)*size(Gx,2));Gx = Gx - mean;Gx = Gx/sht1;
%     sht1 = sum(Gy(:)); mean = sht1/(size(Gy,1)*size(Gy,2));Gy = Gy - mean;Gy = Gy/sht1;
%     sht1 = sum(Gxx(:)); mean = sht1/(size(Gxx,1)*size(Gxx,2));Gxx = Gxx - mean;Gxx = Gxx/sht1;
%     sht1 = sum(Gyy(:)); mean = sht1/(size(Gyy,1)*size(Gyy,2));Gyy = Gyy - mean;Gyy = Gyy/sht1;
%     sht1 = sum(Gxy(:)); mean = sht1/(size(Gxy,1)*size(Gxy,2));Gxy = Gxy - mean;Gxy = Gxy/sht1;
%%%%%%%%

for ii=1:size(COEF0{2},2)/2
%     count=1;E0=E;sizecount=1;EC={};siz=[];t00=[];
%     for jj=2:size(COEF1,2)
%        a1= mod(jj-1,2);
%        if a1==0 && jj~=2
%            count=2*count;end
%        ll=ii*count;
%        for k=1:count
%            E0{jj}{ll}=COEF1{jj}{ll};
%            
%            ll=ll-1;
%        end
%      
%       
%     end
   %%%%%%%%%%I_0=zeros(size(f));I_0 = ifdct_wrapping(E0,1,256,256);%figure,imagesc(I_0) 
%    Ibot = imbothat(I_0, se(:,:,ii));
%  Itop = imtophat(I_0, se(:,:,ii));
%  I11 = imsubtract(imadd(Itop, I_0), Ibot);
   %%%%%%%%teta=((2*pi/8-(ii-1)*pi/8)+(2*pi/8-ii*pi/8))/2;
   %Q_1(:,:,ii,kk) = celldrive( I_0,sigma,teta,f,msk,se(:,:,ii));
  % Q_1(:,:,ii,kk) = celldrive2( I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii));
   %Q_1(:,:,ii,kk) = celldrive3( I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back);
  %Q_1(:,:,ii,kk)=derivative1(I_0,teta,sigma,se(:,:,ii+4),se(:,:,ii),msk);
   %Q_1(:,:,ii,kk) =medial3(I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back,mean_back1,std_back1);
   %Q_1(:,:,ii,kk) =medial1(I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back);
 % [ Q_1(:,:,ii,kk),Q_0(:,:,ii,kk),Q_00(:,:,ii,kk)] =medial2(I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back,mean_back1,std_back1,ii);
  %[ Q_1(:,:,ii,kk),Q_0(:,:,ii,kk)] =medial4(I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back,mean_back1,std_back1,ii);
  [ Q_00(:,:,ii,kk),Q_01(:,:,ii,kk),Q_02(:,:,ii,kk),img_direct1(:,:,ii),H22(:,:,ii,kk)] =medial2eslah_img(I_0(:,:,ii),sigma,teta(ii),f,ii);
  %[ Q_00(:,:,ii,kk),Q_01(:,:,ii,kk),Q_02(:,:,ii,kk),img_direct1(:,:,ii),H22(:,:,ii,kk)] =medial2eslah_color(I_0,sigma,teta,f,msk,ii);
%[inten(:,:,ii,kk)] = medial2eslah_intensitylength( I_0,sigma,teta,f,msk,se(:,:,ii+4),se(:,:,ii),mean_back,back,estimation_back,mean_back1,std_back1,ii );
end

end 

for ii=1:size(COEF0{2},2)/2
    
    
for i=1:size(f,1)
    for j=1:size(f,2)
           
            QP(i,j,ii)=max(Q_02(i,j,ii,:));
            H_final(i,j,ii)=max(H22(i,j,ii,:));
           
    end
end
QP(:,:,ii) = [QP(:,:,ii)~=0].*H_final(:,:,ii);%bwmorph([QP(:,:,ii)~=0],'bridge').*H_final(:,:,ii);

%%%%%%%%%%%%%%%%%%%%%QP(:,:,ii)=(QP_eslahfinalangio(QP(:,:,ii),f,teta(ii),H_final(:,:,ii))~=0).*H_final(:,:,ii);

%QP(:,:,ii)=QP_eslahfinal_color(QP(:,:,ii),f,teta).*H_final(:,:,ii);
%QP(:,:,ii)=bwareaopen(QP(:,:,ii)~=0,4).*img_direct1(:,:,ii);
%QP(:,:,ii)=[QP(:,:,ii)>0].*QP(:,:,ii);
%QP(:,:,ii)=bwareaopen(QP(:,:,ii),7);
%QQ(:,:,ii)=eslahimagedirect(QP(:,:,ii),ii,teta);
%QQ1(:,:,ii)=edgedetector(QP(:,:,ii),se(:,:,ii));

%QP10(:,:,ii) = fast_direct( fast1(:,:,ii),QP(:,:,ii),ii );
%QPP1(:,:,ii)=eslah(QP(:,:,ii),ii);
%QPP(:,:,ii) = hazfcomponent( QP(:,:,ii),ii );
end

AS1=f;

%[RESfinal5] = valadationnewnew( AS1,f,QP,img_direct1,se_siz7,F_norm);
[vessel_long,vessel_thin,RESfinal5] = valadation_direct_malihe( AS1,f,QP,img_direct1,se_siz7);
RESfinal5=(RESfinal5~=0);


