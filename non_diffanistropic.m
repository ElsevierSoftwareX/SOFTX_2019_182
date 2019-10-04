function QPnon = non_diffanistropic( QP1 )
scale_min=.4;scale_max=2;
Ns=5;
for i=1:Ns
sig(1,i)=exp(log(scale_min)+((i-1)*(log(scale_max)-log(scale_min)))/(Ns-1));
end
COEF=1/sqrt(2);
for kk=1:length(sig)
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
    end
 Ixx=conv2(QP1,Gxx,'same');
Ixy=conv2(QP1,Gxy,'same');
Iyy=conv2(QP1,Gyy,'same');
for m=1:size(QP1,1)
    for n=1:size(QP1,2)
        [eigvector,eigvalue] = eig([Ixx(m,n) Ixy(m,n);Ixy(m,n) Iyy(m,n)]);
        landa1=eigvalue(2,2);eigvector1=eigvector(:,2);
        landa2=eigvalue(1,1);eigvector2=eigvector(:,1);
        w=max(abs(landa1),abs(landa2));
        C(m,n,kk)=(landa1^2+landa2^2);

    end
end
%C(:,:,kk)=C(:,:,kk)/max(max(C(:,:,kk)));
end
for i=1:size(QP1,1)
    for j=1:size(QP1,2)
        C00(i,j)=max(C(i,j,:));
    end
end
C00=QP1;
h1=1;[M N]=size(QP1);Q=20;delt=.125;
K=.2;alfa=2;
for repit=1:Q
       for m=h1+1:M-h1
          for n=h1+1:N-h1
    dC1=C00(m-h1,n-h1)-C00(m,n);
    dC2=C00(m,n-h1)-C00(m,n);
    dC3=C00(m+h1,n-h1)-C00(m,n);
     dC4=C00(m-h1,n)-C00(m,n);
     dC5=C00(m+h1,n)-C00(m,n);
     dC6=C00(m-h1,n+h1)-C00(m,n);
     dC7=C00(m,n+h1)-C00(m,n);
     dC8=C00(m+h1,n+h1)-C00(m,n);
     DD(1)=exp(-abs(dC1/K))^alfa;%DD(1)=1/(1+abs(DDQP1/K)^alfa);%
     DD(2)=exp(-abs(dC2/K))^alfa;%DD(2)=1/(1+abs(DDI2/K)^alfa);%
    DD(3)=exp(-abs(dC3/K))^alfa;% DD(3)=1/(1+abs(DDI3/K)^alfa);%
     DD(4)=exp(-abs(dC4/K))^alfa;%DD(4)=1/(1+abs(DDI4/K)^alfa);%
     DD(5)=exp(-abs(dC5/K))^alfa;%DD(5)=1/(1+abs(DDI5/K)^alfa);%
    DD(6)=exp(-abs(dC6/K))^alfa;% DD(6)=1/(1+abs (DDI6/K)^alfa);%
     DD(7)=exp(-abs(dC7/K))^alfa;%DD(7)=1/(1+abs(DDI7/K)^alfa);%
     DD(8)=exp(-abs(dC8/K))^alfa;%
      C00(m,n)=C00(m,n)+delt*DD*[COEF*dC1,dC2,COEF*dC3,dC4,dC5,COEF*dC6,dC7,COEF*dC8]';
          end
       end
end
% QPnon=QP1;
QPnon=C00;%uint8(QP1);




