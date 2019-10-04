function MVpx = mulvarLapGaussPDFGlob(x,mu_x,sigma_x,sigma_n, Wsize)
% MVpx is the multivariate pdf of y1 = x1 + n1 where x1 is Laplacian with mean mu_x and std
% sigma_x  and n1 is a Gaussian with variance sigma_n^2
% Wsize is the size of window around of each pixel. it should be an odd
% number.
% for example if we consider Wsiz=3 (i.e. a 3x3 square that y6 is located
% in its center) we have: pdf(y6)= p(y1).p(y2).P(y3).p(y4).p(y5).p(y6).p(y7).p(y8).p(y9)

%% produce an extended image with mirroring first columns and rows
S1=size(x,1);
S2=size(x,2);
Exw=Wsize-1;
Extendx= zeros(S1+Exw,S2+Exw);
Extendx(((Exw/2)+1):end-Exw/2,((Exw/2)+1):end-Exw/2)=x;
for ii=1:Exw/2
    Extendx(ii,(Exw/2)+1:end-Exw/2)=x((Exw/2)-ii+1,:);
    Extendx(end-ii+1,(Exw/2)+1:end-Exw/2)=x(end-(Exw/2)+ii,:);
end
for jj=1:Exw/2
    Extendx(:,jj)=Extendx(:,Exw-jj+1);
    Extendx(:,end-jj+1)=Extendx(:,end-Exw+jj);
end

%% calculate the pdf for each pixel
sigma_x=sigma_x.*ones(size(Extendx,1),size(Extendx,2));
%------------------ needed constraints for sigma_x and sigma_n
% if sigma_n <= 0.17
%     if sigma_x(1,1) < 0.75  % for Global sigma_x
%         sigma_x= 0.75.*ones(size(x,1),size(x,2));end
% elseif sigma_n < 0.4 && sigma_n > 0.17
%         if sigma_x(1,1) < 0.5  % for Global sigma_x
%         sigma_x= 0.5.*ones(size(x,1),size(x,2));end
% end
%--------------------------------------------------------
M = sigma_n./max(sigma_x,eps);
Extendx1=Extendx-mu_x;
Extendx2 = Extendx1./(sqrt(2)*sigma_n);
Beta=exp(-Extendx2.^2)./(2*sqrt(2)*max(sigma_x,eps));

px= Beta.* (erfcx(M-Extendx2) + erfcx(M+Extendx2));

for i=1:size(Extendx,1)
    for j=1:size(Extendx,2)
        % first limiter
        L1=(M(i,j)+26).*(sqrt(2)*sigma_n);
        if (Extendx1(i,j))> L1 % first limiter
            px(i,j)=eps;
        elseif (Extendx1(i,j))< -L1
                px(i,j)=eps;
        end
    end
end

%% produce multivariate pdf using window around each pixel
PWindow=zeros(Wsize,Wsize); % pdf on each pixel of window
for i1=(Exw/2)+1:size(Extendx,1)-(Exw/2)
    for j1=(Exw/2)+1:size(Extendx,2)-(Exw/2)
        PWindow(:,:)=px(i1-(Exw/2):i1+(Exw/2),j1-(Exw/2):j1+(Exw/2));
        MVpx(i1-(Exw/2),j1-(Exw/2))= prod(prod(PWindow));
    end
end
end

