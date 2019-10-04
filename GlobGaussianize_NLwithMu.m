
% written by Zahra Amini 1.13.2015

function image2=GlobGaussianize_NLwithMu(image1,muL,sigmaL,sigma_n,dir)

%  Gaussianization by converting Normal-Laplac formula to Gaussian
%'dir' is +1 for fwd transform (Gaussianization) and -1 for bwd.



% if nargin<5
%     dir = 1;
% end

if dir==1
%     muG=muL;
    muG=0;
    
    M = sigma_n./max(sigmaL,eps);
    image2= zeros(size(image1,1),size(image1,2));
    for i=1:size(image1,1)
        for j=1:size(image1,2)            
            y = image1(i,j)-muL;

            a=sqrt(2)/sigmaL; % alpha
            PHI=1+erf(y/(sqrt(2)*sigma_n));
            term1=exp((a/2)*(-2*y+a*sigma_n^2));
            if term1==inf
                term1=exp(709);
            elseif term1==0
                term1=exp(-745);end
            T1=erfc((y-a*sigma_n^2)/(sqrt(2)*sigma_n));
            if T1==0
                T1=erfc(27);end
            term2=(1/2)*term1*T1;
            term3=exp((a/2)*(2*y+a*sigma_n^2));
            if term3==inf
                term3=exp(709);
            elseif term3==0
                term3=exp(-745);end
            T2=erfc((y+a*sigma_n^2)/(sqrt(2)*sigma_n));
            if T2==0
                T2=erfc(27);end
            term4=(1/2)*term3*T2;
            cdfNL(i,j)=(1/2)*(PHI+(-term1+term2)+term4);
            
            sigmaG=sqrt(max(sigmaL,eps)^2+sigma_n^2);
            image2(i,j)=muG+(sigmaG*sqrt(2))* erfinv(2*cdfNL(i,j)-1);
                
        end
    end  
     
else
    
%     % Reconstrustion (inverse Gaussianization) dir=-1
%     muL=0;
    muG=muL;
    sigmaG=sqrt(max(sigmaL,eps)^2+sigma_n^2);
    image2= zeros(size(image1,1),size(image1,2));
    for i=1:size(image1,1)
        for j=1:size(image1,2)
            xx=image1(i,j);
            if xx<=muG
                image2(i,j)=muL+(sigmaL/sqrt(2)).* log(1+(erf((xx-muG)/(sigmaG*sqrt(2)))));
            else
                image2(i,j)=muL-(sigmaL/sqrt(2)).* log(1-(erf((xx-muG)/(sigmaG*sqrt(2)))));
            end
        end 
    end
   
end
end