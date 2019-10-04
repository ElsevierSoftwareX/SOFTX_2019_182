function [Mu,Sigma,Coef] = MVEstLapGaussMix_Glob(x,sigma_n,Num_component,Wsize)
% Use Expectation-Maximization algorithm to find parameters of
% p(x) = K1.*MVLapGaussPDF(x,sigma1) + K2.*MVLapGaussPDF(x,sigma2)+...+Kn.*MVLapGaussPDF(x,sigman)
% where sigma_n is std of Gaussian noise.
% Wsize= size of window for building (Wsize x Wsize) multivariate pdf for
% each MVLapGaussPDF

% Minimum std for each component to prohibit numerical errors
N=size(x,1);
M=size(x,2);
MinSigma = sigma_n/20;
mux=mean(x(:));
sigmax=std(x(:));

% Initialize EM algorithm
Lcnt=Num_component;

% select initial values based on figure properties
mu(1)=.71;%*mux;
mu(2)=.28;%*mux;
mu(3)= .72;%*mux;
mu(4)=.67;%*mux;
mu(5)= .35;%*mux;
mu(6)= .24;%*mux;
mu(7)= .49;%*mux;
mu(8)= .31;%*mux;
mu(9)= .66;%*mux;
mu(10)= .58;%*mux;
mu(11)= .61;%*mux;

sigma(1)=0.1*sigmax;
sigma(2)=0.2*sigmax; 
sigma(3)=0.3*sigmax;
sigma(4)=0.4*sigmax;
sigma(5)=0.5*sigmax; 
sigma(6)=0.6*sigmax;
sigma(7)=0.7*sigmax; 
sigma(8)=0.8*sigmax; 
sigma(9)=0.9*sigmax; 
sigma(10)=0.3*sigmax; 
sigma(11)= 0.6*sigmax;

coef(1)=0.15;
coef(2)=0.28; 
coef(3)=0.3;
coef(4)=0.12;
coef(5)=0.15; 
coef(6)=0.2;
coef(7)=0.4; 
coef(8)=0.3; 
coef(9)=0.17; 
coef(10)=0.23; 
coef(11)= 0.34;
%Normalize coefs
Scoef=sum(coef(1:Lcnt));
for layer= 1:Lcnt
    coef(layer)=coef(layer)/Scoef;
end

% EM iterations
for iter = 1:5 %increase the iterations to 10 (or 20) to obtain better results
    D=0;
% E-step
    for layer= 1:Lcnt
        f(:,:,layer)= coef(layer).* max(mulvarLapGaussPDFGlob(x,mu(layer),sigma(layer),sigma_n,Wsize),eps);
        D=f(:,:,layer)+D;
    end
    
    for layer= 1:Lcnt
        r(:,:,layer)= f(:,:,layer)./D;  % responsibility factor
        % M-step
        resp(:,:)= r(:,:,layer);
        coef(layer)=max(sum(resp(:)),eps)/(size(x,1)*size(x,2));
        
        % for Gaussian Mixture Model
        q1=resp.*(2*(x-mu(layer)).^2);
        Wsig=sum(q1(:))./sum(resp(:));
        sigma(layer) =sqrt(max(Wsig,MinSigma));
        
%         q2=resp.*x;
%         mu(layer)=sum(q2(:))./sum(resp(:));
%------------------------------------------------
        % for Laplace Mixture Model:
%         q1=resp.*(sqrt(2)*(x-mu(layer)));
%         Wsig=sum(q1(:))./sum(resp(:));
%         sigma(layer) =(max(Wsig,MinSigma));
       
        q2=resp.*x./abs(x-mu(layer));
        q3=resp./abs(x-mu(layer));
        mu(layer)=sum(q2(:))./sum(q3(:));
%------------------------------------------------       
%         Pcoef(layer,iter)=coef(layer);
%         Pmu(layer,iter)=mu(layer);
%         Psigma(layer,iter)=sigma(layer);
    end
end

Mu=mu(1:Lcnt);
Sigma=sigma(1:Lcnt);
Coef=coef(1:Lcnt);