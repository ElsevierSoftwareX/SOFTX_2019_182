function z = MVGaussianizeDenoise(x,mu,sigma,coef,sigma_n,Num_component,method,Wsize)
% A non-linear rule to estimate x from noisy observation y = x + noise.
% It is assumed that the distribution of x is a global Laplacian mixture model and the noise is a Gaussian.
% So the distribution of Y is LapGauss in each component of our Mix. model
% Wsize= size of window for building (Wsize x Wsize) multivariate pdf for
% each MVLapGaussPDF
% we firstly Gaussianize each component to enhance its contrast, then
% denoise it.
% method determines our rule for denoising;
            %if method=0 => without denoising
            %if method=1 => soft threshold            
            % if method=2 => wiener

N=size(x,1);
M=size(x,2);
Lcnt=Num_component;
% Compute the mixture components of noisy observation
D_final=0;
z1=0;

for layer= 1:Lcnt
        f_final(:,:,layer)= coef(layer).* max(mulvarLapGaussPDFGlob(x,mu(layer),sigma(layer),sigma_n,Wsize),eps);
        D_final=f_final(:,:,layer)+D_final;
% end
%  
% for layer= 1:Lcnt
        % Compute the Gaussianized data z
        Gaussmatch(:,:,layer)=GlobGaussianize_NLwithMu(x,mu(layer),sigma(layer),sigma_n,1);
         
        % Select denoising method
        if method==0
                S(:,:,layer)=Gaussmatch(:,:,layer);
%         elseif method==1
%                 sigmaWhole=sigma(layer).*ones(N,M);
%                 muWhole=mu(layer).*ones(N,M);
%                 S(:,:,layer) = soft2(Gaussmatch(:,:,layer),muWhole,sigma_n,sigmaWhole);
%         else
%                 muWhole=mu(layer).*ones(N,M);
%                 varWhole=(sigma(layer)^2).*ones(N,M);
%                 S(:,:,layer) = wiener_mmse2(Gaussmatch(:,:,layer),muWhole,varWhole,sigma_n);
        end
       
z1=z1+ f_final(:,:,layer)./D_final.*S(:,:,layer);
end

% Four following commands are used to prohibit the Inf and Nan values
Wsig = conv2([1:7]/7,[1:7]/7,(x).^2,'same');
Ssig = sqrt(max(Wsig-sigma_n.^2,eps));
z2=soft(x,sqrt(2)*sigma_n.^2./Ssig);
z=max((1-(isinf(z1)|isnan(z1))).*abs(z1),eps).*(max(sign(z1)+2,eps)-2)+(isinf(z1)|isnan(z1)).*z2;

