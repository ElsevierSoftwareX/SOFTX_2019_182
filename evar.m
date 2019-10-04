function noisevar = evar(y)

%EVAR   Noise variance estimation.
%   Assuming that the deterministic function Y has additive Gaussian noise,
%   EVAR(Y) returns an estimated variance of this noise.
%
%   Note:
%   ----
%   A thin-plate smoothing spline model is used to smooth Y. It is assumed
%   that the model whose generalized cross-validation score is minimal can
%   provide the variance of the additive noise. A few tests showed that
%   EVAR works very well with "not too irregular" functions.
%
%   Examples:
%   --------
%   % 1D signal
%   n = 1e6; x = linspace(0,100,n);
%   y = cos(x/10)+(x/50);
%   var0 = 0.02; % noise variance
%   yn = y + sqrt(var0)*randn(size(y));
%   evar(yn) % estimated variance
%
%   % 2D function
%   [x,y] = meshgrid(0:.01:1);
%   f = exp(x+y) + sin((x-2*y)*3);
%   var0 = 0.04; % noise variance
%   fn = f + sqrt(var0)*randn(size(f));
%   evar(fn) % estimated variance
%
%   % 3D function
%   [x,y,z] = meshgrid(-2:.05:2);
%   f = x.*exp(-x.^2-y.^2-z.^2);
%   var0 = 0.6; % noise variance
%   fn = f + sqrt(var0)*randn(size(f));
%   evar(fn) % estimated variance
%
%   % Other examples
%   Click <a href="matlab:web('http://www.biomecardio.com/matlab/evar.html')">here</a> for more examples
%
%   Note:
%   ----
%   EVAR is only adapted to evenly-gridded 1-D to N-D data.
%
%   See also VAR, STD, SMOOTHN
%
%   -- Damien Garcia -- 2008/04, revised 2010/03
%   website: <a
%   href="matlab:web('http://www.biomecardio.com')">www.BiomeCardio.com</a>

error(nargchk(1,1,nargin));

d = ndims(y);
siz = size(y);
S = zeros(siz);
for i = 1:d
    siz0 = ones(1,d);
    siz0(i) = siz(i);
    S = bsxfun(@plus,S,cos(pi*(reshape(1:siz(i),siz0)-1)/siz(i)));
end
S = 2*(d-S(:));

% N-D Discrete Cosine Transform of Y
if exist('dctn','file')
    y = dctn(y);
    y = y(:);
else
    error('MATLAB:evar:MissingFunction',...
        ['DCTN is required. Download <a href="matlab:web(''',...
        'http://www.biomecardio.com/matlab/dctn.html'')">DCTN</a>.'])
end

%
S = S.^2; y = y.^2;

% Upper and lower bounds for the smoothness parameter
N = sum(siz~=1); % tensor rank of the y-array
hMin = 1e-6; hMax = 0.99;
sMinBnd = (((1+sqrt(1+8*hMax.^(2/N)))/4./hMax.^(2/N)).^2-1)/16;
sMaxBnd = (((1+sqrt(1+8*hMin.^(2/N)))/4./hMin.^(2/N)).^2-1)/16;

% Minimization of the GCV score
fminbnd(@func,log10(sMinBnd),log10(sMaxBnd),optimset('TolX',.1));

function score = func(L)
    % Generalized cross validation score
    M = 1-1./(1+10^L*S);
    noisevar = mean(y.*M.^2);
    score = noisevar/mean(M)^2;
end

end

