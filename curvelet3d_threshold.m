function C=curvelet3d_threshold(varargin)
% curvelet3d_threshold - Threshold the 3D curvelet coefficients
%
% Call          D = curvelet3d_threshold(C, P, [lvl])
%                   [] are optional inputs and default options
% 
% Inputs
%     C         Curvelet coefficients under the form :
%               C{scale_number}{band_number} is a 3D array of coefficients
%     P         a structure containing the parameters of the transform
%               obtained by P=curvelet3d_params()
%    [lvl]      the threshold to apply, if set. else P.NSigma*P.SigmaNoise
%
% Output        
%     D         thresholded coefficients 
%				

% interface function to the non-standard mex curvelet3d_threshold_mex function
% There is a way to use it without memory copy if the output is as the input : suppress the hard_copy functions

	n=length(varargin);

	if(n<2 || n>3)
		error('D=curvelet3d_threshold(C, P, [lvl]);');
	end
	if(n==2) % MEX way : threshold(C,P)
		C=hard_copy(varargin{1});
		curvelet3d_threshold_mex(C,varargin{2});
	else % not MEX way : threshold(C,P,lvl)
		P=varargin{2};
		P.NSigma=1;
		P.SigmaNoise=varargin{3};
		C=hard_copy(varargin{1});
		curvelet3d_threshold_mex(C,P);
	end
