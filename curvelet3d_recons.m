% curvelet3d_recons - Backward 3D curvelet transform
%
% Call          Y = curvelet3d_recons(C, P)
% 
% Inputs
%     C         Curvelet coefficients under the form 
%               C{scale_number}{band_number} is a 3D array of coefficients
%     P         a structure containing the parameters of the transform
%               obtained by P=curvelet3d_params()
%
% Output
%     Y        a simple precision 3D datacube
%
