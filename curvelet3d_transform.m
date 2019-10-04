% curvelet3d_transform - Forward 3D curvelet transform
%
% Call          C = curvelet3d_transform(X, P)
% 
% Inputs
%     X         a simple or double precision 3D datacube
%     P         a structure containing the parameters of the transform
%               obtained by P=curvelet3d_params(X)
%
% Outputs
%     C         Curvelet coefficients under the form 
%               C{scale_number}{band_number} is a 3D array of coefficients
%               
