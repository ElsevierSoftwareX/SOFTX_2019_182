% curvelet3d_params - Creates a 3D curvelet parameter structure
%
% Call          P=curvelet3d_params(X)
% 
% Inputs
%     X         A simple or double precision 3D datacube
%               Initializes P with the size of X
%
% Output        
%     P         a structure containing the parameters of the transform
%               obtained by P=curvelet3d_params()
%                
% Parameters with their default/possible values
%
%                  Nx: 0        Data x size
%                  Ny: 0        Data y size
%                  Nz: 0        Data z size
%            NbrScale: 3        Number of scales
%            NbrDir2d: 16       Number of directions in 2D equivalent
%             no_fine: 0        Set the finest scale to 0
%           no_coarse: 0        Set the coarsest scale to 0
%    threshold_coarse: 0        Threshold the coarsest scale
%          FilterType: 1        Filter the data with
%                                    1 : Hard K-Sigma thresholding
%                                    2 : Soft K-Sigma thresholding
%          SigmaNoise: -1       Noise standard deviation
%              NSigma: 3        Thresholding at NSigma*SigmaNoise
%
