% [Ha Ga Gr] = ON_WAVELET_FILTERS(Hr)
% Given scaling filter Hr, this function returns the corresponding wavelet
% filter and lowpass and highpass filters on the analysis side of the
% binary orthonormal wavelet tree. Remember that the scaling filter is the
% lowpass filter on the analysis side. The filters are determined according
% to what is given in the refence below. Hence, if Hr is causal, Ha, Ga and
% Gr are noncausal. While implementing orthonormal wavelet transform, the
% delays are to be taken care accordingly.
% 
% Input Arguments:
% Hr ==> Scaling filter of a wavelet (Lowpass filter on synthesis side).
% Output Arguments:
% Ha ==> Lowpass filter on the analysis side.
% Ga ==> Highpass filter on the analysis side.
% Gr ==> Wavelet filter (Highpass filter on the synthesis side).
% 
% Hr may be obtained using MakeONFilter from wavelab toolbox.
% 
% Reference: A Wavelet Tour of Signal Processing, Stephane Mallat, 2nd Ed.
% 
% Author: Sandeep Palakkal (sandeep.dion@gmail.com)
% Orgn: IIT Madras.
% Created on: March 2, 2010.


function [ha ga gr] = on_wavelet_filters(hr)

gr = fliplr(((-1).^(0:length(hr)-1)).*hr);
ha = fliplr(hr);
ga = fliplr(gr);