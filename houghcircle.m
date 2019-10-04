function [y0detect,x0detect,Accumulator] = houghcircle(Imbinary,r,thresh)
%HOUGHCIRCLE - detects circles with specific radius in a binary image.
%
%Comments:
%       Function uses Standard Hough Transform to detect circles in a binary image.
%       According to the Hough Transform for circles, each pixel in image space
%       corresponds to a circle in Hough space and vise versa. 
%       upper left corner of image is the origin of coordinate system.
%
%Usage: [y0detect,x0detect,Accumulator] = houghcircle(Imbinary,r,thresh)
%
%Arguments:
%       Imbinary - a binary image. image pixels that have value equal to 1 are
%                  interested pixels for HOUGHLINE function.
%       r        - radius of circles.
%       thresh   - a threshold value that determines the minimum number of
%                  pixels that belong to a circle in image space. threshold must be
%                  bigger than or equal to 4(default).
%
%Returns:
%       y0detect    - row coordinates of detected circles.
%       x0detect    - column coordinates of detected circles. 
%       Accumulator - the accumulator array in Hough space.
%
%Written by :
%       Amin Sarafraz
%       Photogrammetry & Computer Vision Devision
%       Geomatics Department,Faculty of Engineering
%       University of Tehran,Iran
%       sarafraz@geomatics.ut.ac.ir
%
%May 5,2004         - Original version
%November 24,2004   - Modified version,faster and better performance


if nargin == 2
    thresh = 4;
elseif thresh < 4
    error('treshold value must be bigger or equal to 4');
    return
end

%Voting
Accumulator = zeros(size(Imbinary));
[yIndex xIndex] = find(Imbinary);
for cnt = 1:length(xIndex)
    low=xIndex(cnt)-r;
    high=xIndex(cnt)+r;
    if (low<1) low=1; end
    if (high>size(Imbinary,2)) high=size(Imbinary,2); end
    for x0 = low:high
        y01 = yIndex(cnt)-sqrt(r^2-(xIndex(cnt)-x0)^2);
        y02 = yIndex(cnt)+sqrt(r^2-(xIndex(cnt)-x0)^2);
        y01 = round(y01); y02 = round(y02);
        if y01 < size(Imbinary,1) & y01 >= 1
            Accumulator(y01,x0) = Accumulator(y01,x0)+1;
        end
        if y02 < size(Imbinary,1) & y02 >= 1
            Accumulator(y02,x0) = Accumulator(y02,x0)+1;
        end
    end
end

% Finding local maxima in Accumulator
y0detect = []; x0detect = [];
AccumulatorbinaryMax = imregionalmax(Accumulator);
[Potential_y0 Potential_x0] = find(AccumulatorbinaryMax == 1);
Accumulatortemp = Accumulator - thresh;
for cnt = 1:length(Potential_y0)
    if Accumulatortemp(Potential_y0(cnt),Potential_x0(cnt)) >= 0
        y0detect = [y0detect;Potential_y0(cnt)];
        x0detect = [x0detect;Potential_x0(cnt)];
    end
end

