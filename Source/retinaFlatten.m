function flatImage = retinaFlatten(img_data,shifts,interpMethod)

[Y, X, Z] = ndgrid(1:size(img_data,1),1:size(img_data,2),1:size(img_data,3));
Ym = bsxfun(@plus,Y,shiftdim(shifts,-1));
    
if size(img_data,3) == 1
    flatImage = interp2(X,Y,img_data,X,Ym,interpMethod,0);
else    
    if strcmp(version,'7.13.0.564 (R2011b)')
        % Bug in this version
        flatImage = interp3(double(img_data),X,Ym,Z,interpMethod);
    else
        F = griddedInterpolant(double(img_data),interpMethod);
        flatImage = F(Ym,X,Z);
    end
    flatImage(isnan(flatImage)) = 0;
    flatImage = cast(flatImage,class(img_data));
end
