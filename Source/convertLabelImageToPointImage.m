function ptImage = convertLabelImageToPointImage(labelImg)

labelImg = int8(labelImg);

% Make sure label at the top and bottom of image are the same
if labelImg(1,2,1) ~= labelImg(end,2,1)
    labelImg(end,:,:) = 0;
end

numLabels = max(labelImg(:))+1;

% Non-zero first partial derivatives give us boundaries
label_diff = cat(1,diff(int8(labelImg),1),zeros(1,size(labelImg,2),size(labelImg,3)));
% Diff of last boundary gives 0-maxLabel
label_diff(label_diff == -(numLabels-1)) = 1;

% Find diff values >1 (when a layer disappears)
label_diff(label_diff > 1) = 1;

ptImage = uint8((labelImg+1).*label_diff);
