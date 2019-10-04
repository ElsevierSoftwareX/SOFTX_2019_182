function pts = convertLabelImageToPoints(labelImg)

labelImg = int8(labelImg);

% Make sure label at the top and bottom of image are the same
if labelImg(1,2,1) ~= labelImg(end,2,1)
    labelImg(end,:,:) = 0;
end

% Convert to boundary point image
ptImg = convertLabelImageToPointImage(labelImg);
numLabels = max(ptImg(:));

% Extract the points from the boundary point image
pts = zeros(size(labelImg,2),size(labelImg,3),numLabels);
prev_inds = zeros(size(labelImg,2),size(labelImg,3));
for i = 1:numLabels
    [~, inds] = max(ptImg==i,[],1);
    inds = squeeze(inds);
    % If there is no point (due to overlap), max will be the first point;
    % lets change it to the previous layer's point
    inds(inds == 1) = prev_inds(inds == 1);
    prev_inds = squeeze(inds);
    pts(:,:,i) = squeeze(inds);
end
