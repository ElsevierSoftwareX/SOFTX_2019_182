function mask=msk(RGB);
RGB=imresize(RGB,[256,256]);
RGB=imresize(RGB,[256,256]);
RED=RGB(:,:,1);
k=ones(256,256);
for i=1:256
    for j=1:256
        if RED(i,j)<35
            k(i,j)=0;
        end
    end
end
se = strel('square',3);
I_opened = imopen(k,se);
IM2 = imclose(I_opened,se);
mask = imerode(IM2,se);
figure
imshow(RGB)
figure
imshow(mask)
% im=RGB;
% im1=rgb2gray(RGB);
% im1=medfilt2(im1,[3 3]); %Median filtering the image to remove noise%
% BW = edge(im1,'sobel'); %finding edges 
% [imx,imy]=size(BW);
% msk=[0 0 0 0 0 0 0 0 0 0 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 1 1 1 1 1 1 1 1 1 0;
%      0 0 0 0 0 0 0 0 0 0 0;];
% B=conv2(double(BW),double(msk));
% B=imresize(B,[256,256]);
% se = strel('line',11,90);
% B = imdilate(B,se);
% % figure,imshow(B);
% mask=ones(256,256);
% for i=1:256
%     j=0;
%     d=B(i,j+1);
%     while  (d==0) && (j<150)
%            j=1+j;
%            mask(i,j)=0;
%            d=B(i,j);
%     end
%     
% end
% for i=1:256
%     j=257;
%     d=B(i,j-1);
%     while  (d==0) && (j>200)
%            j=j-1;
%            mask(i,j)=0;
%            d=B(i,j);
%     end
%     
% end
% se = strel('disk',20);
% mask = imerode(mask,se);
% mask(1:4,1:256)=0;
% mask(252:256,1:256)=0;
% figure




