function msk=mask_malihe(fluorangio)
%fluorangio=imresize(fluorangio,[256,256]);
%RGB=imresize(RGB,[256,256]);
%RED=RGB(:,:,1);
k=ones(size(fluorangio));
for i=1:size(fluorangio,1)
    for j=1:size(fluorangio,2)
        if fluorangio(i,j)<=15 
            k(i,j)=0;
           
        end
    end
end
% for i=90:170
%     for j=60:200
%         if fluorangio(i,j)<18
%             k(i,j)=0;
%         end
%     end
% end
se = strel('square',5);
I_opened = imopen(k,se);
IM2 = imclose(I_opened,se);
%se1=strel('disk',3);
msk = imerode(IM2,se);
