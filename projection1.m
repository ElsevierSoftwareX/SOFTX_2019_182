function[projection]=projection1(d3)
%d4=d3(160:550,:,:);
[m1,m2,m3]=size(d3);
s=zeros(m3,m2);
s1=d3(:,:,76);
% s1=d4(:,:,1);
% s1=uint8(s1);
%imshow(s1(160:550,:),[]);
% s2=d3(1,:,:);
% figure
% imshow(uint8(s2))
for k=1:m3
    s1=(d3(:,:,k));
    %s1=uint8(s1);
 s2=sum(s1)./m1;
% s2=uint8(s2);
%  s2=max(s1);
%  s3=min(s1);
%  s2=s2+s3./2
s(k,:)=s2;

end
%s= imresize(s,[512,512]);
% s=uint8(s);
projection=fliplr(s);
% s=imresize(s,[1000,1000]);
figure
imshow(projection,[]);
end