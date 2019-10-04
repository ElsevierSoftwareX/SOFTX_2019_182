function[projection_oct]=projection_oct(d3)
%d4=d3(160:550,:,:);
% d3=d3_register;
[v1,v2,v3]=size(d3);
s=zeros(v3,v2);
% s1=d3(:,:,76);
% s1=d4(:,:,1);
% s1=uint8(s1);
%imshow(s1(160:550,:),[]);
% s2=d3(1,:,:);
% figure
% imshow(uint8(s2))
for k=1:v3
    s1=(d3(:,:,k));
    %s1=uint8(s1);
 s2=sum(s1)./v1;
% s2=uint8(s2);
%  s2=max(s1);
%  s3=min(s1);
%  s2=s2+s3./2
s(k,:)=s2;

end
%s= imresize(s,[512,512]);
% s=uint8(s);
projection_oct=fliplr(s);
% s=imresize(s,[1000,1000]);
% figure
% imshow(projection_oct,[]);
end