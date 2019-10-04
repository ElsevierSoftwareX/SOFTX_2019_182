function [result10,result9] = Filling( RESRES,fluorangio)
clc
se=strel('disk',2);fluorangio=double(fluorangio);

img(1:30,1:30)=0;s_c=strel('disk',1);
vec=[1 2 3 4 5 6 7 8];%img=(img/max(img(:)))*256;
for i=1:size(vec,2)
    s_o=strel('disk',vec(i));
g1=zeros(size(fluorangio));g1=imopen(imclose(fluorangio,s_c),s_o);
g2(:,:,i)=fluorangio-min(g1,fluorangio);end
l=1;
for i=1:2:7
    g3(:,:,l)=(g2(:,:,i)+g2(:,:,i+1))./2;
%     for ii=1:256
%         for jj=1:256
%     g3(:,:,l)=geomean([g2(ii,jj,i),g2(ii,jj,i+1)]);%(g2(:,:,i)+g2(:,:,i+1))./2;
%         end
%     end
    l=l+1;
end
% for i=1:4
%     g3(:,:,i)=(g3(:,:,i)/max(max(g3(:,:,i))))*256;
% end
f3(:,:,1)=bitget(uint8(g3(:,:,1)),6)+bitget(uint8(g3(:,:,1)),7)+bitget(uint8(g3(:,:,1)),5)+bitget(uint8(g3(:,:,1)),4)+bitget(uint8(g3(:,:,1)),3);
f3(:,:,2)=bitget(uint8(g3(:,:,2)),6)+bitget(uint8(g3(:,:,2)),7)+bitget(uint8(g3(:,:,2)),5);%+bitget(uint8(g3(:,:,2)),4);%+
f3(:,:,3)=bitget(uint8(g3(:,:,3)),6)+bitget(uint8(g3(:,:,3)),7)+bitget(uint8(g3(:,:,3)),8)+bitget(uint8(g3(:,:,3)),5);
f3(:,:,4)=bitget(uint8(g3(:,:,4)),7)+bitget(uint8(g3(:,:,4)),8)+bitget(uint8(g3(:,:,4)),6);%+bitget(uint8(g3(:,:,4)),5);
ff3(:,:,1)=[f3(:,:,1)~=0].*double(g3(:,:,1));
ff3(:,:,2)=[f3(:,:,2)~=0].*double(g3(:,:,2));
ff3(:,:,3)=[f3(:,:,3)~=0].*double(g3(:,:,3));
ff3(:,:,4)=[f3(:,:,4)~=0].*double(g3(:,:,4));
RES5=RESRES;result6=zeros(size(RES5));
%%region growing
for kk=1:4
stop=0;
V=[];V=RES5;UU=[];UU=ff3(:,:,kk);result3=zeros(size(RES5));result3=V;result4=[];result4=RES5;test=zeros(size(RES5));l=0;
seed=[];seed=V;
while stop==0 
seed=[];seed=V;result3=result4;seednew=zeros(size(seed));l=l+1;
for i=2:256-1
    for j=2:256-1
        if seed(i,j)~=0 &  numel(find(UU(i-1:i+1,j-1:j+1)~=0))>1 
            win=[];win=UU(i-1:i+1,j-1:j+1);win(2,2)=0;win=double(win);x=[];y=[];[x y]=find(win>(mean(win(:))));%+.4*(std(win(:))));
            xx=[];yy=[];xx=[x-2+i];yy=[y-2+j];xy=[];xy=[xx yy];xy_ind=[];xy_ind=sub2ind(size(seed),xy(:,1),xy(:,2));
            seednew(xy_ind)=[result3(xy_ind)==0];end
    end
end
test=seed;
result4=result4+seednew;
V=seednew;
m=[];n=[];[m n]=find(V~=0);
stop=[size(m,1)==0];
end
result6=result4+result6;  
end
%%end region growing
result6=[result6~=0];
result60=result6;

for i=2:256-1
    for j=2:256-1
        if result6(i,j)==0 & result6(i-1,j)*result6(i+1,j)*result6(i,j-1)*result6(i,j+1)~=0
            result60(i-1,j)=0;result60(i+1,j)=0;result60(i,j-1)=0;result60(i,j+1)=0;end
    end
end
result7=result60-bwareaopen(result60,25);
se=strel('square',2);
result8=result7-imopen(result7,se);
result9=result8+bwareaopen(result60,25);
result10=bwareaopen(result9,9);



