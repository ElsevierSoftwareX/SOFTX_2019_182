function [ rgbimgl,vessel_fundusl,proj,ptransform1,a2,b2,thetarotate ] = register_function(RGBimgl ,vessel_fundusl,projection,vessel_projectionl)
%a2and b2=the position of projection in fundus
% thetarotate is rotate;
c=2;
d=47;
vessel_fundusl=imresize(vessel_fundusl,[620+c,620+d]);
i=1;
% % % [c,d]=find(D==max(max(D)));
% D1= zeros(11,11);
for isize=260:280
    j=1;
   for jsize=260:280
       vessel_projectionli= imresize(vessel_projectionl,[isize,jsize]);
  YY=xcorr2(double(vessel_fundusl),double(vessel_projectionli));
   D1(i,j)=max(max(YY));
   j=j+1;
   end
i=i+1;
end
% 
[c1,d1]=find(D1==max(max(D1)));
c11=c1(1,1);
d11=d1(1,1)
vessel_projectionl=imresize(vessel_projectionl,[260+c11,260+d11]);
figure(1)
imshowpair(double(vessel_fundusl),double(vessel_projectionl));
% 
p=double(vessel_projectionl);
[p1,p2]=size(vessel_projectionl);
f=double(vessel_fundusl);
[f1,f2]=size(vessel_fundusl);
%
corr_l=xcorr2(double(vessel_fundusl),double(vessel_projectionl));
figure(2)
imshow(corr_l,[])
[al,bl]=find(corr_l==max(max(corr_l)));
ptransform=zeros(f1,f2);
ptransform(al-p1+1:al,bl-p2+1:bl)=p;
figure(3)
imshowpair(f,ptransform);title('scaling and transform')
%phase3 ratation
C=zeros(1,30);
i=1;
while i<=30
for theta2=-15:15
    vessel_octi=imrotate(vessel_projectionl,theta2);
    XX=xcorr2(double(vessel_octi),f);
    C(i)=max(max(XX));
    i=i+1;
end
end
[a,b]=find(C==max(C));
thetarotate=-16+b;
ptrans_rotat=imrotate(vessel_projectionl,thetarotate);
ptrans_rotat=imresize(ptrans_rotat,[p1,p2]);
corr_ll=xcorr2(double(vessel_fundusl),double(ptrans_rotat));
figure(4)
imshow(corr_ll,[]),title('corr-rotate')
% hold on
% plot(b11,a11,'r*')
[a2,b2]=find(corr_ll==max(max(corr_ll)));
[a3,b3]=size(ptrans_rotat);
ptransform1=zeros(f1,f2);
ptransform1(a2-a3+1:a2,b2-b3+1:b2)=ptrans_rotat;
figure(5)
imshowpair(f,ptransform1),title('scaling transform rotate')

%%%%
[f1,f2]=size(f);
rgbimgl=imresize(RGBimgl,[f1,f2]);
[p1,p2]=size(p);
projection=imresize(projection,[p1,p2]);
projection=imrotate(projection,thetarotate);
projection=imresize(projection,[p1,p2]);
proj=zeros(f1,f2);
proj(al-p1+1:al,bl-p2+1:bl)=projection;
figure(6)
imshowpair(rgbimgl,proj);


end

