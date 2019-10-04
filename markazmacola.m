 close all
vessel_fundusr=load('D:\mydata\oct-mokhtary\mokhtari\data_oct\OD-mat\v1452rf.mat');
 load 'D:\mydata\oct-mokhtary\code oct\marzieh\1450dr.mat';
 
macola_image=rgbimgl(y_m-100:y_m+100,x_m-100:x_m+100);

figure, imshow(macola_image,[]);
vessel_fundusr=vessel_fundusr.Vessels_oct;
%  vessel_fundusr=fliplr(vessel_fundusr);
vessel=imresize(vessel_fundusr,[622,667]);
vessel_macola=vessel;
% vessel_macola=vessel_macola(y_m-100:y_m+100,x_m-150:x_m+150);
% vessel_macola=edge(vessel_macola,'canny');
v=zeros(size(vessel_macola,1)+4,size(vessel_macola,2)+4);
v(3:size(vessel_macola,1)+2,3:size(vessel_macola,2)+2)=vessel_macola;
vessel_macola1=v;

% vessel_macola=zeropad(vessel_macola);
[a,b]=find(vessel_macola1==1);
A=[a,b];
i=size(A,1);
clear s;
for e=1:i;
[c,d]= find(vessel_macola1(A(e,1)-2:A(e,1)+2,A(e,2)-2:A(e,2)+2)==1);
s(e,1)=size(c,1);
end
[s1,s2]=find(s<=9 & s>=8);

for i=1:size(s1,1);
if (((A(s1(i,1),2)-x_m)^2)+((A(s1(i,1),1)-y_m)^2))<=10000;
g=[A(s1(i,1),1),A(s1(i,1),2)];
else
    g=[0,0];
end
B(i,:)=g;
end
%
figure, imshow(vessel_macola1,[]);
hold on
plot(A(s1,2),A(s1,1),'*b');
hold off
figure, imshow(rgbimgl,[]);
hold on
plot(A(s1,2),A(s1,1),'*b');
hold off
%
figure, imshow(vessel_macola1,[]);
hold on
plot(B(:,2),B(:,1),'*b');
hold off
figure, imshow(rgbimgl,[]);
hold on
plot(B(:,2),B(:,1),'*b');
hold off