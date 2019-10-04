function [x_d,y_d,x_m,y_m] = findleft( RGBimgl )
rgbimgl=RGBimgl;
[y_d,x_d]=find(rgbimgl==max(max(rgbimgl(250:350,460:590))));
[aaa,bbb]=find((x_d<460|x_d>590));
x_d(aaa,:)=[];
y_d(aaa,:)=[];
[aa,bb]=find((y_d<250|y_d>350));
x_d(aa,:)=[];
y_d(aa,:)=[];
x_d=mean(x_d);
y_d=mean(y_d);
hold on
plot(x_d,y_d,'b*');
[y_m,x_m]=find(rgbimgl==min(min(rgbimgl(260:350,270:400))));
[aaa,bbb]=find((x_m<270|x_m>400));
x_m(aaa,:)=[];
y_m(aaa,:)=[];
[aa,bb]=find((y_m<260|y_m>350));
x_m(aa,:)=[];
y_m(aa,:)=[];
x_m=mean(x_m);
y_m=mean(y_m);
hold on
plot(x_m,y_m,'r*');
end

