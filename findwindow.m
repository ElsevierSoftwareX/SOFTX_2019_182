 function [A]=findwindow(x_position,y_position, THETA, radial, rangY)
% clear Y_2;
% clear all
% x_position=30
% y_position=50;
% THETA=-3.75*(45-1)+90;
% rangX=10; 
% rangY=10;
for i=-1*radial:radial
    for j=-rangY:rangY
    X_1(i+radial+1,j+rangY+1)=x_position+(i*cos((THETA/180)*pi));
%   Y(i)=y_position-(r*sin((THETA/180)*pi));
    Y_1(i+radial+1,j+rangY+1)=y_position-(i*sin((THETA/180)*pi))+j;
    end
end
Y_1=Y_1';
A1=Y_1(:);
X_1=X_1';
A2=X_1(:);
A=(round([A1,A2]));
end