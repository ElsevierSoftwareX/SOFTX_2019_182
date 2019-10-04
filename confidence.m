function [ M1,M2 ] = Untitled2( mean,std,n,z )
M1=mean-z*std/sqrt(n);
M2=mean+z*std/sqrt(n);
end

