% vessels_projection1=imfill(vessels_projection,'holes');
function [Vessel2]=vessel_extraction_slo(slo,matrix);
vessels_slo= vessel_extraction( slo );
% % % % % % % % % % % % %  [ slo1 ] = histogramequlaizer( uint8(slo) );
% % % % % % % % % % % % % % slo1=slo;
% % % % % % % % % % % % % I=double(slo1);
% % % % % % % % % % % % % m = false(size(I,1),size(I,2));
% % % % % % % % % % % % % m(matrix(1,8)-70:matrix(1,8)+70,matrix(1,7)-70:matrix(1,7)+70) = true;
% % % % % % % % % % % % % [seg,phi] = localized_seg(I, m, 200);  %-- run segmentation
% % % % % % % % % % % % % S=edge(seg,'Canny');
% % % % % % % % % % % % % [R,C]=find(S==1);
% % % % % % % % % % % % % Vessel1=vessels_slo;
% % % % % % % % % % % % % % [r,c]=find(seg==1);
% % % % % % % % % % % % %  seg2=erfcinv(double(seg));
% % % % % % % % % % % % % Vessel2=Vessel1.*seg2;
% % % % % % % % % % % % % % Vessel2(r,c)=0;
% % % % % % % % % % % % % % figure, imshow(Vessel1,[])
% % % % % % % % % % % % % % [r,c]=find(seg==1);
% % % % % % % % % % % % % % Vessel2(r,c)=0;
Vessel2=vessels_slo;
figure, imshow(Vessel2,[]);
end