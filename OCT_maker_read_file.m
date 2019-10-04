% function [slo ,matrix,OCT ] = data_xml(adress )function [slo ,matrix,OCT ] = data_xml(adress )
% for i=1:
% A=BAGHERI.id_75661.id_133425.id_1248300.bscans{1,i}.img;
% j=2
   function [slo ,matrix,OCT ] = OCT_maker_read_file(j,B1)
% B1=strcat('DERAK01A_OS'); % 
DATA_A  = readbin(strcat(B1,'.octbin'));
A_cell = struct2cell(DATA_A);
A_cell2_1=struct2cell(A_cell{2,1});
A_cell2_2=struct2cell(A_cell2_1{2,1});
% Data1=A_cell2_2{2,1};
% Data2=A_cell2_2{2,2};
% j=2;
for i=1:size(A_cell2_2{j,1}.bscans,2)
   OCT_1(:,:,i)=A_cell2_2{j,1}.bscans{1,i}.img;
   
   slo=A_cell2_2{j,1}.slo.img;
   X_start(i)=A_cell2_2{j,1}.bscans{1,i}.data.start_mm.x;
   Y_start(i)=A_cell2_2{j,1}.bscans{1,i}.data.start_mm.y;
   X_stop(i)=A_cell2_2{j,1}.bscans{1,i}.data.end_mm.x;
   Y_stop(i)=A_cell2_2{j,1}.bscans{1,i}.data.end_mm.y;
   X_scale=A_cell2_2{j,1}.bscans{1,1}.data.scaleFactor.x;
      Y_scale=A_cell2_2{j,1}.bscans{1,1}.data.scaleFactor.y;
         Z_scale=A_cell2_2{j,1}.bscans{1,1}.data.scaleFactor.z;
   X_scale_slo=A_cell2_2{j,1}.slo.scaleFactor.x;
      Y_scale_slo=A_cell2_2{j,1}.slo.scaleFactor.y;
% startpoint(i,1)=(X_start(i)/X_scale_slo)+A_cell2_2{j,1}.slo.shift_px.x;
% startpoint(i,2)=(Y_start(i)/Y_scale_slo)+A_cell2_2{j,1}.slo.shift_px.y;
% stoppoint(i,1)= (X_stop(i)/X_scale_slo)+A_cell2_2{j,1}.slo.shift_px.x;
% stoppoint(i,2)=(Y_stop(i)/Y_scale_slo)+A_cell2_2{j,1}.slo.shift_px.y;
matrix(i,2)=(X_start(i)/X_scale_slo)+A_cell2_2{j,1}.slo.shift_px.x;% start_point
matrix(i,3)=(Y_start(i)/Y_scale_slo)+A_cell2_2{j,1}.slo.shift_px.y;% start point
matrix(i,4)=(X_stop(i)/X_scale_slo)+A_cell2_2{j,1}.slo.shift_px.x;% stop point
matrix(i,5)=(Y_stop(i)/Y_scale_slo)+A_cell2_2{j,1}.slo.shift_px.y;%stop point
%  center_point=[A_cell2_2{j,1}.slo.shift_px.x,A_cell2_2{j,1}.slo.shift_px.y];
center_point=[(matrix(i,2)+matrix(i,4))/2,(matrix(i,3)+matrix(i,5))/2];
matrix(i,7:8)=center_point;% center point
matrix(i,6)=((((matrix(i,2)-matrix(i,4))^2)+((matrix(i,3)-matrix(i,5))^2))^ 0.5); % radial
% figure(1),
% subplot(1,2,1), imshow(A_cell2_2{j,1}.bscans{1,i}.img,[]);
% subplot(1,2,2),imshow(A_cell2_2{j,1}.slo.img,[]);
% hold on
% slo1 = plot([startpoint(i,1),stoppoint(i,1)],[startpoint(i,2),stoppoint(i,2)],'color','y');
% pause
% end
OCT=OCT_1;
for ii=1
% for i=1:size(A_cell2_2{3,1}.bscans,2)
%    
%  OCT_2(:,:,i)=A_cell2_2{3,1}.bscans{1,i}.img;
% end
% for i=1:size(A_cell2_2{4,1}.bscans,2)
%    
%  OCT_3(:,:,i)=A_cell2_2{4,1}.bscans{1,i}.img;
% end
% for i=1:size(A_cell2_2{5,1}.bscans,2)
%    
%  OCT_4(:,:,i)=A_cell2_2{5,1}.bscans{1,i}.img;
% end
% for i=1:size(OCT_1,3)
%  figure(1),imshow(OCT_1(:,:,i),[]),pause,
% end;
end
end% x_1 =A_cell2_2{2,1}.bscans{1,4}.data.start_mm.x;
  end