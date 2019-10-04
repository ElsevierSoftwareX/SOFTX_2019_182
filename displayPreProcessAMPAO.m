function displayPreProcessAMPAO(handles)

output_PreProcess = handles.data.PreProcess;
output = handles.data;

cla(handles.axes1,'reset');
cla(handles.axes2,'reset');

for k =  output_PreProcess.NumBScans
    cla(handles.axes1)
    axes(handles.axes1);
    if strcmp(output.format , 'vol')
        imshow(output.BScan(:,:,k).^.25,[]);
    else
        imshow(output.BScan(:,:,k),[]);
    end
    title(sprintf('Original Image %i \n',k));
    
    cla(handles.axes2)
    axes(handles.axes2);
    imshow(output_PreProcess.BScan(:,:,k),[]);
    title(sprintf('Processed Image %i \n',k));
    
    pause(1);
end

% output.BScan = Im_den;
% output.NumBScans = count;
%
% if output.format == 'vol'
%     cla(handles.axes1)
%     axes(handles.axes1);
%     imshow(output.slo,[]);
%     hold on;
%     for k = 1 : output.NumBScans
%         axes(handles.axes1);
%         plot([output.StartX(k) /output.ScaleX,output.EndX(k)/output.ScaleX],[output.StartY(k)/output.ScaleX,output.EndY(k)/output.ScaleX],'color','y');
%         hold on;
%         %         oct=OCT3D(:,k);
%         %         oct=reshape(oct,M,N,Q);
%         axes(handles.axes2);
%         imshow(output.BScan(:,:,k).^.25);
%         title(sprintf('ID = %i \n',k));
%
%         pause(0.17);
%         hold off;
%     end
% elseif output.format == 'xml'
%
%     output = handles.data;
%
%     axes(handles.axes1);
%     imshow(output.slo,[]);
%     hold on;
%     for k = 1 : output.NumBScans
%         axes(handles.axes1);
%         plot([output.StartX(k) /output.ScaleX,output.EndX(k)/output.ScaleX],[output.StartY(k)/output.ScaleX,output.EndY(k)/output.ScaleX],'color','y');
%         hold on;
%         %         oct=OCT3D(:,k);
%         %         oct=reshape(oct,M,N,Q);
%         axes(handles.axes2);
%         imshow(output.BScan(:,:,k));
%         title(sprintf('ID = %i \n',k));
%
%         pause(0.17);
%         hold off;
%     end
% else %mat format
%
%     output = handles.data;
%
%
%     for k = 1 : output.NumBScans
%         axes(handles.axes1);
%         imshow(zeros(512,512),[]);
%         hold on
%         axes(handles.axes2);
%         imshow(output.BScan(:,:,k),[]);
%         title(sprintf('ID = %i \n',k));
%
%         pause(0.17);
%         hold off;
%     end
% end