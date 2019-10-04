function displayThicknessMapAMPAO(boundaryPoints,header,circ_rad,showgrid,handles)
%
output = handles.data;
cla(handles.axes1,'reset');
cla(handles.axes2,'reset');

if strcmp(output.format , 'vol')
    thicknessMapAMPAO2(boundaryPoints,header,circ_rad,showgrid,handles)
    colormap winter
    
    cla(handles.axes1)
    axes(handles.axes1);
    imshow(output.slo,[]);
    hold on;
    for k = 1 : output.NumBScans
        axes(handles.axes1);
        plot([output.StartX(k) /output.ScaleX,output.EndX(k)/output.ScaleX],[output.StartY(k)/output.ScaleX,output.EndY(k)/output.ScaleX],'color','y');
        hold on;
        %         oct=OCT3D(:,k);
        %         oct=reshape(oct,M,N,Q);
        %         axes(handles.axes2);
        %         imshow(output.BScan(:,:,k).^.25,[]);
        %         title(sprintf('ID = %i \n',k));
        
        pause(0.17);
        colormap winter
        
    end
    if ~strncmp(header.ScanPosition,'OD',2)
        warndlg('ETDRS is fliped in left eye for consistency of "nasal" to the right side')
    end
    %
    % elseif strcmp(output.format,'xml')
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
    %         imshow(output.BScan(:,:,k),[]);
    %         title(sprintf('ID = %i \n',k));
    %
    %         pause(0.17);
    %         hold off;
    %     end
    % else %mat format
    %
    %     if strcmp(output.format,'matInside')
    %         output = handles.data;
    %         cla(handles.axes1)
    % %         axes(handles.axes1);
    % %         imshow(output.slo,[]);
    % %         hold on;
    %         for k = 1 : output.NumBScans
    % %             axes(handles.axes1);
    % %             plot([output.StartX(k) /output.ScaleX,output.EndX(k)/output.ScaleX],[output.StartY(k)/output.ScaleX,output.EndY(k)/output.ScaleX],'color','y');
    % %             hold on;
    %             %         oct=OCT3D(:,k);
    %             %         oct=reshape(oct,M,N,Q);
    %             axes(handles.axes2);
    %             imshow(output.BScan(:,:,k),[]);
    %             title(sprintf('ID = %i \n',k));
    %
    %             pause(0.17);
    %             hold off;
    %         end
    %
    %     else
    %         output = handles.data;
    %
    %         for k = 1 : output.NumBScans
    %             axes(handles.axes1);
    %             imshow(zeros(512,512),[]);
    %             hold on
    %             axes(handles.axes2);
    %             imshow(output.BScan(:,:,k),[]);
    %             title(sprintf('ID = %i \n',k));
    %
    %             pause(0.17);
    %             hold off;
    %         end
    %     end
    
end