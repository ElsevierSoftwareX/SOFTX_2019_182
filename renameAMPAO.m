function output = renameAMPAO(FileName,PathName,handles)
if FileName(end-2:end) == 'vol'
    [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (strcat(PathName,FileName));
    BScans(BScans>1e+38)=0;
    %     output.BScan = BScans.^.25;
    output.BScan = BScans;
    output.slo = slo;
    output.ScaleX= header.ScaleX;
    output.ScaleZ= header.ScaleZ;
    output.NumBScans=header.NumBScans;
    output.StartX = BScanHeader.StartX;
    output.EndX = BScanHeader.EndX;
    output.StartY = BScanHeader.StartY;
    output.EndY = BScanHeader.EndY;
    output.format = 'vol';
    
    output.PID = header.PID;
    output.SizeX = size(BScans,2);
    output.SizeZ = size(BScans,1);
    output.Distance = header.Distance;
    output.ExamTime = header.ExamTime;
    output.ScanPosition = header.ScanPosition;
    %     output.ScanType = header.ScanType;
    bd_pts = zeros(512,header.NumBScans,11);
    bd_pts(:,:,1)= BScanHeader.Boundary_1';
    bd_pts(:,:,2)= BScanHeader.Boundary_3';
    bd_pts(:,:,3)= BScanHeader.Boundary_4';
    bd_pts(:,:,4)= BScanHeader.Boundary_5';
    bd_pts(:,:,5)= BScanHeader.Boundary_6';
    bd_pts(:,:,6)= BScanHeader.Boundary_7';
    bd_pts(:,:,7)= BScanHeader.Boundary_9';
    bd_pts(:,:,8)= BScanHeader.Boundary_15';
    bd_pts(:,:,9)= BScanHeader.Boundary_16';
    bd_pts(:,:,10)= BScanHeader.Boundary_17';
    bd_pts(:,:,11)= BScanHeader.Boundary_2';
    output.bd_pts = bd_pts;
    output.headerThickness = header;
    %     output.sloRotation = 0;
elseif FileName(end-2:end) == 'xml'
    docNode = xmlread(strcat(PathName,FileName));
    data = docNode.getElementsByTagName('ExamURL');
    im_slo = char(data.item(0).getFirstChild.getData);
    loc = strfind(im_slo,'\') ;
    im_slo(1:loc(end)) =[];
    slo = imread(im_slo);
    i = 1;
    matrix = ones(19,5);
    for k = 1 : 19
        docNode = xmlread(strcat(PathName,FileName));
        data = docNode.getElementsByTagName('ExamURL');
        im_oct = char(data.item(k).getFirstChild.getData);
        loc = strfind(im_oct,'\') ;
        im_oct(1:loc(end)) =[];
        
        l = docNode.getElementsByTagName('ID');
        ID = char(l.item(k+3).getFirstChild.getData);
        matrix(k,1)=str2num(ID);
        x = docNode.getElementsByTagName('X');
        X1 = char(x.item(i).getFirstChild.getData);
        matrix(k,2)=str2num(X1);
        x1_1 = str2num(X1);
        y = docNode.getElementsByTagName('Y');
        y1 = char(y.item(i).getFirstChild.getData);
        matrix(k,3)=str2num(y1);
        y1_1 = str2num(y1);
        x = docNode.getElementsByTagName('X');
        X2 = char(x.item(i+1).getFirstChild.getData);
        matrix(k,4)=str2num(X2);
        x2_2 = str2num(X2);
        y = docNode.getElementsByTagName('Y');
        y2 = char(y.item(i+1).getFirstChild.getData);
        matrix(k,5)=str2num(y2);
        y2_2 = str2num(y2);
        dot(1,:,k)=[x1_1,x2_2,y1_1,y2_2];
        %  axes(handles.axes1);
        %  plot([x2_2/sclae_x,x1_1/sclae_x],[y2_2/sclae_x,y1_1/sclae_x],'color','y');
        % hold on;
        oct = imread(im_oct);
        %         [M,N,Q]=size(oct);
        %         OCT3D(:,k)=oct(:);
        % axes(handles.axes2);
        % imshow(oct)
        %         [M,N,Q]=size(oct);
        %         save oct oct
        OCT3D(:,:,k)=oct(:,:,1);
        % title(sprintf('ID = %i \n',k));
        %
        % pause(0.17);
        i = i+2;
        % hold off;
    end
    output.BScan = OCT3D;
    output.slo = slo;
    output.ScaleX= 0.0117;
    output.ScaleZ= 0.0117;
    output.NumBScans=19;
    output.StartX = dot(1,2,:);
    output.EndX = dot(1,1,:);
    output.StartY = dot(1,4,:);
    output.EndY = dot(1,3,:);
    output.format = 'xml';
elseif FileName(end-2:end) == 'img'
    [img_vol,vol_info] = octCirrusReader(filename);
    output.BScan = img_vol;
    %output.slo = slo;
    output.ScaleX= vol_info.vol_res(2);
    output.ScaleZ= vol_info.vol_res(1);
    output.NumBScans=size(img_vol,3);
    %     output.StartX = dot(1,2,:);
    %     output.EndX = dot(1,1,:);
    %     output.StartY = dot(1,4,:);
    %     output.EndY = dot(1,3,:);
    output.format = 'img';
    
    output.PID = vol_info.pid;
    output.SizeX = size(img_vol,2);
    output.SizeZ = size(img_vol,1);
    output.Distance = vol_info.vol_res(3);
    output.ExamTime = vol_info.scan_date;
    output.ScanPosition = vol_info.eye_side;
    output.ScanType = vol_info.scan_type;
    
    
    %             header.PID = vol_info.pid;
    %             header.SizeX = size(img_vol,2);
    %             header.NumBScans = ;
    %             header.SizeZ = size(img_vol,1);
    %             header.ScaleX = ;
    %             header.Distance = vol_info.vol_res(3);
    %             header.ScaleZ = ;
    %             header.ExamTime = vol_info.scan_date;
    %             header.ScanPosition = vol_info.eye_side;
    %             header.ScanType = vol_info.scan_type;
    
else %(for mat file)
    choice = questdlg('Is selected .mat produced in AMPAO?', ...
        'Mat Selection','Yes', 'No','Cancel');
    % Handle response
    switch choice
        case 'Yes'
            load (strcat(PathName,FileName));
            %         ss = load (strcat(PathName,FileName));
            %         name1 = fieldnames(ss);
            %         %     OCT3D = evalin('base', name1{1, 1});
            %         OCT3D = eval(name1{1, 1});
            %         output.BScan = OCT3D;
            %         output.NumBScans=size(OCT3D,3);
            output.format = 'matInside';
        case 'No'
            load (strcat(PathName,FileName));
            ss = load (strcat(PathName,FileName));
            name1 = fieldnames(ss);
            %     OCT3D = evalin('base', name1{1, 1});
            OCT3D = eval(name1{1, 1});
            output.BScan = OCT3D;
            output.NumBScans=size(OCT3D,3);
            output.format = 'matOutside';
        case 'Cancel'
    end
    
    
%     choice = choosedialog_Mat;
%     if strcmp(choice,'Yes')
%         load (strcat(PathName,FileName));
%         %         ss = load (strcat(PathName,FileName));
%         %         name1 = fieldnames(ss);
%         %         %     OCT3D = evalin('base', name1{1, 1});
%         %         OCT3D = eval(name1{1, 1});
%         %         output.BScan = OCT3D;
%         %         output.NumBScans=size(OCT3D,3);
%         output.format = 'matInside';
%         
%     else
%         
%         load (strcat(PathName,FileName));
%         ss = load (strcat(PathName,FileName));
%         name1 = fieldnames(ss);
%         %     OCT3D = evalin('base', name1{1, 1});
%         OCT3D = eval(name1{1, 1});
%         output.BScan = OCT3D;
%         output.NumBScans=size(OCT3D,3);
%         output.format = 'matOutside';
%     end
end
% handles.data = output;
% output