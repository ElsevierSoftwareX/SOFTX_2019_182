function octandsloshow(PathName,FileName, handles)

if FileName(end-2:end) == 'vol'
    [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (strcat(PathName,FileName));
    %     OCT3D = BScans;
    sclae_x = header.ScaleX;
    BScans(BScans>1e+38)=0;
    num=header.NumBScans;
    axes(handles.axes1);
    imshow(slo,[]);
    hold on;
    for k = 1 : num
        axes(handles.axes1);
        plot([BScanHeader.StartX(k) /sclae_x,BScanHeader.EndX(k)/sclae_x],[BScanHeader.StartY(k)/sclae_x,BScanHeader.EndY(k)/sclae_x],'color','y');
        hold on;
        %         oct=OCT3D(:,k);
        %         oct=reshape(oct,M,N,Q);
        axes(handles.axes2);
        imshow(BScans(:,:,k).^.25);
        title(sprintf('ID = %i \n',k));
        
        pause(0.17);
        hold off;
    end
    
elseif FileName(end-2:end) == 'xml'
    docNode = xmlread(strcat(PathName,FileName));
    i = 1;
    matrix = ones(19,5);
    data = docNode.getElementsByTagName('ExamURL');
    im_slo = char(data.item(0).getFirstChild.getData);
    loc = strfind(im_slo,'\') ;
    im_slo(1:loc(end)) =[];
    slo = imread(im_slo);
    
    % axes(handles.axes1);
    % imshow(slo,[]);
    sclae_x = 0.0117;
    % hold on;
    
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
        % axes(handles.axes2);
        % imshow(oct)
        [M,N,Q]=size(oct);
        OCT3D(:,k)=oct(:);
        
        % title(sprintf('ID = %i \n',k));
        %
        % pause(0.17);
        i = i+2;
        % hold off;
    end
    
    
    
    
    num=19;
    axes(handles.axes1);
    imshow(slo,[]);
    hold on;
    for k = 1 : num
        axes(handles.axes1);
        plot([dot(1,2,k)/sclae_x,dot(1,1,k)/sclae_x],[dot(1,4,k)/sclae_x,dot(1,3,k)/sclae_x],'color','y');
        hold on;
        oct=OCT3D(:,k);
        oct=reshape(oct,M,N,Q);
        axes(handles.axes2);
        imshow(oct);
        title(sprintf('ID = %i \n',k));
        
        pause(0.17);
        hold off;
    end
else %(for mat file)
    ss = load (strcat(PathName,FileName));
    name1 = fieldnames(ss);
    OCT3D = evalin('base', name1{1, 1});
    for k = 1 : size(OCT3D,3)
        oct=OCT3D(:,:,k);
        axes(handles.axes2);
        imshow(oct,[]);
        title(sprintf('ID = %i \n',k));        
        pause(0.17);
        hold off;
    end
end