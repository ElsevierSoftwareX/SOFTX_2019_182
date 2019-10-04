function conversion(FileName,PathName,handles)

uu = get(handles.conversion,'Value');
if uu==1
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
        % axes(handles.axes2);
        % imshow(oct)
        %         [M,N,Q]=size(oct);
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
elseif uu == 2
%     VOLWrite(strcat(PathName,FileName));
%     %     Volname = uiputfile('voldata.vol','Save file name')
%     %     save(Volname)
%     %     handles.PathName = PathName;
%     %     handles.FileName = FileName;
% else
    [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (strcat(PathName,FileName));
    BScans(BScans>1e+38)=0;
    output.BScan = BScans.^.25;
    output.slo = slo;
    output.ScaleX= header.ScaleX;
    output.ScaleZ= header.ScaleZ;
    output.NumBScans=header.NumBScans;
    output.StartX = BScanHeader.StartX;
    output.EndX = BScanHeader.EndX;
    output.StartY = BScanHeader.StartY;
    output.EndY = BScanHeader.EndY;
    
    
end
handles.data = output;
uisave('output')