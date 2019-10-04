%%%%%%%%%%%%%%%%%%%%%%%%%%
%(c)  Rahele Kafieh 2011
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% THE MODIFICATIONS ARE DONE ON 2012/12
function output = curvCorrectionAMPAO(handles)

%------------- Get Slice Numbers
prompt={'Enter value \range of slices to be aligned'};
name = 'Scan Number(s)';
defaultans = {'1'};
options.Interpreter = 'none';
answer = inputdlg(prompt,name,[1 40],defaultans,options);
scan_num=str2num(answer{1,1});
%-------------
if strcmp(handles.data.format , 'vol')
    data = double(handles.data.BScan).^.25;
else 
    data = double(handles.data.BScan);
end

for counter = scan_num
    
    this_scan = counter

    lImage1 = data(:,:,counter);
    lImage1 = imrotate(lImage1,-90);

    [mm,nn] = size(lImage1);
    X = zeros(mm,1);
    width = 10;
    smoothii=100;
    device = 'HEIDELBERG';
    vert_change = 10;
    % end
    
    x_old =0;
    for count = 1:mm./width
        lImage = lImage1((count-1)*width+1:count*width,:);
        
        a=imrotate(lImage, 90);
        %     figure, imagesc(a),colormap('gray')
        fun = @(x) mean(mean(x));
        I = blkproc(a,[6 6],fun);
        j = imresize(I, size(a));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%
        im = imrotate(j, 90);
        
        im =- double(im);
        [m,n] = size(im);
        c = zeros(m,n);
        p = zeros(m,n,'int8');  % save memory by using 8bit integers
        c(1,:) = im(1,:);
        
        for i = 2:m
            c0 = c(i-1,:);
            d = repmat( im(i,:), 3, 1 ) + [c0; c0(2:end) c0(end); c0(1) c0(1:(end-1))];
            [c(i,:),p(i,:)] = min(d);
        end
        
        x = zeros(m,1);
        [cost,xpos] = min( c(m,:) );
        for i = m:-1:2
            x(i) = xpos;
            if p(i,xpos)==2 && xpos<n
                xpos = xpos+1;
            elseif p(i,xpos)==3 && xpos>1
                xpos = xpos-1;
            end
        end
        x(1) = xpos;
        
        x=x(end:-1:1);
        x = smooth(x,0.3,'rloess');
        if count~=1
            if abs(mean(x)-mean(x_old))>vert_change
                x= x_old;
            end
        end
        x_old =x;
        X((count-1)*width+1:count*width)=x;
    end
    X1 = smooth(X,smoothii,'rloess');
    A=imrotate(lImage1, 90);
    
    a=imrotate(lImage1, 90);
    %     figure, imagesc(a),colormap('gray')
    fun = @(x) mean(mean(x));
    I = blkproc(a,[6 6],fun);
    j = imresize(I, size(a));
    % j=a;
    %     figure; imagesc(j)   ,colormap('gray') ;
    %  subplot 211,imshow(a,[]),title('original image')
    %  subplot 212,imshow(j,[]),title('low image')
    %%%%%%%%%%%%%%%%%%%%%%%%%
    im = imrotate(j, 90);
    im =- double(im);
    [m,n] = size(im);
    c = zeros(m,n);
    p = zeros(m,n,'int8');  % save memory by using 8bit integers
    c(1,:) = im(1,:);
    
    
    for i = 2:m
        c0 = c(i-1,:);
        d = repmat( im(i,:), 3, 1 ) + [c0; c0(2:end) c0(end); c0(1) c0(1:(end-1))];
        [c(i,:),p(i,:)] = min(d);
    end
    
    x = zeros(m,1);
    [cost,xpos] = min( c(m,:) );
    for i = m:-1:2
        x(i) = xpos;
        if p(i,xpos)==2 && xpos<n
            xpos = xpos+1;
        elseif p(i,xpos)==3 && xpos>1
            xpos = xpos-1;
        end
    end
    x(1) = xpos;
    x=x(end:-1:1);
    X3 = x;
    
    %get seed starts
    %%%%%%%%%%%%%%%%%%
    
    X2 = smooth(X3,0.3,'rloess');
    
    %%%%%%%%%%%%%
    % Alignment
    X_max = max(X2)+1;
    Im_aligned =zeros(size(A));
    diff = X_max - X2;
    for i=1:mm
        leni = A(:,i);
        leni(floor(diff(i)):nn) = leni(1:nn-floor(diff(i))+1);
        leni(1:floor(diff(i))-1) = leni(nn-floor(diff(i)):nn-2);
        Im_aligned(:,i) = leni;
    end
    Im_alig(:,:,counter)=Im_aligned;
end 
output.BScan = Im_alig;
output.NumBScans = scan_num;
% handles.data.PreProcess = output;
% guidata(hObject,handles);

uisave('output')
