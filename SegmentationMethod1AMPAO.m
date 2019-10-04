function  SegmentationMethod1AMPAO(handles)

folderPath = handles.PathName;
% folderPath = 'D:\kafieh\AMPAO';
filterindex = 2;
% folderPath = strcat(folderPath,'\');

image = double(handles.data.BScan);

%------------- Get Slice Numbers
prompt={'Enter value \range of slices to be Segmented (at least 3 slices)'};
name = 'Scan Number(s)';
defaultans = {'1:3'};
options.Interpreter = 'none';
answer = inputdlg(prompt,name,[1 40],defaultans,options);
scan_num=str2num(answer{1,1});

%------------- convert to tiff
 counti = 0;
    for i=scan_num
        counti = counti +1;
        mat2tiff(image(:,:,i),folderPath, strcat(handles.FileName(1:end-4),num2str(i)))
        filename{counti} = strcat(handles.FileName(1:end-4),num2str(i),'.tiff');
    end

%% --------------------------------------Layer Segmentation Casserer


 counti = 0;
 for i = scan_num
     counti = counti +1;
     imagePath{counti} = [folderPath ,filename{counti}];
     
 end

% figure;imshow(imread(imagePath{1}),[])
figure;imagesc(imread(imagePath{1})), colormap gray

title('pick a region of interest to segment for the selected images');
[trsh rect] = imcrop;
xrange = round(rect(1)):round(rect(1)+rect(3));
yrange = round(rect(2)):round(rect(2)+rect(4));



%% Section 2, automatically segments the retinal layers based on graph theory.

for i = 1:numel(imagePath)
    
    display(sprintf('segmenting image %d of %d',i,numel(imagePath)));
    
    % read in the image.
    img = imread(imagePath{i});
    
    % error checking, get one channel from image.
    if size(img,3) > 1
        img = img(:,:,1);
        display('warning: this is probably not an oct image');
    end
    
    % make image type as double.
    img = double(img);
    
    % get size of image.
    szImg = size(img);
    
    %segment whole image if yrange/xrange is not specified.
    if isempty(yrange) && isempty(xrange)
        yrange = 1:szImg(1);
        xrange = 1:szImg(2);
    end
    img = img(yrange,xrange);
    
    % get retinal layers.
    [retinalLayers, params] = getRetinalLayers(img);
    
    % save range of image.
    params.yrange = yrange;
    params.xrange = xrange;
    
    % save data to struct.
    imageLayer(i).imagePath = imagePath{i};
    imageLayer(i).retinalLayers = retinalLayers;
    imageLayer(i).params = params;
    
end

% save segmentation
filename = [imageLayer(1).imagePath '_octSegmentation.mat'];
save(filename, 'imageLayer');
display(sprintf('segmentation saved to %s',filename));

% %%   Section 3, using a GUI, iterate through the segmentation results,
% %              and maually or semi-automatically correct the segmented
% %              retainl layers.
% 
% close all;
% 
% filename = [imagePath{1} '_octSegmentation.mat'];
% 
% isReviewSegmentation = 1;
% if isReviewSegmentation,
%     [h,guiParam] = octSegmentationGUI(filename);
%     
%     if guiParam.proceed
%         delete(guiParam.figureOCT);
%         delete(h);
%     else
%         return;
%     end
% end
% 
% 
% %%  Section 4, calculate and print out retinal thickness (in pixels)
% 
uisave( 'imageLayer','SegmentedData')

calculateRetinalThickness(imagePath) %scan_num is scans (more than 3) selected by user for segmentation
% 
  %% ------------------------------------- 
%     Im_seg(:,:,count)=imfi;
% end
% output.BScan = Im_den;
% output.NumBScans = scan_num;
% % handles.data.PreProcess = output;
% guidata(hObject,handles);

% uisave('output')
% save Im_den