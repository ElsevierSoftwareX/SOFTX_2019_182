function  VOLWrite(address)
% clc;
% clear all;
% close all;

path =  'voldata.vol';
fid = fopen(path);
visu=1;
 
%--------------------------------------------------------------------------
% File header read

header.Version = fread( fid, 12, '*int8' );
header.SizeX = fread( fid, 1, '*int32' );
header.NumBScans = fread( fid, 1, '*int32' );
header.SizeZ = fread( fid, 1, '*int32' );
header.ScaleX = fread( fid, 1, '*double' );
header.Distance = fread( fid, 1, '*double' );
header.ScaleZ = fread( fid, 1, '*double' );
header.SizeXSlo = fread( fid, 1, '*int32' );
header.SizeYSlo = fread( fid, 1, '*int32' );
header.ScaleXSlo = fread( fid, 1, '*double' );
header.ScaleYSlo = fread( fid, 1, '*double' );
header.FieldSizeSlo = fread( fid, 1, '*int32' );
header.ScanFocus = fread( fid, 1, '*double' );
header.ScanPosition = char(fread( fid, 4, '*uchar' )');
header.ExamTime = fread( fid, 1, '*int64' );
header.ScanPattern = fread( fid, 1, '*int32' );
header.BScanHdrSize = fread( fid, 1, '*int32' );
header.ID = char(fread( fid, 16, '*uchar' )');
header.ReferenceID = char(fread( fid, 16, '*uchar' )');
header.PID = fread( fid, 1, '*int32' );
header.PatientID = char(fread( fid, 21, '*uchar' )');
header.Padding = fread( fid, 3, '*int8' );
header.DOB = fread( fid, 1, '*double' );
header.VID = fread( fid, 1, '*int32' );
header.VisitID = char(fread( fid, 24, '*uchar' )');
header.VisitDate = fread( fid, 1, '*double' );
header.GridType = fread( fid, 1, '*int32');
header.GridOffset = fread( fid, 1, '*int32');
header.Spare = fread( fid, 1832, '*int8' );
 

    disp(['---------------------------------------------']);
    disp(['           Version: ' char(header.Version')]);
    disp(['             SizeX: ' num2str(header.SizeX)]);
    disp(['         NumBScans: ' num2str(header.NumBScans)]);
    disp(['             SizeZ: ' num2str(header.SizeZ)]);
    disp(['            ScaleX: ' num2str(header.ScaleX) ' mm']);
    disp(['          Distance: ' num2str(header.Distance) ' mm']);
    disp(['            ScaleZ: ' num2str(header.ScaleZ) ' mm']);
    disp(['          SizeXSlo: ' num2str(header.SizeXSlo)]);
    disp(['          SizeYSlo: ' num2str(header.SizeYSlo)]);
    disp(['         ScaleXSlo: ' num2str(header.ScaleXSlo) ' mm']);
    disp(['         ScaleYSlo: ' num2str(header.ScaleYSlo) ' mm']);
    disp(['FieldSizeSlo (FOV): ' num2str(header.FieldSizeSlo) 'deg']);
    disp(['         ScanFocus: ' num2str(header.ScanFocus) ' dpt']);
    disp(['      ScanPosition: ' char(header.ScanPosition)]);
    disp(['          ExamTime: ' datestr(double(header.ExamTime(1)/(1e7*60*60*24)+584755+(2/24)))]);
    disp(['       ScanPattern: ' num2str(header.ScanPattern)]);
    disp(['      BScanHdrSize: ' num2str(header.BScanHdrSize) ' bytes']);
    disp(['                ID: ' char(header.ID)]);
    disp(['       ReferenceID: ' char(header.ReferenceID)]);
    disp(['               PID: ' num2str(header.PID)]);
    disp(['         PatientID: ' char(header.PatientID)]);
    disp(['               DOB: ' datestr(header.DOB+693960)]);
    disp(['               VID: ' num2str(header.VID)]);
    disp(['           VisitID: ' char(header.VisitID)]);
    disp(['         VisitDate: ' datestr(double(header.VisitDate+693960))]);
    disp(['          GridType: ' num2str(header.GridType)]);
    disp(['        GridOffset: ' num2str(header.GridOffset)]);
    disp(['---------------------------------------------']);


status = fseek( fid, 2048, -1 );
 
%--------------------------------------------------------------------------
% SLO image read

mfile = memmapfile('voldata.vol','format','uint8','writable',true);

docNode = xmlread(address);
data = docNode.getElementsByTagName('ExamURL');
im_slo = char(data.item(0).getFirstChild.getData);
loc = strfind(im_slo,'\') ; 
im_slo(1:loc(end)) =[];
SLO = rgb2gray(imread(im_slo));
SLO = SLO(:);

mfile.Data(2048:591871) = SLO;

    slo = fread(fid, header.SizeXSlo*header.SizeYSlo, '*uint8');
    slo = reshape(slo, header.SizeXSlo, header.SizeYSlo);
    slo = imrotate(slo, 90);
    slo = flipud(slo);
    slo = slo';
 
%--------------------------------------------------------------------------
% Display the image

if visu==1
    scrsz = get(0,'ScreenSize');
    figure('Position',[1 0 scrsz(3) scrsz(4)-70]);
    subplot(1,2,1);
    
    imshow(slo);
end
 
%--------------------------------------------------------------------------
% B-scans read (including image data and scanning position)

status = fseek( fid, 2048+(header.SizeXSlo*header.SizeYSlo), -1 );
 
    BScans=zeros(header.SizeZ, header.SizeX ,header.NumBScans, 'single');

BScanHeader.StartX = zeros(1, header.NumBScans, 'double');
BScanHeader.StartY = zeros(1, header.NumBScans, 'double');
BScanHeader.EndX = zeros(1, header.NumBScans, 'double');
BScanHeader.EndY = zeros(1, header.NumBScans, 'double');
BScanHeader.NumSeg = zeros(1, header.NumBScans, 'int32');
BScanHeader.Quality = zeros(1, header.NumBScans, 'single');
BScanHeader.Shift = zeros(1, header.NumBScans, 'int32');
BScanHeader.ILM = zeros(header.NumBScans,header.SizeZ, 'single');
BScanHeader.RPE = zeros(header.NumBScans,header.SizeZ, 'single');
BScanHeader.NFL = zeros(header.NumBScans,header.SizeZ, 'single');
 
fid1 = fopen(path);
for zz = 1:header.NumBScans
    ii = zz - 1;
    
     data = docNode.getElementsByTagName('ExamURL');
     im_oct = char(data.item(zz).getFirstChild.getData);
     loc = strfind(im_oct,'\') ; 
     im_oct(1:loc(end)) =[];
     oct = double(imread(im_oct));

     oct = oct(:,:,1)';
     
     oct = oct(:); 

    %----------------------------------------------------------------------
    %BScan position in SLO image
        
    status = fseek( fid1, 16+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    StartX = fread( fid1, 1, '*double' ); 
    StartY = fread( fid1, 1, '*double' );  
    EndX = fread( fid1, 1, '*double' );
    EndY = fread( fid1, 1, '*double' );  
    NumSeg = fread( fid1, 1, '*int32' );
    OffSeg = fread( fid1, 1, '*int32' );
    Quality = fread( fid1, 1, '*float32' );
    Shift = fread( fid1, 1, '*int32' );
    
    BScanHeader.StartX(zz) = StartX;
    BScanHeader.StartY(zz) = StartY;
    BScanHeader.EndX(zz) = EndX;
    BScanHeader.EndY(zz) = EndY;
    BScanHeader.NumSeg(zz) = NumSeg;
    BScanHeader.Shift(zz) = Shift;
    BScanHeader.Quality(zz) = Quality;
 
    %----------------------------------------------------------------------
    % Display the images
    
    if visu==1
        StartX_px = round(StartX/header.ScaleXSlo); %StartX in pixels
        StartY_px = round(StartY/header.ScaleYSlo); %StartY in pixels
        EndX_px = round(EndX/header.ScaleXSlo); %EndX in pixels
        EndY_px = round(EndY/header.ScaleYSlo); %EndY in pixels
        subplot(1,2,1); hold on, line([StartX_px EndX_px],[StartY_px EndY_px],'color','b');
    end
    
    %----------------------------------------------------------------------
    
BsBlkOffs=header.BScanHdrSize+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4));
    mfile = memmapfile('voldata.vol','Offset',BsBlkOffs,'format','single','writable',true);
    mfile.Data(1 : 1+253951) = oct;
      
    % BScan reading
    
        status = fseek( fid, header.BScanHdrSize+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
        oct = fread( fid, header.SizeX*header.SizeZ, '*float32' );
        oct = reshape( oct, header.SizeX, header.SizeZ );
        BScans(:,:,zz)=oct';
     
    %----------------------------------------------------------------------
    % Display the images
    
    if visu==1
        subplot(1,2,2);
        imshow(BScans(:,:,zz),[]);
        drawnow
    end
 
    %----------------------------------------------------------------------
    % Segmentation reading
    
    status = fseek( fid, 256+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    seg = (fread( fid, NumSeg*header.SizeZ, '*float' ))';
    BScanHeader.ILM(zz,:) = seg(1:header.SizeZ);
    BScanHeader.RPE(zz,:) = seg(header.SizeZ+1:header.SizeZ*2);
    if NumSeg == 3
        BScanHeader.NFL(zz,:) = seg(header.SizeZ*2+1:header.SizeZ*3);
    end
end

%----------------------------------------------------------------------
% This is for ILM and RPE INVALID values interpolation using median of 3x3 
% neigborhood around INVALID value.

BScanHeader.ILM(BScanHeader.ILM>1e6) = nan;
BScanHeader.RPE(BScanHeader.RPE>1e6) = nan;
if NumSeg == 3
    BScanHeader.NFL(BScanHeader.NFL>1e6) = nan;
end
fclose( fid );

end