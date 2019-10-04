function [header, BScanHeader, slo, BScans, ThicknessGrid] = open_vol (path, options)
%OPENVOL Read Heidelberg Engineering (HE) OCT raw files (VOL ending)
% [HEADER, BSCANHEADER, SLO, BSCANS] = OPENVOL(PATH, OPTIONS)
% This function performs volume OCT data (xxxx.vol) reading. 
% HEADER: Header information as described by HE. Struct with each entry
%   named by the HE conventions. 
% BSCANHEADER: B-scan header information. Struct with each entry named by 
%   the HE conventions. The entries consist of vectors/matrices, with one 
%   field per B-Scan. 
% SLO: Slo image as unsigned integers.
% BScans: BScans. 3D matrix with floating point data of the B-Scans.
% PATH: Filename of the VOL-file to read, with ending.
%
% OPTIONS: Various output possibilites, 
%   written in the options variable as string text, i.e. 'visu writemeta'
%   Possibilities: 
%       visu: The read-in is visualized.
%       visuseg: The HE segmentation is also visualized.
%       metawrite: A metafile with the header and B-scan header information 
%           is written out (if it does not exist already)
%       metawriteforce: An existing metafile is replaced
%       header: Only the Header and BScansHeader Information is read, 
%            not the image data  
%       nodisp: nothing is diplayed during read in
%
% Originally writen by Radim Kolar, Brno University, Czech Republic
% Modified by Markus Mayer, Pattern Recognition Lab, University of
% Erlangen-Nuremberg
%
% The segmented data and thickness grid reading parts have been expanded by
% Seyedamirhosein Motamedi, Charite Universitätmedizin Berlin, August 2016
% 
% Currently up to Version HSF-XXX-102 is supported, without Thickness Grid
%
% First final Version: March 2010
 
%% If only one argument is defined
if nargin==1,
    options = '';
end
%% Open file
visu = 0;
visuseg = 0;

if numel(strfind(options, 'visu')) ~= 0
    visu = 1;
end

if numel(strfind(options, 'visuseg')) ~= 0
    display('hello');
    visuseg = 1;
end

fid = fopen(path);
 
%% File header read
header.Version = char(fread( fid, 12, 'uchar' ));
header.SizeX = fread( fid, 1, 'int32' );
header.NumBScans = fread( fid, 1, 'int32' );
header.SizeZ = fread( fid, 1, 'int32' );
header.ScaleX = fread( fid, 1, 'double' );
header.Distance = fread( fid, 1, 'double' );
header.ScaleZ = fread( fid, 1, 'double' );
header.SizeXSlo = fread( fid, 1, 'int32' );
header.SizeYSlo = fread( fid, 1, 'int32' );
header.ScaleXSlo = fread( fid, 1, 'double' );
header.ScaleYSlo = fread( fid, 1, 'double' );
header.FieldSizeSlo = fread( fid, 1, 'int32' );
header.ScanFocus = fread( fid, 1, 'double' );
header.ScanPosition = char(fread( fid, 4, 'uchar' )');
header.ExamTime = datestr(fread( fid, 1, 'uint64' )/(1e7*60*60*24) + datenum('1 January 1601 00:00:00'));
header.ScanPattern = fread( fid, 1, 'int32' );
header.BScanHdrSize = fread( fid, 1, 'int32' );
header.ID = char(fread( fid, 16, 'uchar' )');
header.ReferenceID = char(fread( fid, 16, 'uchar' )');
header.PID = fread( fid, 1, 'int32' );
header.PatientID = char(fread( fid, 21, 'uchar' )');
header.Padding = fread( fid, 3, 'int8' );
header.DOB = datestr(fread( fid, 1, 'double' ) + datenum('30 December 1899 00:00:00'));
header.VID = fread( fid, 1, 'int32' );
header.VisitID = char(fread( fid, 24, 'uchar' )');
header.VisitDate = datestr(fread( fid, 1, 'double' ) + datenum('30 December 1899 00:00:00'));
header.GridType = fread( fid, 1, 'int32');
header.GridOffset = fread( fid, 1, 'int32');
header.Spare = fread( fid, 1832, 'int8' );
 
% if numel(strfind(options, 'nodisp')) == 0
%     disp(['---------------------------------------------']);
%     disp(['           Version: ' char(header.Version')]);
%     disp(['             SizeX: ' num2str(header.SizeX)]);
%     disp(['         NumBScans: ' num2str(header.NumBScans)]);
%     disp(['             SizeZ: ' num2str(header.SizeZ)]);
%     disp(['            ScaleX: ' num2str(header.ScaleX) ' mm']);
%     disp(['          Distance: ' num2str(header.Distance) ' mm']);
%     disp(['            ScaleZ: ' num2str(header.ScaleZ) ' mm']);
%     disp(['          SizeXSlo: ' num2str(header.SizeXSlo)]);
%     disp(['          SizeYSlo: ' num2str(header.SizeYSlo)]);
%     disp(['         ScaleXSlo: ' num2str(header.ScaleXSlo) ' mm']);
%     disp(['         ScaleYSlo: ' num2str(header.ScaleYSlo) ' mm']);
%     disp(['FieldSizeSlo (FOV): ' num2str(header.FieldSizeSlo) '???']);
%     disp(['         ScanFocus: ' num2str(header.ScanFocus) ' dpt']);
%     disp(['      ScanPosition: ' char(header.ScanPosition)]);
%     % OD - right eye; OS - left eye
%     % disp([' ExamTime: ' datestr(header.ExamTime/(1e7*60*60*24)+1601*365.24 - 0.24 + 5 + 1)]);
%     disp(['       ScanPattern: ' num2str(header.ScanPattern)]);
%     disp(['      BScanHdrSize: ' num2str(header.BScanHdrSize) ' bytes']);
%     disp(['                ID: ' char(header.ID)]);
%     disp(['       ReferenceID: ' char(header.ReferenceID)]);
%     disp(['               PID: ' num2str(header.PID)]);
%     disp(['         PatientID: ' char(header.PatientID)]);
%     disp(['               DOB: ' datestr(header.DOB+693975)]);
%     disp(['               VID: ' num2str(header.VID)]);
%     disp(['           VisitID: ' char(header.VisitID)]);
%     disp(['         VisitDate: ' datestr(header.VisitDate+693975)]);
%     disp(['          GridType: ' num2str(header.GridType)]);
%     disp(['        GridOffset: ' num2str(header.GridOffset)]);
%     disp(['---------------------------------------------']);
% end

status = fseek( fid, 2048, -1 );
 
%% SLO image read

if(numel(strfind(options, 'header')) == 0)
    slo = fread( fid, header.SizeXSlo*header.SizeYSlo, 'uint8' );
    slo = reshape( slo, header.SizeXSlo, header.SizeYSlo );
    slo = imrotate( slo, 90 );
    slo = flipud( slo );
else
    slo = [];
end
scrsz = get(0,'ScreenSize');
 
%% Display the image
if visu==1
    figure('Position',[1 0 scrsz(3) scrsz(4)-70]);
    clf;
    subplot(1,2,1);
    imshow(slo,[]);
end
 
%% B-scans and segmentation read (including image data and scanning position)
status = fseek( fid, 2048+(header.SizeXSlo*header.SizeYSlo), -1 );
 
if(numel(strfind(options, 'header')) == 0)
    BScans=zeros(header.SizeZ, header.SizeX ,header.NumBScans);
else
    BScans= [];
end

BScanHeader.StartX = zeros(1, header.NumBScans);
BScanHeader.StartY = zeros(1, header.NumBScans);
BScanHeader.EndX = zeros(1, header.NumBScans);
BScanHeader.EndY = zeros(1, header.NumBScans);
BScanHeader.NumSeg = zeros(1, header.NumBScans);
BScanHeader.Quality = zeros(1, header.NumBScans);
BScanHeader.Shift = zeros(1, header.NumBScans);

for zz = 1:header.NumBScans

    ii = zz - 1;
    
    %Evaluation of BScan position in SLO image
    status = fseek( fid, 2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    Version = char(fread( fid, 12, 'uchar' ));
    BScanHdrSize = fread( fid, 1, 'int32' );
    StartX = fread( fid, 1, 'double' ); 
    StartY = fread( fid, 1, 'double' );  
    EndX = fread( fid, 1, 'double' );
    EndY = fread( fid, 1, 'double' );  
    NumSeg = fread( fid, 1, 'int32' );
    OffSeg = fread( fid, 1, 'int32' );
    Quality = fread( fid, 1, 'float32' );
    Shift = fread( fid, 1, 'int32' );
    Spare = fread(fid, 192, 'int8');
    
    BScanHeader.Version(:, zz) = Version;
    BScanHeader.BScanHdrSize(zz) = BScanHdrSize;
    BScanHeader.StartX(zz) = StartX;
    BScanHeader.StartY(zz) = StartY;
    BScanHeader.EndX(zz) = EndX;
    BScanHeader.EndY(zz) = EndY;
    BScanHeader.NumSeg(zz) = NumSeg;
    BScanHeader.OffSeg(zz) = OffSeg;
    BScanHeader.Shift(zz) = Shift;
    BScanHeader.Quality(zz) = Quality;
    BScanHeader.Spare(:, zz) = Spare;
 
    % Display the images
    if visu==1
        StartX_px=round(StartX/(header.SizeXSlo*header.ScaleXSlo)*header.SizeXSlo); %StartX in pixels
        StartY_px=round(StartY/(header.SizeYSlo*header.ScaleYSlo)*header.SizeYSlo); %StartY in pixels
        EndX_px=round(EndX/(header.SizeXSlo*header.ScaleXSlo)*header.SizeXSlo); %EndX in pixels
        EndY_px=round(EndY/(header.SizeYSlo*header.ScaleYSlo)*header.SizeYSlo); %EndY in pixels
        subplot(1,2,1); line([StartX_px EndX_px],[StartY_px EndY_px]);
    end
    
    %BScan reading
    if(numel(strfind(options, 'header')) == 0)
        status = fseek( fid, header.BScanHdrSize+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
        oct = fread( fid, header.SizeX*header.SizeZ, 'float32' );

        oct = reshape( oct, header.SizeX, header.SizeZ );
        oct = flipud( oct );
        oct = rot90( oct, 3);
        BScans(:,:,zz)=oct;
    end
 
    
    % Segmented data reading
    status = fseek( fid, BScanHeader.OffSeg(zz)+2048+(header.SizeXSlo*header.SizeYSlo)+(ii*(header.BScanHdrSize+header.SizeX*header.SizeZ*4)), -1 );
    seg = (fread( fid, NumSeg*header.SizeX, 'float' ))';
    for i = 1:NumSeg
        BScanHeader.(sprintf('Boundary_%d', i))(zz, :) = seg((i-1)*header.SizeX+1:i*header.SizeX);
    end
    
    
%    if visu==1
%         subplot(1,2,2);
%         imshow(BScans(:,:,zz));
%         hold on; plot(ilm(:),'b*');
%         drawnow
%     end

 
end
    
% % This is for ILM and RPE INVALID values interpolation using median of 3x3 neigborhood around
% % INVALID value.
%  
%     max_ILM = max(max(BScanHeader.ILM));
%     [N,M]=size(BScanHeader.ILM);
%     for j=1:N
%         for i=1:M
%             if BScanHeader.ILM(j,i) == max_ILM
%                 ILM_sec = imcrop(BScanHeader.ILM, [i-1 j-1 3 3]);
%                 BScanHeader.ILM(j,i) = median(find(ILM_sec ~= max_ILM));
%             end
%             BScanHeader.ILM(j,i) = header.SizeZ-BScanHeader.ILM(j,i);
%         end
%     end
%  
%     max_RPE = max(max(BScanHeader.RPE));
%     [N,M]=size(BScanHeader.RPE);
%     for j=1:N
%         for i=1:M
%             if BScanHeader.RPE(j,i) == max_RPE
%                 RPE_sec = imcrop(BScanHeader.RPE, [i-1 j-1 3 3]);
%                 BScanHeader.RPE(j,i) = median(find(RPE_sec ~= max_RPE));
%             end
%             BScanHeader.RPE(j,i) = header.SizeZ-BScanHeader.RPE(j,i);
%         end
%     end
%  
%     if visuseg == 1,
%         figure; 
%         mesh(BScanHeader.ILM(11:end-10,10:end-10));
%         hold on
%         mesh(BScanHeader.RPE(11:end-10,10:end-10));
%     end
    

%% Thickness grid read
if header.GridType ~= 0  % If the GridType is zero, there is no thickness information
    status = fseek(fid, header.GridOffset, -1);
    
    ThicknessGrid.Type = fread(fid, 1, 'int32');
    ThicknessGrid.Diameter = fread(fid, 3, 'double');
    ThicknessGrid.CenterPos = fread(fid, 2, 'double');
    ThicknessGrid.CentralThk = fread(fid, 1, 'float32');
    ThicknessGrid.MinCentralThk = fread(fid, 1, 'float32');
    ThicknessGrid.MaxCentralThk = fread(fid, 1, 'float32');
    ThicknessGrid.TotalVolume = fread(fid, 1, 'float32');
    for i = 1:9
        ThicknessGrid.Sectors(i).Thickness = fread(fid, 1, 'float32');
        ThicknessGrid.Sectors(i).Volume = fread(fid, 1, 'float32');
    end
else % If no ThicknessGrid is there, return ThicknessGrid as an empty structure
    ThicknessGrid = struct;
end

fclose( fid );

    %%
    if numel(strfind(options, 'metawrite')) ~= 0
        metafilename = [path(1:end-3) 'meta'];

        if((exist(metafilename) == 0) || numel(strfind(options, 'metawriteforce')) ~= 0)
            fidW = fopen(metafilename, 'w');

            metaWriteHelper(fidW, 'ExamTime', header.ExamTime);
            metaWriteHelper(fidW, 'ScanFocus', header.ScanFocus);
            metaWriteHelper(fidW, 'ScanPosition', header.ScanPosition, 'str');
            metaWriteHelper(fidW, 'ScanPattern', header.ScanPattern());
            metaWriteHelper(fidW, 'ID', header.ID, 'str');
            metaWriteHelper(fidW, 'ReferenceID', header.ReferenceID, 'str');
            metaWriteHelper(fidW, 'PID', header.PID()  );
            metaWriteHelper(fidW, 'PatientID', header.PatientID, 'str');
            metaWriteHelper(fidW, 'DOB', header.DOB);
            metaWriteHelper(fidW, 'VID', header.VID);
            metaWriteHelper(fidW, 'VisitID', header.VisitID, 'str');
            metaWriteHelper(fidW, 'VisitDate', header.VisitDate);
            metaWriteHelper(fidW, 'OctSize', [header.SizeX header.SizeZ header.NumBScans]);
            metaWriteHelper(fidW, 'OctScale', [header.ScaleX header.ScaleZ header.Distance]);
            metaWriteHelper(fidW, 'SloSize', [header.SizeXSlo header.SizeYSlo]);
            metaWriteHelper(fidW, 'SloScale', [header.ScaleXSlo header.ScaleYSlo]);
            metaWriteHelper(fidW, 'SloFieldSize', header.FieldSizeSlo);
            metaWriteHelper(fidW, 'BScanStartX', BScanHeader.StartX);
            metaWriteHelper(fidW, 'BScanStartY', BScanHeader.StartY);
            metaWriteHelper(fidW, 'BScanEndX', BScanHeader.EndX);
            metaWriteHelper(fidW, 'BScanEndY', BScanHeader.EndY);
            
            fclose( fidW );
        end
    end
end

    function metaWriteHelper(fidW, tag, data, mode)
        if nargin < 4
            mode = 'num';
        end
        if strcmp(mode, 'num')
            str = sprintf('%g ', data);
        else
            str = sprintf('%s ', data);
        end       
        outputstring = [tag ' ' str];
        outputstring = deblank(outputstring);
        fprintf(fidW, '%s\n', outputstring);
    end