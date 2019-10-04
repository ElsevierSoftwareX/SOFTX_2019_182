% function b = writeXml(A,param,filename)
%
% Input A is the array to be written
% Input param is a structure containing the image meta-data
% Input filename is the name of the file to be written
% John Bogovic
% 08/04/2009

% Writes matlab arrays to MIPAV xml-raw pair

function b = writeXml(A,param,filename)

b = 0;

xmlfile = strcat(filename,'.xml');
rawfile = strcat(filename,'.raw');

dim = size(A);
res = ones(size(dim));
slcspc = 1;
slcthk = 0;
numdims = size(dim,2);


orient = 'Unknown';
if(numdims>=3)
    axisorient = {'Unknown','Unknown','Unknown'};
elseif(numdims==2)
    axisorient = {'Unknown','Unknown'};
end

%%%%%%%%%%%%%%%%%%%
% % READ HEADER % %
%%%%%%%%%%%%%%%%%%%

% each row is a data type
% first col is mipav name, second is matlab name
typemap = {'Double','double'; 'Float','float'; 'Integer','int';
    'Short','short';'Byte', 'schar'; 'Long', 'long';
    'Unsigned Integer', 'uint32';'Unsigned Short', 'ushort'; 'Unsigned Byte', 'uint8'};

if(isfield(param,'res'))
    res = param.res;
end
if(isfield(param,'type'))
    type = param.type;
    typeind = find(strcmp({typemap{:,1}},type));
    if(typeind>0)
        mattype = typemap{typeind,2};
    else
       disp('DEFAULTING TO FLOAT TYPE FOR OUTPUT');
       type = 'Float';
        mattype = 'float';
    end
    
else
    mattype = class(A);
    typeind = find(strcmp({typemap{:,2}},mattype));
    if(typeind>0)
        type = typemap{typeind,1};
    else
       disp('DEFAULTING TO FLOAT TYPE FOR OUTPUT');
       type = 'Float';
       mattype = 'float';
    end
end
if(isfield(param,'sliceSpacing'))
    slcspc = param.sliceSpacing;
end
if(isfield(param,'sliceThickness'))
    slcthk = param.sliceThickness;
end
if(isfield(param,'orientation'))
    orient = param.orientation;
end
if(isfield(param,'axisOrientation'))
    axisorient = param.axisOrientation;
end


%%%%%%%%%%%%%%%%%%%%%%
% % write raw file
%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(rawfile,'w');
braw = fwrite(fid,A,mattype);
fclose(fid);

    
%%%%%%%%%%%%%%%%%%%%%%
% % write xml file
%%%%%%%%%%%%%%%%%%%%%%

fid = fopen(xmlfile,'w');

fprintf(fid,'<?xml version="1.0" encoding="UTF-8"?>\n');
fprintf(fid,'<!-- MIPAV header file -->\n');
fprintf(fid,'<image xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" nDimensions="%c">\n',num2str(numdims));
fprintf(fid,'\t<Dataset-attributes>\n');

fprintf(fid,'\t\t<Image-offset>0</Image-offset>\n');
fprintf(fid,'\t\t<Data-type>%s</Data-type>\n',type);
fprintf(fid,'\t\t<Endianess>Little</Endianess>\n');
for i=1:length(dim)
    fprintf(fid,'\t\t<Extents>%d</Extents>\n',dim(i));
end

fprintf(fid,'\t\t<Resolutions>\n');
for i=1:length(res)
    fprintf(fid,'\t\t\t<Resolution>%g</Resolution>\n',res(i));
end
fprintf(fid,'\t\t</Resolutions>\n');

if(numdims>=3)
    fprintf(fid,'\t\t<Slice-spacing>%g</Slice-spacing>\n',slcspc);
    fprintf(fid,'\t\t<Slice-thickness>%g</Slice-thickness>\n',slcthk);
end


for i=1:length(res)
    fprintf(fid,'\t\t<Units>Millimeters</Units>\n');
end
fprintf(fid,'\t\t<Compression>none</Compression>\n');
fprintf(fid,'\t\t<Orientation>%s</Orientation>\n',orient);

for i=1:numdims
    if(i>length(axisorient))
        
    else
        fprintf(fid,'\t\t<Subject-axis-orientation>%s</Subject-axis-orientation>\n',axisorient{i});
    end
    
end

for i=1:length(res)
    fprintf(fid,'\t\t<Origin>0.0</Origin>\n');
end
fprintf(fid,'\t\t<Modality>Unknown Modality</Modality>\n');

fprintf(fid,'\t</Dataset-attributes>\n');
fprintf(fid,'</image>\n');

fclose(fid);
end