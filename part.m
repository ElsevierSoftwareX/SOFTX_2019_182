
function info=part(address)

docNode = xmlread(address);

%PatientID
data = docNode.getElementsByTagName('PatientID');
info.PatientID = char(data.item(0).getFirstChild.getData);
%LastName
data = docNode.getElementsByTagName('LastName');
info.LastName = char(data.item(0).getFirstChild.getData);
%FirstNames
data = docNode.getElementsByTagName('FirstNames');
info.firstname = char(data.item(0).getFirstChild.getData);
%Year
data = docNode.getElementsByTagName('Year');
info.Year = char(data.item(0).getFirstChild.getData);
%Sex
data = docNode.getElementsByTagName('Sex');
info.Sex = char(data.item(0).getFirstChild.getData);
%Laterality
data = docNode.getElementsByTagName('Laterality');
info.Laterality = char(data.item(0).getFirstChild.getData);
%Diameter1
data = docNode.getElementsByTagName('Diameter1');
info.Diameter1 = char(data.item(0).getFirstChild.getData);
%Diameter2
data = docNode.getElementsByTagName('Diameter2');
info.Diameter2 = char(data.item(0).getFirstChild.getData);
%Diameter3
data = docNode.getElementsByTagName('Diameter3');
info.Diameter3 = char(data.item(0).getFirstChild.getData);
%X
data = docNode.getElementsByTagName('X');
info.Coord.X = char(data.item(0).getFirstChild.getData);
%Y
data = docNode.getElementsByTagName('Y');
 info.Coord.Y = char(data.item(0).getFirstChild.getData);
% % CenterPos
% data = docNode.getElementsByTagName('CenterPos');
% info.CenterPos = char(data.item(0).getFirstChild.getData);
%CentralThickness
data = docNode.getElementsByTagName('CentralThickness');
info.CenterPos.CentralThickness = char(data.item(0).getFirstChild.getData);
%MinCentralThickness
data = docNode.getElementsByTagName('MinCentralThickness');
info.CenterPos.MinCentralThickness = char(data.item(0).getFirstChild.getData);
%MaxCentralThickness
data = docNode.getElementsByTagName('MaxCentralThickness');
info.CenterPos.MaxCentralThickness = char(data.item(0).getFirstChild.getData);
%TotalVolume
data = docNode.getElementsByTagName('TotalVolume');
info.CenterPos.TotalVolume = char(data.item(0).getFirstChild.getData);

% C0
data = docNode.getElementsByTagName('Name');
info.Name.C0 = char(data.item(3).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.C0 = char(data.item(0).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.C0 = char(data.item(0).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.C0 = char(data.item(0).getFirstChild.getData);
% N1
data = docNode.getElementsByTagName('Name');
info.Name.N1 = char(data.item(4).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.N1 = char(data.item(1).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.N1 = char(data.item(1).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.N1 = char(data.item(1).getFirstChild.getData);
% N2
data = docNode.getElementsByTagName('Name');
info.Name.N2 = char(data.item(5).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.N2 = char(data.item(2).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.N2 = char(data.item(2).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.N2 = char(data.item(2).getFirstChild.getData);
% S1
data = docNode.getElementsByTagName('Name');
info.Name.S1 = char(data.item(6).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.S1 = char(data.item(3).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.S1 = char(data.item(3).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.S1 = char(data.item(3).getFirstChild.getData);
% S2
data = docNode.getElementsByTagName('Name');
info.Name.S2 = char(data.item(7).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.S2 = char(data.item(4).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.S2 = char(data.item(4).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.S2 = char(data.item(4).getFirstChild.getData);
% T1
data = docNode.getElementsByTagName('Name');
info.Name.T1 = char(data.item(8).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.T1 = char(data.item(5).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.T1 = char(data.item(5).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.T1 = char(data.item(5).getFirstChild.getData);
% T2
data = docNode.getElementsByTagName('Name');
info.Name.T2 = char(data.item(9).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.T2 = char(data.item(6).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.T2 = char(data.item(6).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.T2 = char(data.item(6).getFirstChild.getData);
% I1
data = docNode.getElementsByTagName('Name');
info.Name.I1 = char(data.item(10).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.I1 = char(data.item(7).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.I1 = char(data.item(7).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.I1 = char(data.item(7).getFirstChild.getData);
% I2
data = docNode.getElementsByTagName('Name');
info.Name.I2 = char(data.item(11).getFirstChild.getData);
%AvgThickness
data = docNode.getElementsByTagName('AvgThickness');
info.AvgThickness.I2 = char(data.item(8).getFirstChild.getData);
%Volume
data = docNode.getElementsByTagName('Volume');
info.Volume.I2 = char(data.item(8).getFirstChild.getData);
%ValidPixelPercentage
data = docNode.getElementsByTagName('ValidPixelPercentage');
info.ValidPixelPercentage.I2 = char(data.item(8).getFirstChild.getData);
% NumImages
data = docNode.getElementsByTagName('NumImages');
info.NumImages = char(data.item(0).getFirstChild.getData);
% % for all ID
%Hour
data = docNode.getElementsByTagName('Hour');
info.Hour = char(data.item(0).getFirstChild.getData);
%Minute
data = docNode.getElementsByTagName('Minute');
info.Minute = char(data.item(0).getFirstChild.getData);
%UTCBias
data = docNode.getElementsByTagName('UTCBias');
info.UTCBias = char(data.item(0).getFirstChild.getData);
%ScaleX
data = docNode.getElementsByTagName('ScaleX');
info.ScaleX = char(data.item(0).getFirstChild.getData);
%ScaleY
data = docNode.getElementsByTagName('ScaleY');
info.ScaleY = char(data.item(0).getFirstChild.getData);

for i = 0 : 19
docNode = xmlread(address);
data = docNode.getElementsByTagName('ID');
ID = char(data.item(i+3).getFirstChild.getData);
info.ID(i+1).ID = ID;
% %Second
data = docNode.getElementsByTagName('Second');
Second = char(data.item(i).getFirstChild.getData);
info.ID(i+1).Second = Second;
%NumAve
data = docNode.getElementsByTagName('NumAve');
NumAve = char(data.item(i).getFirstChild.getData);
info.ID(i+1).NumAve = NumAve;
% %Width
data = docNode.getElementsByTagName('Width');
Width = char(data.item(i).getFirstChild.getData);
info.ID(i+1).Width = Width;
%Height
data = docNode.getElementsByTagName('Height');
Height = char(data.item(i).getFirstChild.getData);
info.ID(i+1).Height = Height;
end

%for ID = 0
%Angle
data = docNode.getElementsByTagName('Angle');
info.ID(1).Angle = char(data.item(0).getFirstChild.getData);
%Focus
data = docNode.getElementsByTagName('Focus');
info.ID(1).Focus = char(data.item(0).getFirstChild.getData);
%SensorGain
data = docNode.getElementsByTagName('SensorGain');
info.ID(1).SensorGain = char(data.item(0).getFirstChild.getData);

%%%%%%%%%%%%%%% images %%%%%%%%%%%%%%%%%%%%%%%%%%%
%SLO image
data = docNode.getElementsByTagName('ExamURL');
im_slo = char(data.item(0).getFirstChild.getData);
loc = strfind(im_slo,'\') ; 
im_slo(1:loc(end)) =[];
info.slo = imread(im_slo)

%OCT image
for k = 1 : 19
data = docNode.getElementsByTagName('ExamURL');
im_oct = char(data.item(k).getFirstChild.getData);
loc = strfind(im_oct,'\') ; 
im_oct(1:loc(end)) =[];
info.oct = imread(im_oct);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PatientID = info.PatientID;
LastName = info.LastName;
firstname = info.firstname;
Year = info.Year;
Sex = info.Sex;
Laterality = info.Laterality;
Diameter1 = info.Diameter1;
Diameter2 = info.Diameter2;
Diameter3 = info.Diameter3;
X = info.Coord.X;
Y = info.Coord.Y;
CentralThickness = info.CenterPos.CentralThickness;
MinCentralThickness = info.CenterPos.MinCentralThickness;
MaxCentralThickness = info.CenterPos.MaxCentralThickness;
TotalVolume = info.CenterPos.TotalVolume;
%%C0
Name.C0 = info.Name.C0;
AvgThickness.C0 = info.AvgThickness.C0;
Volume.C0 = info.Volume.C0;
ValidPixelPercentage.C0 = info.ValidPixelPercentage.C0;
%%N1
Name.N1 = info.Name.N1;
AvgThickness.N1 = info.AvgThickness.N1;
Volume.N1 = info.Volume.N1;
ValidPixelPercentage.N1 = info.ValidPixelPercentage.N1;
%%N2
Name.N2 = info.Name.N2;
AvgThickness.N2 = info.AvgThickness.N2;
Volume.N2 = info.Volume.N2;
ValidPixelPercentage.N2 = info.ValidPixelPercentage.N2;
%%S1
Name.S1 = info.Name.S1;
AvgThickness.S1 = info.AvgThickness.S1;
Volume.S1 = info.Volume.S1;
ValidPixelPercentage.S1 = info.ValidPixelPercentage.S1;
%%S2
Name.S2 = info.Name.S2;
AvgThickness.S2 = info.AvgThickness.S2;
Volume.S2 = info.Volume.S2;
ValidPixelPercentage.S2 = info.ValidPixelPercentage.S2;
%%T1
Name.T1 = info.Name.T1;
AvgThickness.T1 = info.AvgThickness.T1;
Volume.T1 = info.Volume.T1;
ValidPixelPercentage.T1 = info.ValidPixelPercentage.T1;
%%T2
Name.T2 = info.Name.T2;
AvgThickness.T2 = info.AvgThickness.T2;
Volume.T2 = info.Volume.T2;
ValidPixelPercentage.T2 = info.ValidPixelPercentage.T2;
%%I1
Name.I1 = info.Name.I1;
AvgThickness.I1 = info.AvgThickness.I1;
Volume.I1 = info.Volume.I1;
ValidPixelPercentage.I1 = info.ValidPixelPercentage.I1;
%%I2
Name.I2 = info.Name.I2;
AvgThickness.I2 = info.AvgThickness.I2;
Volume.I2 = info.Volume.I2;
ValidPixelPercentage.I2 = info.ValidPixelPercentage.I2;
%%%
NumImages = info.NumImages;
Hour = info.Hour;
Minute = info.Minute;
UTCBias = info.UTCBias;
Scalex = info.ScaleX;
ScaleY = info.ScaleY;
ID = info.ID;
SLO = info.slo;
OCT = info.oct;

save('info.mat','-struct','info');