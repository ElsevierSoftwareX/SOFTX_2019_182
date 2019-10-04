function varargout = AMPAO(varargin)
% AMPAO MATLAB code for AMPAO.fig
%      AMPAO, by itself, creates a new AMPAO or raises the existing
%      singleton*.
%
%      H = AMPAO returns the handle to a new AMPAO or the handle to
%      the existing singleton*.
%
%      AMPAO('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in AMPAO.M with the given input arguments.
%
%      AMPAO('Property','Value',...) creates a new AMPAO or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before AMPAO_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to AMPAO_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help AMPAO

% Last Modified by GUIDE v2.5 15-Apr-2019 03:52:14

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @AMPAO_OpeningFcn, ...
    'gui_OutputFcn',  @AMPAO_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before AMPAO is made visible.
function AMPAO_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to AMPAO (see VARARGIN)

% Choose default command line output for AMPAO
handles.output = hObject;
handles.LoadFlag = 0;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes AMPAO wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = AMPAO_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in Load_Data.
function Load_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
clc
% close all
[FileName,PathName] = uigetfile({'*.xml;*.vol;*.mat'},'Select xml, vol, or mat file');
handles.FileName=FileName;
guidata(hObject,handles);
handles.PathName=PathName;
guidata(hObject,handles);
handles.data = renameAMPAO(FileName,PathName,handles);
handles.LoadFlag = 1;
guidata(hObject,handles);
% disp(handles)


function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Convert_to.
function Convert_to_Callback(hObject, eventdata, handles)
% hObject    handle to Convert_to (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on selection change in conversion.
function conversion_Callback(hObject, eventdata, handles)
% hObject    handle to conversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns conversion contents as cell array
%        contents{get(hObject,'Value')} returns selected item from conversion
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    FileName=handles.FileName;
    PathName=handles.PathName;
    conversion(FileName,PathName,handles);
end
% --- Executes during object creation, after setting all properties.
function conversion_CreateFcn(hObject, eventdata, handles)
% hObject    handle to conversion (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Display_Data.
function Display_Data_Callback(hObject, eventdata, handles)
% hObject    handle to Display_Data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% cla(handles.axes1,'reset');
% cla(handles.axes2,'reset');

% PathName=handles.PathName;
% FileName=handles.FileName;
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    displayAMPAO(handles)
end

% --- Executes on button press in Alignment.
function Alignment_Callback(hObject, eventdata, handles)
% hObject    handle to Alignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    output = curvCorrectionAMPAO(handles);
    handles.data.PreProcess = output;
    guidata(hObject,handles);
end

% --- Executes on button press in Denoising.
function Denoising_Callback(hObject, eventdata, handles)
% hObject    handle to Denoising (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% scan_num = 1;
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    output = DenoisingAMPAO(handles);
    handles.data.PreProcess = output;
    guidata(hObject,handles);
end

%
% function ScanNum_denoise_Callback(hObject, eventdata, handles)
% % hObject    handle to ScanNum_denoise (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
%
% % Hints: get(hObject,'String') returns contents of ScanNum_denoise as text
% %        str2double(get(hObject,'String')) returns contents of ScanNum_denoise as a double
% user_entry = str2double(get(hObject,'string'));
% if isnan(user_entry)
% errordlg('You must enter a numeric value','Bad Input','modal')
% end
% handles.ScanNum_denoise = user_entry;
% print(user_entry)

% --- Executes during object creation, after setting all properties.
function ScanNum_denoise_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ScanNum_denoise (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Display_PreProcess.
function Display_PreProcess_Callback(hObject, eventdata, handles)
% hObject    handle to Display_PreProcess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    displayPreProcessAMPAO(handles)
    
end

% --- Executes on button press in SegMethod1.
function SegMethod1_Callback(hObject, eventdata, handles)
% hObject    handle to SegMethod1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    SegmentationMethod1AMPAO(handles);
    % handles.data.PreProcess = output;
    % guidata(hObject,handles);
end

% --- Executes on button press in SegMethod2.
function SegMethod2_Callback(hObject, eventdata, handles)
% hObject    handle to SegMethod2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    uiwait(msgbox('The operation may take a while, Please wait...','Warning','modal'));
    
    SegmentationMethod2AMPAO(handles);
    guidata(hObject,handles);
end

% --- Executes on button press in thicknessMap.
function thicknessMap_Callback(hObject, eventdata, handles)
% hObject    handle to thicknessMap (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Construct a questdlg with three options

if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    choice = questdlg('Select Thickness map Calculation Method', ...
        'Select Thickness map Calculation Method', ...
        'Saved Segmentation result','Segmentated boundaries by device','Cancel, need segmentation Block','Cancel, need segmentation Block');
    % Handle response
    switch choice
        case 'Saved Segmentation result'
            %         disp([choice ' coming right up.'])
            [FileName,PathName] = uigetfile('*.mat','Select the Boundary Points');
            load(strcat(PathName,FileName));
            % thicknessMapAMPAO(bd_pts,header,[0.5 1.5 2.5]*1000,1);
            thicknessMapAMPAO(bd_pts,header,[0.5 1.5 3]*1000,0,handles);
            handles.data.bd_pts_Segmented = bd_pts;
            handles.data.headerThickness_Segmented = header;
            guidata(hObject,handles);
        case 'Segmentated boundaries by device'
            %         disp([choice ' coming right up.'])
            thicknessMapAMPAO(handles.data.bd_pts,handles.data.headerThickness,[0.5 1.5 3]*1000,0,handles);
        case 'Cancel, need segmentation Block'
            %         disp('I''ll bring you your check.')
    end
end



% --- Executes on button press in DisplayETDRS.
function DisplayETDRS_Callback(hObject, eventdata, handles)
% hObject    handle to DisplayETDRS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
    
    
else
    choice = questdlg('Select Thickness map Calculation Method', ...
        'Select Thickness map Calculation Method', ...
        'Saved Segmentation result','Segmentated boundaries by device','Cancel, need segmentation Block','Cancel, need segmentation Block');
    % Handle response
    switch choice
        case 'Saved Segmentation result'
            warndlg('Be sure to run Thickness Map in (Saved Segmentation result) mode  before this stage','!! Warning !!')
            
            %         disp([choice ' coming right up.'])
            displayThicknessMapAMPAO(handles.data.bd_pts_Segmented,handles.data.headerThickness_Segmented,[0.5 1.5 3]*1000,1,handles);
        case 'Segmentated boundaries by device'
            %         disp([choice ' coming right up.'])
            displayThicknessMapAMPAO(handles.data.bd_pts,handles.data.headerThickness,[0.5 1.5 3]*1000,1,handles);
        case 'Cancel, need segmentation Block'
            %         disp('I''ll bring you your check.')
    end
end
% thicknessMapAMPAO2(handles.data.bd_pts,handles.data.headerThickness,[0.5 1.5 3]*1000,1,handles);


% --- Executes on button press in ONHLoc.
function ONHLoc_Callback(hObject, eventdata, handles)
% hObject    handle to ONHLoc (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    if ~strcmp(handles.data.format,'vol')
        errordlg('Method II works only on VOL data','File Error');
        return
    end
    
    SLO=handles.data.slo;
    
    if handles.data.ScanPosition(:,1:2)=='OS'
        SLO = fliplr(SLO);
        
    end
    SLO=uint8(SLO);
    
    [center00,radius0,x0,y0,x00,y00]=fovea_onh_readimage(SLO);
    vab=SLO;
    a1=[x0 x00];
    b1=[y0 y0];
    a2=[x0 x00];
    b2=[y0 y00+39];
    
    diff = ( - atan((b2(2)-b2(1))/(a2(2)-a2(1)))) * 180/pi;
    %     uisave({'FoveaLoc','ONHLoc'})
    if handles.data.ScanPosition(:,1:2)=='OS'
        diff=-diff;
    end
    handles.temporarysloRotation = diff;
    guidata(hObject,handles);
end




% --- Executes on button press in ONHSeg.
function ONHSeg_Callback(hObject, eventdata, handles)
% hObject    handle to ONHSeg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.LoadFlag==0
    warndlg('You should load your data in Convert Block','!! Warning !!')
else
    
    choice = questdlg('If loaded data is NOT from ONH, you will encounter an ERROR; Do you want to continue?', ...
        'ONH data WARNING','Yes', 'No','Cancel');
    % Handle response
    switch choice
        case 'Yes'
            uiwait(msgbox('The operation may take a while, Please wait...','Warning','modal'));
            
            output = ONHSegAMPAO(handles);
            handles.data.ONHSeg = output;
            guidata(hObject,handles);
        case 'No'
            
        case 'Cancel'
    end
    
    
    
end

% --- Executes on button press in SLOalignment.
function SLOalignment_Callback(hObject, eventdata, handles)
% hObject    handle to SLOalignment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Construct a questdlg with three options
% choice = questdlg('Are you sure to rotate thickness maps?', ...
% 	'Yes', ...
% 	'No');
% % Handle response
% switch choice
%     case 'Yes'
%         handles.data.sloRotation = handles.temporarysloRotation;
%     case 'No'
% end
