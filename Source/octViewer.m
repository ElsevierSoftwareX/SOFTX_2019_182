function octViewer(octvol,segpts)
% octViewer.m
% 
%   GUI to display OCT data and segmentation results.
%
%   To display an OCT volume:
%       - Run octViewer without any arguments and go to File->Load Volume
%         and choose the desired format. Compatible formats include the
%         .vol format from Spectralis scanners, the .img format from
%         Cirrus (filename format should be XX_cube_z.img or
%         XX_cube_raw.img), or a saved MATLAB .mat file containing the
%         variable 'img_vol' representing the data as a 3 dimensional
%         array.
%       - Run octViewer with the first argument being a 3 dimensional array
%         representing the data (ex. a 496x1024x49 array containing 496
%         pixels per A-scan, 1024 A-scans per B-scan, and 49 B-scans per
%         volume).
%   To display a segmentation result:
%       - Run octViewer without any arguments and go to Segmentation->Load
%         File and choose the desired format. Compatible formats include a
%         .txt file from the OCT Toolbox manual segmentation tool, or a
%         .mat file stored after running the OCTLayerSegmentation function
%         in MATLAB.
%       - Run octViewer with the first argument being a 3 dimensional array
%         containing the position of each boundary at every A-scan (ex. a
%         1024x49x9 array containing the boundary position along each of
%         1024 A-scans per B-scan, 49 B-scans, and 9 boundaries).

% TODO:
%   Features:
%       * layer segmentation tools (including choroid):
%           * manual layer segmentation (load and save function, tag with initials of rater)
%           * automatic layer segmentation (multiple methods)
%       * cyst tools:
%           * manual segmentation
%           * automatic segmentation
%           * fundus density
%       * fundus image (en-face and SLO)
%       * Image processing tools:
%           * intensity normalization
%           * filtering
%       * display gradient image
%       * flatten
%       * fix aspect ratio
%       * keyboard shortcuts - change image, toggle segmentation
%
%   Bug fixes:
%       * prevent Pan tool from moving outside image
%       * setting "FontUnits 'normalized'" to all uicontrols breaks linux
%
% Andrew Lang (alang9@jhu.edu)
% Last updated: 12/4/2013

%% Input handling

startx = 1;
starty = 1;
if nargin == 2
    % Second argument is boundary segmentation
    if ~isempty(segpts)
        szv = size(octvol); % ex. [496 1024 49]
%         if size(segpts,1) ~= szv(2)
%             segpts = permute(segpts,[2 1 3 4]);
%         end
%         sz = size(segpts); % ex. [1024 49 9]
%         if (sz(1) ~= szv(2)) || (sz(2) ~= szv(3))
%             error('Invalid size for segmentation points')
%         end
        
        % Segmentation gets centered on the data
        startx = round(size(octvol,2)/2) - floor((size(segpts,1)-1)/2);
        starty = round(size(octvol,3)/2) - floor((size(segpts,2)-1)/2);
    end
end
if nargin < 2
    segpts = [];
end
if nargin < 1
    octvol = [];
end

%%  Initialization tasks

figure_name = 'OCT Viewer v0.1';

% Build figure
hf = figure('Visible','on',...
            'Name',figure_name,...
            'Resize','on',...
            'MenuBar','none',...
            'NumberTitle','off',...
            'Toolbar','figure',...
            'Units','normalized',... % normalized - beware - different size on different monitors
            'Position',[0.1 0.1 0.8 0.8]); % [left bottom width height]
    
control_pts = [];
voxel_size = [];
% startx = 1;
% starty = 1;
        
% Handles to plot objects
h_img = [];
h_pts = zeros(11,1);

%%  Construct the components
ha = axes('Parent',hf,...
          'Position',[0.05 0.2 0.9 0.75],...
          'XTick',[],...
          'YTick',[],...
          'box','on',...
          'DrawMode','fast');
      
% Context menu for right clicking on image
% uicontextmenu

% Text below image to say 'image x of y'
uicontrol('Style','text',...
          'Units','normalized',...
          'Position',[.87 .175 .03 .02], ...
          'FontSize',8, ...
          'BackgroundColor',get(hf,'Color'),...
          'String','image');      
      
% Edit box to control the current bscan index
hedit_sl = uicontrol('Style','edit', ...
                      'Units','normalized',...
                      'Position',[.9 .175 .02 .02], ...
                      'FontSize',8, ...
                      'BackgroundColor','white', ...
                      'String','-');

uicontrol('Style','text',...
          'Units','normalized',...
          'Position',[.918 .175 .01 .02], ...
          'FontSize',8, ...
          'BackgroundColor',get(hf,'Color'),...
          'String','of');

% Text box displaying total number of images in the volume
htxt_sl = uicontrol('Style','text', ...
                    'Units','normalized',...
                    'Position',[.93 .175 .02 .02], ...
                    'FontSize',8, ...
                    'BackgroundColor','white', ...
                    'String','-');
                
% Check box to hide the segmentation boundaries
hcb_seg = uicontrol('Style','checkbox',...
                    'String','Hide segmentation',...
                    'Units','normalized',...
                    'Position',[0.05 0.15 0.09 0.02],...
                    'FontSize',8,...
                    'BackgroundColor','white',...
                    'Enable','off');
                      
% Text box for debug output
htxt_db = uicontrol('Style','text',...
                    'Units','normalized',...
                    'Position',[.4 .175 .2 .02], ...
                    'FontSize',8, ...
                    'BackgroundColor',get(hf,'Color'),...
                    'String','');

% Toolbar - only want Zoom In, Zoom Out, Pan
htb = findall(hf,'type','UIToolbar');
hb = findobj(allchild(htb),'-not','ToolTipString','Zoom In',...
                           '-not','ToolTipString','Zoom Out',...
                           '-not','ToolTipString','Pan');
set(hb,'Visible','off')

% Menu bar
hmf = uimenu(hf,'Label','File');
hmf_load = uimenu(hmf,'Label','Load volume');
uimenu(hmf_load,'Label','Spectralis (.vol) file',...
                'Callback',@load_vol,...
                'Tag','vol')
uimenu(hmf_load,'Label','Cirrus (.img) file',...
                'Callback',@load_vol,...
                'Tag','img')
uimenu(hmf_load,'Label','.mat file',...
                'Callback',@load_vol,...
                'Tag','mat')

hms = uimenu(hf,'Label','Segmentation','Enable','off');
hms_load = uimenu(hms,'Label','Load file');
uimenu(hms_load,'Label','OCT Toolbox file',...
                'Callback',@load_seg,...
                'Tag','otb')
uimenu(hms_load,'Label','.mat file',...
                'Callback',@load_seg,...
                'Tag','mat')

%%  Initialization tasks

if ~isempty(octvol)
    octvol = double(octvol);
    octvol = (octvol - min(octvol(:)))/(max(octvol(:))-min(octvol(:)));
    
    curimg = round(size(octvol,3)/2);
    
    
    
    set(htxt_sl,'String',num2str(size(octvol,3)));
    
    set(hms,'Enable','on');
    
    updateimage();
    
    if ~isempty(segpts)
        set(hcb_seg,'Enable','on')

        control_pts = cell(size(segpts,2),11);
        for i = 1:size(segpts,2)
            for j = 1:size(segpts,3)
                control_pts{i,j} = [1:size(segpts,1); segpts(:,i,j)']';
            end
            control_pts{i,10} = [];
            control_pts{i,11} = [];
        end
        
        updateseg();
    end
end


%%  Callbacks for MYGUI

set(hf,'WindowScrollWheelFcn',@scrollfcn);
set(hedit_sl,'Callback',@change_image);
set(hcb_seg,'Callback',@hide_seg);

% uiwait(hf)



%%  Utility functions for MY

    % ------------------------------------------------------------------- %
    % Callback for mouse scrollwheel to change the current image
    function scrollfcn(hfig,events)
    % Callback for mouse scrollwheel
    %   hfig: handle to figure
    %   events: event structure, containing members
    %       VerticalScrollCount  - +1 for scroll down, -1 for scroll up
    %       VerticalScrollAmount - lines to scroll per count (unused here)
    % Called when scroll wheel is adjusted in figure hfig.
        cp = get(hfig,'CurrentPoint');  % Where mouse is, in figure units
        obj_pos = get(ha,'Position');   % Axes bounds

        % Is mouse pointer within axes?
        if cp(1) > obj_pos(1) && cp(1) < obj_pos(1) + obj_pos(3) && ...
           cp(2) > obj_pos(2) && cp(2) < obj_pos(2) + obj_pos(4)

            sc = events.VerticalScrollCount;
            sa = events.VerticalScrollAmount;
            
            if ~isempty(octvol)
                % Change current B-scan image
                curimg = curimg - sc;
                
                updateimage();
                if ~isempty(control_pts)
                    updateseg();
                end
            end
        else                      
            % Out of axis
        end    
    end
    
    % ------------------------------------------------------------------- %
    % Update the image (display a new slice)
    function updateimage
        if ~isempty(octvol)
            if curimg < 1
                curimg = 1;
            elseif curimg > size(octvol,3)
                curimg = size(octvol,3);
            end
            
            % Show image
            if isempty(h_img)
                h_img = imagesc(repmat(octvol(:,:,curimg).^0.25,[1 1 3]),'Parent',ha,'HitTest','off');
%                 colormap gray
                axis(ha,'off')
                axis(ha,'tight')
                hold(ha,'on')
            else
                set(h_img,'CData',repmat(octvol(:,:,curimg).^0.25,[1 1 3]));
            end
            
            % Display current bscan index
            set(hedit_sl,'String',num2str(curimg));
        end        
    end

    % ------------------------------------------------------------------- %
    % Update all boundary points on an image
    function updateseg
        hideseg = (get(hcb_seg,'Value') == get(hcb_seg,'Max'));
        
        % Remove boundaries from plot
        delete(h_pts(h_pts>0));
        h_pts(h_pts>0) = 0;
        
        % Add new points
        if ~hideseg && curimg >= starty && curimg+starty-1 <= size(octvol,3)
            for i = 1:size(control_pts,2)
                ctrl_pts = control_pts{curimg-starty+1,i};
                if isempty(ctrl_pts)
                    continue
                end
                
                % Draw lines
                if size(ctrl_pts,1) > 3
                    [~,inds] = sort(ctrl_pts(:,1));
                    ctrl_pts = ctrl_pts(inds,:);
                    if size(ctrl_pts,1) ~= size(octvol,2)
                        pts = interpolateCtrlPts(ctrl_pts);
                    else
                        pts = ctrl_pts;
                    end
                    
                    h_pts(i) = line(pts(:,1)+startx-1,pts(:,2),'Color','r','Parent',ha,'HitTest','off');
                end                
            end
        end
    end

    % ------------------------------------------------------------------- %
    function change_image(hObject,eventdata)
        user_entry = str2double(get(hObject,'string'));
        if isnan(user_entry)
            set(hObject,'String',curimg)
        else
            curimg = user_entry;            
            updateimage();
            updateseg();
        end
    end

    % ------------------------------------------------------------------- %
    % Click "Hide segmentation" button
    function hide_seg(hObject, eventdata)                
        % Update display
        updateseg();
    end
    
    % ------------------------------------------------------------------- %
    % Load OCT volume
    function load_vol(hObject, eventdata)
        set(hf,'Name',figure_name);
        
        delete(h_pts(h_pts>0));
        h_pts(h_pts>0) = 0;
        control_pts = [];
        
        tag = get(hObject,'Tag');
        pathname = '.';
        switch tag
            case 'vol'
                fileext = '*.vol';
                if ispref('octviewer','vol_path')
                    pathname = getpref('octviewer','vol_path');
                else
                    addpref('octviewer','vol_path','.');
                end
            case 'img'
                fileext = '*.img';
                if ispref('octviewer','img_path')
                    pathname = getpref('octviewer','img_path');
                else
                    addpref('octviewer','img_path','.');
                end
            case 'mat'
                fileext = '*.mat';
                if ispref('octviewer','mat_vol_path')
                    pathname = getpref('octviewer','mat_vol_path');
                else
                    addpref('octviewer','mat_vol_path','.');
                end
        end
        [filename, pathname] = uigetfile(fileext,'Load OCT data',pathname);
        
        if isequal(filename,0)
           return
        end
        filefull = fullfile(pathname, filename);        
        
        set(htxt_db,'String','Loading...','ForegroundColor','r','FontWeight','bold')
        drawnow
        
        switch tag
            case 'vol'
                setpref('octviewer','vol_path',pathname)
                
                [header, BScanHeader, slo, octvol] = openVolFast(filefull,'nodisp');
                octvol = double(octvol).^0.25;
                octvol(octvol>1) = 0;
                
                voxel_size = [header.ScaleZ header.ScaleX header.Distance];
            case 'img'
                setpref('octviewer','img_path',pathname)
                
                [octvol,vol_info] = octCirrusReader(filefull);
                voxel_size = vol_info.vol_res;
            case 'mat'
                setpref('octviewer','mat_vol_path',pathname)
                
                S = load(filefull);
                octvol = double(S.img_vol).^0.25;
                voxel_size = [S.header.ScaleZ S.header.ScaleX S.header.Distance];
        end
        
        % Return focus to figure
        figure(hf);
        
        % Set image data
        octvol = double(octvol);
        octvol = (octvol - min(octvol(:)))/(max(octvol(:))-min(octvol(:)));
        
        % Set current image index
        curimg = round(size(octvol,3)/2);
        set(htxt_sl,'String',num2str(size(octvol,3)));
        
        if ~isempty(h_img)
            delete(h_img)
            h_img = [];
        end
        
        updateimage();
        
        set(hcb_seg,'Enable','on')
        
        set(hf,'Name',[figure_name ' - ' filename]);
        set(hms,'Enable','on');
        set(htxt_db,'String','')
    end
 
    % ------------------------------------------------------------------- %
    % Load segmentation file
    function load_seg(hObject, eventdata)
        tag = get(hObject,'Tag');
        pathname = '.';
        switch tag
            case 'otb'
                fileext = '*.txt';
                if ispref('octviewer','seg_path_otb')
                    pathname = getpref('octviewer','seg_path_otb');
                else
                    addpref('octviewer','seg_path_otb','.');
                end
            case 'mat'
                fileext = '*.mat';
                if ispref('octviewer','seg_path_mat')
                    pathname = getpref('octviewer','seg_path_mat');
                else
                    addpref('octviewer','seg_path_mat','.');
                end
        end
        [filename, pathname] = uigetfile(fileext,'Load segmentation file',pathname);
        
        if isequal(filename,0)
           return
        end
        filefull = fullfile(pathname, filename);
        
        % Return focus to figure
        figure(hf);
        
        switch tag
            case 'otb'
                setpref('octviewer','seg_path_otb',pathname)
                
                control_pts = readControlPoints(filefull,[],false);
                
                for i = 1:size(control_pts,1)
                    for j = 1:size(control_pts,2)
                        % OCT toolbox indexing starts at 0
                        control_pts{i,j} = control_pts{i,j} + 1;
                    end
                end
            case 'mat'
                setpref('octviewer','seg_path_mat',pathname)
                
                S = load(filefull);
                
                if isfield(S,'bd_pts')
                    bd_pts = S.bd_pts;
                    
                    control_pts = cell(size(bd_pts,2),11);
                    for i = 1:size(bd_pts,2)
                        for j = 1:size(bd_pts,3)
                            control_pts{i,j} = [1:size(bd_pts,1); bd_pts(:,i,j)']';
                        end
                        control_pts{i,10} = [];
                        control_pts{i,11} = [];
                    end
                end
                
                % Segmentation gets centered on the data
                startx = round(size(octvol,2)/2) - floor((size(bd_pts,1)-1)/2);
                starty = round(size(octvol,3)/2) - floor((size(bd_pts,2)-1)/2);
        end
        
        updateseg();
    end
end

%%
function pts = interpolateCtrlPts(ctrl_pts)
% Interpolate a vector of control points to get values between the first and
% last point
    
    % Assume sorted
    xmin = ctrl_pts(1,1);
    xmax = ctrl_pts(end,1);
    
    xpts = (xmin:xmax)';
    
    % remove duplicate points 
    [pt1 m n] = unique(ctrl_pts(:,1));
    pt2 = ctrl_pts(m,2);
    ctrl_pts = [pt1, pt2];
    
%     pts = interp1(ctrl_pts(:,1),ctrl_pts(:,2),xpts,'cubic',nan);
    pts = cubicSplineInterp(ctrl_pts(:,1),ctrl_pts(:,2),xpts);
    pts = cat(2,xpts,pts);
end

%% 
function splPts = cubicSplineInterp(xpts,ypts,xptsOut)

    xpts = xpts-xpts(1)+1;
    n = length(xpts);

    % X = calculateNaturalCubicSpline(n-1,xpts);
    Y = calculateNaturalCubicSpline(n-1,ypts);

    splPts = zeros(size(xptsOut));
    for i = 1:n-1
        splPts(xpts(i)) = ypts(i);

        t = (((xpts(i)+1):(xpts(i+1)-1))-xpts(i))/(xpts(i+1)-xpts(i));
        splPts((xpts(i)+1):(xpts(i+1)-1)) = Y(i,1) + Y(i,2)*t + Y(i,3)*t.*t + Y(i,4)*t.*t.*t;

    %     for j = (xpts(i)+1):(xpts(i+1)-1)
    %         t = (j - xpts(i))/(xpts(i+1)-xpts(i));
    %         y = getSplineValue(Y(i,:),t);
    %         splPts(j) = y;
    %     end    
    end
    splPts(xpts(n)) = ypts(n);

end

function C = calculateNaturalCubicSpline(n,x)

    gamma = zeros(n+1,1);
    delta = zeros(n+1,1);
    D = zeros(n+1,1);

    gamma(1) = 0.5;
    for i = 2:n
        gamma(i) = 1/(4-gamma(i-1));
    end
    gamma(end) = 1/(2-gamma(end-1));

    delta(1) = 3*(x(2)-x(1))*gamma(1);
    for i = 2:n
        delta(i) = (3*(x(i+1)-x(i-1))-delta(i-1))*gamma(i); 
    end
    delta(end) = (3*(x(end)-x(end-1))-delta(end-1))*gamma(end);

    D(end) = delta(end);
    for i = n:-1:1
        D(i) = delta(i) - gamma(i)*D(i+1);
    end

    C = zeros(n,4);
    for i = 1:n
        C(i,:) = [x(i),D(i),3*(x(i+1)-x(i))-2*D(i)-D(i+1),2*(x(i)-x(i+1))+D(i)+D(i+1)];
    end

end


function v = getSplineValue(c,x)

    v = c(1) + c(2)*x + c(3)*x^2 + c(4)*x^3;

end