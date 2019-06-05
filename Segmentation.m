function varargout = Segmentation(varargin)
% SEGMENTATION M-file for Segmentation.fig
%      SEGMENTATION, by itself, creates a new SEGMENTATION or raises the existing
%      singleton*.
%
%      H = SEGMENTATION returns the handle to a new SEGMENTATION or the handle to
%      the existing singleton*.
%
%      SEGMENTATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SEGMENTATION.M with the given input arguments.
%
%      SEGMENTATION('Property','Value',...) creates a new SEGMENTATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before segmentation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to segmentation_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Author: Zeynep Akalin Acar, SCCN, 2008

% Copyright (C) 2007 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


% Edit the above text to modify the response to help segmentation

% Last Modified by GUIDE v2.5 01-Nov-2011 07:48:07

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @segmentation_OpeningFcn, ...
                   'gui_OutputFcn',  @segmentation_OutputFcn, ...
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

% ---------- Display Function -----------
function view_slice_Display(handles)
% This function displays handles.volume
vol = handles.volume;

% XXX mask = mask.*filteredvol;
if get(handles.MRradiobutton, 'Value') == 1
    vol = handles.volume;
elseif get(handles.Scalpradiobutton, 'Value') == 1
%    vol = handles.segm.scalpmask .* handles.filteredvol;
    vol = handles.segm.scalpmask;
elseif get(handles.Filteredradiobutton, 'Value') ==1
    vol = handles.filteredvol;
elseif get(handles.Brainradiobutton, 'Value') == 1
%    vol = handles.segm.brainmask .* handles.filteredvol;
    vol = handles.segm.brainmask;
elseif get(handles.OSradiobutton, 'Value') == 1
%    vol = handles.segm.outerskullmask .* handles.filteredvol;
    vol = handles.segm.outerskullmask;
elseif get(handles.ISradiobutton, 'Value') == 1
%    vol = handles.segm.innerskullmask .* handles.filteredvol;
    vol = handles.segm.innerskullmask;
end

sz = size(vol);

pos1 = ceil(get(handles.slider1,'Value'));
pos2 = ceil(get(handles.slider2,'Value'));
pos3 = ceil(get(handles.slider3,'Value'));
pos_val = sprintf('( %d, %d, %d )', pos3, pos1, pos2);
set(handles.text_xyz,'String',pos_val);

v1 = reshape(vol(pos1,:,:), sz(2), sz(3));
v2 = reshape(vol(:,pos2,:), sz(1), sz(3));
v3 = vol(:,:,pos3);

b1 = get(handles.axes1, 'ButtonDownFcn');
b2 = get(handles.axes2, 'ButtonDownFcn');
b3 = get(handles.axes3, 'ButtonDownFcn');
col = [0 1 0]; % green

imagesc(v1, 'Parent',handles.axes1); colormap pink;
set(handles.axes1,'YDir','normal');  % new 12/16/2008
line([pos3 pos3], [1 sz(2)],'Parent',handles.axes1, 'color', col);
line([1 sz(3)], [pos2 pos2],'Parent',handles.axes1, 'color', col);

imagesc(v2,'Parent',handles.axes2); colormap pink;
line([pos3 pos3], [1 sz(1)],'Parent',handles.axes2, 'color', col);
line([1 sz(3)], [pos1 pos1],'Parent',handles.axes2, 'color', col);

imagesc(v3,'Parent',handles.axes3); colormap pink;
line([1 sz(2)], [pos1 pos1], 'Parent',handles.axes3, 'color', col);
line([pos2 pos2], [1 sz(1)], 'Parent',handles.axes3, 'color', col);

set(handles.axes1, 'ButtonDownFcn', b1);
set(get(handles.axes1,'Children'), 'ButtonDownFcn',b1);
set(handles.axes2, 'ButtonDownFcn', b2);
set(get(handles.axes2,'Children'), 'ButtonDownFcn',b2);
set(handles.axes3, 'ButtonDownFcn', b3);
set(get(handles.axes3,'Children'), 'ButtonDownFcn',b3);

set(handles.text33,'String','y');
set(handles.text32,'String','x');
set(handles.text35,'String','z');
set(handles.text37,'String','y');
set(handles.text36,'String','x');
set(handles.text34,'String','z');

% --- Executes just before segmentation is made visible.
function segmentation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to segmentation (see VARARGIN)

% Parse arguments and set handles as necessary
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'subjectdir')
        i = i + 1;
        handles.OutputFolder = varargin{i};
    elseif strcmp(varargin{i}, 'subject')
        i = i + 1;
        handles.arg_subject = varargin{i};
    elseif strcmp(varargin{i}, 'session')
        i = i + 1;
        handles.arg_session = varargin{i};
    end
end

% Choose default command line output for segmentation
handles.output = hObject;
guidata(hObject, handles);
if isfield(handles,'OutputFolder')
    set(handles.textFolder, 'String', handles.OutputFolder);
end
%if strcmp(get(handles.figure1,'Visible'),'off')
%    view_slice_Display(handles);
%end

function handles = segmentation_SetVolume(handles, volume)
% Sets a volume for display
handles.volume = volume;
sz = size(handles.volume);
pos = ceil(sz/2);

set(handles.slider1, 'Max', sz(1));
set(handles.slider1, 'SliderStep', [1/sz(1) 1/sz(1)]);
set(handles.slider2, 'Max', sz(2));
set(handles.slider2, 'SliderStep', [1/sz(2) 1/sz(2)]);
set(handles.slider3, 'Max', sz(3));
set(handles.slider3, 'SliderStep', [1/sz(3) 1/sz(3)]);

set(handles.slider1,'Value', pos(1));
set(handles.slider2,'Value', pos(2));
set(handles.slider3,'Value', pos(3));

% set segmented volumes to 0
handles.segm.scalpmask = [];
handles.segm.brainmask = [];
handles.segm.outerskullmask = [];
handles.segm.innerskullmask = [];
handles.filteredvol = [];


% initialize current operation
handles.CurrentOperation = 1;
setup_operation(handles);

% Update handles structure
guidata(handles.figure1, handles);

update_view_selector(handles);





function update_view_selector(handles)

if not(isempty(handles.volume))
    set(handles.MRradiobutton, 'Enable', 'on')
else
    set(handles.MRradiobutton, 'Enable', 'off')
end

if not(isempty(handles.filteredvol))
    set(handles.Filteredradiobutton, 'Enable', 'on')
    set(handles.pushbuttonSaveFiltered, 'Enable', 'on')
else
    set(handles.Filteredradiobutton, 'Enable', 'off')
    set(handles.pushbuttonSaveFiltered, 'Enable', 'off')
end

if not(isempty(handles.segm.scalpmask))
    set(handles.Scalpradiobutton, 'Enable', 'on')
    set(handles.pushbuttonSaveSegm, 'Enable', 'on')
else
    set(handles.Scalpradiobutton, 'Enable', 'off')
    set(handles.pushbuttonSaveSegm, 'Enable', 'off')
end

if not(isempty(handles.segm.brainmask))
    set(handles.Brainradiobutton, 'Enable', 'on')
else
    set(handles.Brainradiobutton, 'Enable', 'off')
end

if not(isempty(handles.segm.outerskullmask))
    set(handles.OSradiobutton, 'Enable', 'on')
else
    set(handles.OSradiobutton, 'Enable', 'off')
end

if not(isempty(handles.segm.innerskullmask))
    set(handles.ISradiobutton, 'Enable', 'on')
else
    set(handles.ISradiobutton, 'Enable', 'off')
end
    

% --- Outputs from this function are returned to the command line.
function varargout = segmentation_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axes1);
cla;

popup_sel_index = get(handles.popupmenu1, 'Value');
switch popup_sel_index
    case 1
        plot(rand(5));
    case 2
        plot(sin(1:0.01:25.99));
    case 3
        bar(1:.5:10);
    case 4
        plot(membrane);
    case 5
        surf(peaks);
end


% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, path] = uigetfile('*.hdr');
if ~isequal(file, 0) && length(file) > 5
    handles.data.filename = file(1:length(file)-4);
    handles.data.filepath = path;
    handles.parameters.MRfile = file;
    handles.parameters.MRpath = path;
    % remove .hdr extension
    file = [path file(1:length(file)-4)];
    [x,y,z] = segm_readanalyze(file);
    set(handles.pbinhomog,'Enable','on');
    if get(handles.checkboxLRflip, 'Value') == 1
        % flip image left-right (image is flipped during MR acquisition)
        % image is saggital
        [K,L,M]=size(x); x1=zeros(K,L,M);
        for i=1:M
            x1(:,:,i)=x(:,:,M-i+1);
        end
        x=x1; clear x1;
    end
    handles.parameters.LRflip = get(handles.checkboxLRflip,'Value');
    handles = segmentation_SetVolume(handles, x);
    view_slice_Display(handles);
end


% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1


% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc
    set(hObject,'BackgroundColor','white');
else
    set(hObject,'BackgroundColor',get(0,'defaultUicontrolBackgroundColor'));
end

set(hObject, 'String', {'plot(rand(5))', 'plot(sin(1:0.01:25))', 'bar(1:.5:10)', 'plot(membrane)', 'surf(peaks)'});



% --- Executes during object creation, after setting all properties.
function axes2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes2


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

view_slice_Display(handles);


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
view_slice_Display(handles);


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
view_slice_Display(handles);


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes during object creation, after setting all properties.
function figure1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




% --- Executes on mouse press over axes background.
function axes1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=ceil(get(handles.axes1, 'CurrentPoint'));

%set(handles.slider1,'Value', pos(1));
set(handles.slider2,'Value', x(1,2));
set(handles.slider3,'Value', x(1,1));
view_slice_Display(handles);




% --- Executes on mouse press over axes background.
function axes2_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=ceil(get(handles.axes2, 'CurrentPoint'));

%set(handles.slider1,'Value', pos(1));
set(handles.slider1,'Value', x(1,2));
set(handles.slider3,'Value', x(1,1));
view_slice_Display(handles);


% --- Executes on mouse press over axes background.
function axes3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to axes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
x=ceil(get(handles.axes3, 'CurrentPoint'));

%set(handles.slider1,'Value', pos(1));
set(handles.slider1,'Value', x(1,2));
set(handles.slider2,'Value', x(1,1));
view_slice_Display(handles);




% --- Executes during object creation, after setting all properties.
function text1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3


% --- Executes on button press in checkbox4.
function checkbox4_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox4


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% --- Executes on button press in Runbutton.
function Runbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Runbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.CurrentOperation == 1
    iter = str2num(get(handles.editNumberofIter, 'String'));
%    ts   = str2num(get(handles.editTimeStep, 'String'));
    ts = 0.0625;
    cond = str2num(get(handles.editConductance, 'String'));
    handles.parameters.filter.iter = iter;
    handles.parameters.filter.cond = cond;
    set(handles.textStatus,'String','Filtering...'); pause(1);
    handles.filteredvol = segm_aniso_filtering(handles.volume, iter, ts, cond);
    set(handles.textStatus,'String','Image is filtered!'); pause(0.1);
elseif handles.CurrentOperation == 2
    % scalp
    set(handles.textStatus,'String','Segmenting scalp...'); pause(1);
    handles.segm.scalpmask = segm_scalp(handles.filteredvol);
    set(handles.textStatus,'String','Scalp segmented!')
elseif handles.CurrentOperation == 3
    % brain
    sli = str2num(get(handles.editSlice, 'String'));
    WMp(1) = str2num(get(handles.editWMp1, 'String'));
    WMp(2) = str2num(get(handles.editWMp2, 'String'));
    WMp(3) = str2num(get(handles.editWMp3, 'String'));
    sl = str2num(get(handles.editSetLevel, 'String'));
    st = str2num(get(handles.editSetThres, 'String'));
    handles.parameters.brain.slice = sli;
    handles.parameters.brain.WMp = WMp;
    handles.parameters.brain.filllevel = sl;
    handles.parameters.brain.threshold = st;
    set(handles.textStatus,'String','Segmenting brain...'); pause(1);
    handles.segm.brainmask = segm_brain(handles.filteredvol,handles.segm.scalpmask, sli, WMp, sl, st);
    set(handles.textStatus,'String','Brain segmented!')
elseif handles.CurrentOperation == 4
    % outer skull
    sli_eyes = str2num(get(handles.editSliceForEyes, 'String'));
    handles.parameters.skull.sli_eyes = sli_eyes;
    set(handles.textStatus,'String','Segmenting skull...'); pause(1);
    [handles.segm.outerskullmask, handles.X_dark, thr] = segm_outer_skull(handles.filteredvol, handles.segm.scalpmask, handles.segm.brainmask, sli_eyes);
    handles.parameters.skull.thr = thr; % 1/29/2013
    set(handles.textStatus,'String','Skull segmented!')
elseif handles.CurrentOperation == 5
    % inner skull
    set(handles.textStatus,'String','Segmenting CSF...'); pause(1);
    WMp(1) = str2num(get(handles.editWMp1, 'String'));
    WMp(2) = str2num(get(handles.editWMp2, 'String'));
    WMp(3) = str2num(get(handles.editWMp3, 'String'));
    handles.segm.innerskullmask = segm_inner_skull(handles.filteredvol, handles.segm.outerskullmask, handles.X_dark, handles.segm.brainmask,WMp);
    set(handles.textStatus,'String','Correcting skull and scalp...'); pause(1);
    [handles.segm.scalpmask, handles.segm.outerskullmask]=segm_final_skull(handles.segm.scalpmask, handles.segm.outerskullmask, handles.segm.innerskullmask, WMp);
    set(handles.textStatus,'String','Segmentation complete!'); pause(1);
end


% Enable and select Next button
set(handles.Nextbutton, 'Enable', 'on');
uicontrol(handles.Nextbutton);
    
% Update handles structure
guidata(handles.figure1, handles);

update_view_selector(handles);

function setup_operation(handles)
% setup current operation
col_on = [1 0 0];
col_off = [0 0 0.322];
set(handles.Prevbutton, 'Enable', 'on');
set(handles.Nextbutton, 'Enable', 'off');
if handles.CurrentOperation == 1
    set(handles.Prevbutton, 'Enable', 'off');
    set(handles.Filttext2, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'TooltipString','Run anisotropic filtering');
else
    set(handles.Filttext2, 'ForegroundColor',col_off);
end

if handles.CurrentOperation == 2
    set(handles.Scalptext2, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'TooltipString','Run scalp segmentation');
else
    set(handles.Scalptext2, 'ForegroundColor',col_off);
end
    
if handles.CurrentOperation == 3
    set(handles.Braintext2, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'TooltipString','Run brain segmentation');
else
    set(handles.Braintext2, 'ForegroundColor',col_off);
end

if handles.CurrentOperation == 4
    set(handles.OStext2, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'TooltipString','Run outer skull segmentation');
else
    set(handles.OStext2, 'ForegroundColor',col_off);
end

if handles.CurrentOperation == 5
    set(handles.IStext2, 'ForegroundColor',col_on);
    set(handles.Runbutton, 'TooltipString','Run inner skull segmentation');
else
    set(handles.IStext2, 'ForegroundColor',col_off);
end

uicontrol(handles.Runbutton);


% --------------------------------------------------------------------
function uipanel3_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% get(handles.uipanel3);
view_slice_Display(handles);


% --- Executes on button press in Prevbutton.
function Prevbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Prevbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.CurrentOperation > 1
    handles.CurrentOperation = handles.CurrentOperation - 1;
    setup_operation(handles);
    
    % Update handles structure
    guidata(handles.figure1, handles);
end


% --- Executes on button press in Nextbutton.
function Nextbutton_Callback(hObject, eventdata, handles)
% hObject    handle to Nextbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.CurrentOperation < 5
    handles.CurrentOperation = handles.CurrentOperation + 1;
    setup_operation(handles);
    
    % Update handles structure
    guidata(handles.figure1, handles);
end



function editNumberofIter_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofIter as text
%        str2double(get(hObject,'String')) returns contents of editNumberofIter as a double


% --- Executes during object creation, after setting all properties.
function editNumberofIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editTimeStep_Callback(hObject, eventdata, handles)
% hObject    handle to editTimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editTimeStep as text
%        str2double(get(hObject,'String')) returns contents of editTimeStep as a double


% --- Executes during object creation, after setting all properties.
function editTimeStep_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editTimeStep (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editConductance_Callback(hObject, eventdata, handles)
% hObject    handle to editConductance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editConductance as text
%        str2double(get(hObject,'String')) returns contents of editConductance as a double


% --- Executes during object creation, after setting all properties.
function editConductance_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editConductance (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editWMp1_Callback(hObject, eventdata, handles)
% hObject    handle to editWMp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWMp1 as text
%        str2double(get(hObject,'String')) returns contents of editWMp1 as a double


% --- Executes during object creation, after setting all properties.
function editWMp1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWMp1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWMp2_Callback(hObject, eventdata, handles)
% hObject    handle to editWMp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWMp2 as text
%        str2double(get(hObject,'String')) returns contents of editWMp2 as a double


% --- Executes during object creation, after setting all properties.
function editWMp2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWMp2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editWMp3_Callback(hObject, eventdata, handles)
% hObject    handle to editWMp3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editWMp3 as text
%        str2double(get(hObject,'String')) returns contents of editWMp3 as a double


% --- Executes during object creation, after setting all properties.
function editWMp3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editWMp3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSlice_Callback(hObject, eventdata, handles)
% hObject    handle to editSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSlice as text
%        str2double(get(hObject,'String')) returns contents of editSlice as a double


% --- Executes during object creation, after setting all properties.
function editSlice_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbuttonGetSlice.
function pushbuttonGetSlice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGetSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos2 = ceil(get(handles.slider2,'Value'));
set(handles.editSlice, 'String', num2str(pos2));


% --- Executes on button press in pushbuttonGetWMp.
function pushbuttonGetWMp_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGetWMp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

pos1 = ceil(get(handles.slider1,'Value'));
pos2 = ceil(get(handles.slider2,'Value'));
pos3 = ceil(get(handles.slider3,'Value'));
set(handles.editWMp1, 'String', num2str(pos1));
set(handles.editWMp2, 'String', num2str(pos2));
set(handles.editWMp3, 'String', num2str(pos3));




function editSliceForEyes_Callback(hObject, eventdata, handles)
% hObject    handle to editSliceForEyes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSliceForEyes as text
%        str2double(get(hObject,'String')) returns contents of editSliceForEyes as a double


% --- Executes during object creation, after setting all properties.
function editSliceForEyes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSliceForEyes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonGetEyeSlice.
function pushbuttonGetEyeSlice_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGetEyeSlice (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
pos2 = ceil(get(handles.slider2,'Value'));
set(handles.editSliceForEyes, 'String', num2str(pos2));




% --- Executes on button press in pushbuttonSaveFiltered.
function pushbuttonSaveFiltered_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveFiltered (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if isfield(handles, 'arg_subject')
    if isempty(handles.arg_subject)
        f = handles.data.filename; % if subject name is not given save the output 
    else                           % with MR image's name.
        f = handles.arg_subject;
    end
else
    f = handles.data.filename;
end
    
set(handles.textStatus,'String',['Saving filtered image as ' f '_filtered.mat']); pause(0.5);

filt_im = handles.filteredvol;
% convert to axial (from saggital) before saving!
ax = filt_im; clear filt_im
[K,L,M] = size(ax);
for i = 1:M
   filt_im(:,i,:) = (reshape(ax(:,:,i), K, L));
end
[K,L,M] = size(filt_im);
for i = 1:M
   filt_im(:,:,i) = rot90(filt_im(:,:,i),3);
end


OutputFolder = handles.OutputFolder;
if isempty(OutputFolder)
    error('You must enter output folder...')
end
lof = length(OutputFolder);
if OutputFolder(lof) ~= '/'
    OutputFolder(lof+1) = '/';
end
p = OutputFolder;


%save([p f '_filtered.mat'],'filt_im');

filt_im = filt_im/max(max(max(filt_im)));
[K,L,M] = size(filt_im);
mri.dim = [K L M];
mri.xgrid = [1:K];
mri.ygrid = [1:L];
mri.zgrid = [1:M];
mri.anatomy = filt_im;
mri.transform = eye(4);
mri.hdr = handles.parameters.MRfile;
save([p f '_mri'],'mri') 

clear filt_im ax;
set(handles.textStatus,'String',['Filtered image saved as ' f '_filtered.mat'])

% --- Executes on button press in pushbuttonSaveSegm.
function pushbuttonSaveSegm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSaveSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if isfield(handles, 'arg_subject')
    if isempty(handles.arg_subject)
        f = handles.data.filename; % if subject name is not given save the output 
    else                           % with MR image's name.
        f = handles.arg_subject;
    end
else
    f = handles.data.filename;
end

OutputFolder = handles.OutputFolder;
if isempty(OutputFolder)
    error('You must enter output folder...')
end
lof = length(OutputFolder);
if OutputFolder(lof) ~= '/'
    OutputFolder(lof+1) = '/';
end
p = OutputFolder;

set(handles.textStatus,'String',['Saving segmentation as ' f '_segments.mat']); pause(0.1);


%p = handles.data.filepath; % output folder will be entered
Segm = handles.segm;
% convert to axial (from saggital) before saving!
ax = Segm.scalpmask; 
[K,L,M] = size(ax);
for i = 1:M
   Sca(:,i,:) = (reshape(ax(:,:,i), K, L));
end
Segm.scalpmask = Sca;
clear Sca

[K1,L1,M1] = size(Segm.scalpmask);
for i=1:M1;    Segm.scalpmask(:,:,i) = rot90(Segm.scalpmask(:,:,i),3); end

ax = Segm.brainmask; 
if ~isempty(ax)
    for i = 1:M
        Sca(:,i,:) = (reshape(ax(:,:,i), K, L));
    end
    Segm.brainmask = Sca;
    clear Sca
    for i=1:M1;    Segm.brainmask(:,:,i) = rot90(Segm.brainmask(:,:,i),3); end
end

ax = Segm.outerskullmask; 
if ~isempty(ax)
    for i = 1:M
        Sca(:,i,:) = (reshape(ax(:,:,i), K, L));
    end
    Segm.outerskullmask = Sca;
    clear Sca
    for i=1:M1;    Segm.outerskullmask(:,:,i) = rot90(Segm.outerskullmask(:,:,i),3); end
end

ax = Segm.innerskullmask; 
if ~isempty(ax)
    for i = 1:M
        Sca(:,i,:) = (reshape(ax(:,:,i), K, L));
    end
    Segm.innerskullmask = Sca;
    clear Sca
    for i=1:M1;    Segm.innerskullmask(:,:,i) = rot90(Segm.innerskullmask(:,:,i),3); end

end

Segm.parameters = handles.parameters;

save([p f '_segments.mat'],'Segm');
clear Segm;
set(handles.textStatus,'String',['Segmentation saved as ' f '_segments.mat']);



function editSetLevel_Callback(hObject, eventdata, handles)
% hObject    handle to editSetLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSetLevel as text
%        str2double(get(hObject,'String')) returns contents of editSetLevel as a double


% --- Executes during object creation, after setting all properties.
function editSetLevel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSetLevel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSetThres_Callback(hObject, eventdata, handles)
% hObject    handle to editSetThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSetThres as text
%        str2double(get(hObject,'String')) returns contents of editSetThres as a double


% --- Executes during object creation, after setting all properties.
function editSetThres_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSetThres (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function editFolder_Callback(hObject, eventdata, handles)
% hObject    handle to editFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editFolder as text
%        str2double(get(hObject,'String')) returns contents of editFolder as a double


% --- Executes during object creation, after setting all properties.
function editFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.OutputFolder = uigetdir;
set(handles.textFolder, 'String', handles.OutputFolder);
% Update handles structure
guidata(handles.figure1, handles);


% --- Executes on button press in checkboxLRflip.
function checkboxLRflip_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLRflip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLRflip


% --- Executes on button press in pbinhomog.
function pbinhomog_Callback(hObject, eventdata, handles)
% hObject    handle to pbinhomog (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

vol = handles.volume;
[K,L,M] = size(vol);
I0 = vol(:,:,round(M/2));
set(handles.Runbutton, 'Enable', 'off');
set(handles.textStatus,'String','Calculating inhomogeneity...'); pause(1);
[B, B0, Ic, ImaskClean] = segm_inhomogeneity_correction(I0);
set(handles.Runbutton, 'Enable', 'on');
set(handles.textStatus,'String','Inhomogeneity calculated'); pause(0.5);
(var(B(:)))
if var(B(:))>0.06
    msgbox('Inhomogeneity in this image may result in incorrect segmentation. Please load a corrected volume. The NFT manual gives instructions for using Freesurfer inhomogeneity correction.','Inhomogeneity check')
else
    msgbox('No need for inhomogeneity correction!','Inhomogeneity check')
end


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function axes1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: place code in OpeningFcn to populate axes1


% --------------------------------------------------------------------
function uipanel3_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
