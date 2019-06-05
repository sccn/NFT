function varargout = Coregistration(varargin)
% COREGISTRATION M-file for Coregistration.fig
%      COREGISTRATION, by itself, creates a new COREGISTRATION or raises the existing
%      singleton*.
%
%      H = COREGISTRATION returns the handle to a new COREGISTRATION or the handle to
%      the existing singleton*.
%
%      COREGISTRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in COREGISTRATION.M with the given input arguments.
%
%      COREGISTRATION('Property','Value',...) creates a new COREGISTRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Coregistration_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Coregistration_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Coregistration

% Last Modified by GUIDE v2.5 31-Mar-2010 20:55:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Coregistration_OpeningFcn, ...
                   'gui_OutputFcn',  @Coregistration_OutputFcn, ...
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


% --- Executes just before Coregistration is made visible.
function Coregistration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Coregistration (see VARARGIN)


% Parse arguments and set handles as necessary
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'subjectdir')
        i = i + 1;
        handles.MeshFolder = varargin{i};
    elseif strcmp(varargin{i}, 'subject')
        i = i + 1;
        handles.arg_subject = varargin{i};
    elseif strcmp(varargin{i}, 'session')
        i = i + 1;
        handles.arg_session = varargin{i};
    end
end

% Choose default command line output for Mesh_generation
handles.output = hObject;
if isfield(handles,'MeshFolder')
    set(handles.textMeshFolder, 'String', handles.MeshFolder);
end

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Coregistration wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Coregistration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbuttonMeshFolder.
function pushbuttonMeshFolder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMeshFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MeshFolder = uigetdir;
set(handles.textMeshFolder, 'String', handles.MeshFolder);
% Update handles structure
guidata(handles.figure1, handles);

if ~isempty(handles.MeshFolder) %& ~isempty(handles.elocfn)
    set(handles.pushbuttonManualCoreg, 'Enable', 'on')
end




% --------------------------------------------------------------------
% function OpenEloc_Callback(hObject, eventdata, handles)
% % hObject    handle to OpenEloc (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% [file, path] = uigetfile('*.elp');
% if ~isequal(file, 0) && length(file) > 1
%     handles.eloc = readlocs([path file]);   % subject's electrode locations
%     handles.elocfn = [path file];
%     handles.sensorpath = path;
%     handles.sensorname = file(1:length(file)-4);
% end
% % Update handles structure
% guidata(handles.figure1, handles);


% --- Executes on button press in pushbuttonManualCoreg.
function pushbuttonManualCoreg_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonManualCoreg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% manual coregistration
[C1,E1] = mesh_readsmf([handles.MeshFolder '/Scalp.smf'],0,0,0,1); % subject's scalp mesh

handles.Coord = C1;
handles.Elem = E1;

% subtract mean of the mesh for co-registration
ma_mesh = mean(C1(:,2:4));
C1(:,2:4) = C1(:,2:4) -  ones(length(C1),1) * ma_mesh;

eloc = handles.eloc;
ne = length(eloc);
for i = 1:ne
    elo(i,1) = eloc(i).X;
    elo(i,2) = eloc(i).Y;
    elo(i,3) = eloc(i).Z;
end
a1 = max(elo)-min(elo);
a2 = max(C1(:,2:4))-min(C1(:,2:4));
rat = mean(a2./a1);
% make the same scale with the mesh
if rat>500
    elo = elo * 1000; 
elseif rat>50
    elo = elo * 100; 
elseif rat>5
    elo = elo * 10;
end
ne = size(elo,1);

[d, elo] = warping_distafterwarping([0 0 0 0 0 90], elo, elo); % arrange orientation
%[d, elo] = warping_distafterwarping([0 0 0 15 0 0], elo, elo); % arrange orientation

% find coarse translation
a1 = max(elo);        b1 = min(elo);
a2 = max(C1(:,2:4));  b2 = min(C1(:,2:4));
d1 = (a2-b2-(a1-b1)) / 2;

tr = max(elo) - max(C1(:,2:4));
tr(1) = tr(1) + d1(1);
tr(2) = tr(2) + d1(2);
elo2 = elo - ones(length(elo),1) * tr;

eloc2 = eloc;
for i=1:ne
    eloc2(i).X = elo2(i,1);
    eloc2(i).Y = elo2(i,2);
    eloc2(i).Z = elo2(i,3);
end

mesh.TRI1 = E1(:,2:4);
mesh.POS = C1(:,2:4);
[ch, tr] = coregister(eloc2,handles.elocfn,'mesh',mesh,'autoscale','off');
set(handles.textTransformation,'String',num2str(tr(1:3)));
set(handles.text11,'String',num2str(tr(4:6)));

% add the mean
ch.pnt = ch.pnt +  ones(length(ch.pnt),1) * ma_mesh;

handles.ch = ch;
handles.tr = tr;

ind = warping_scalp_eloc_index(handles.ch.pnt, handles.Coord,1:ne);

elo2 = handles.ch.pnt(ind,:);

[elo1, dm1] = warping_distmeshafterwarping([0 0 0 0 0 0], elo2, handles.Coord, handles.Elem);
[elox, dm1] = warping_distmeshafterwarping([0 0 0 0 0 0], handles.ch.pnt, handles.Coord, handles.Elem);
h = plotmesh(handles.Elem(:,2:4),handles.Coord(:,2:4));hold
set(gcf, 'Name', 'Figure: Co-registered electrode locations (initial)', 'NumberTitle', 'off', 'Color', [0.925 0.956 1])
plot3(elo1(:,1),elo1(:,2),elo1(:,3),'r.'); view(180,0)
plot3(elox(:,1),elox(:,2),elox(:,3),'go'); view(180,0)

handles.init_sensors = elo1;
handles.ind = ind;
handles.init_param = tr;
% Update handles structure
guidata(handles.figure1, handles);

set(handles.pushbutton4, 'Enable', 'on')
set(handles.pushbutton7, 'Enable', 'on')

% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% automatic coregistration

ind = warping_scalp_eloc_index(handles.ch.pnt, handles.Coord, 1:length(handles.ch.pnt));

elo2 = handles.ch.pnt(ind,:);
options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun',1e-6);
Xo = [0 0 0 0 0 0];

% find the electrode locations that are close to the scalp to use with
% co-registration
[elox, dm] = warping_distmeshafterwarping(Xo, elo2, handles.Coord, handles.Elem);
mdm = median(dm); sdm = std(dm);
kdm = find((dm < mdm + sdm) & (dm > mdm - sdm)); % whos kdm
elo_tr = elo2(kdm,:);

set(handles.textStatus,'String','Computing translation and rotation parameters...'); pause(1);

% find translation and rotation for electrode positions - elo_tr
X = fminsearch(@(X) coreg_funrst(X, elo_tr, handles.Coord, handles.Elem), Xo, options);
set(handles.textStatus,'String','Translation and rotation parameters are computed.'); pause(1);

% apply translation and rotation to elo2
[elo3, dm1] = warping_distmeshafterwarping(X, elo2, handles.Coord, handles.Elem);
set(handles.textTrans2,'String',num2str(X(1:3)));
set(handles.text12,'String',num2str(X(4:6)));

h = plotmesh(handles.Elem(:,2:4),handles.Coord(:,2:4));hold
set(gcf, 'Name', 'Figure: Co-registered electrode locations', 'NumberTitle', 'off', 'Color', [0.925 0.956 1])
plot3(elo3(:,1),elo3(:,2),elo3(:,3),'r.'); view(180,0)

handles.comp_sensors = elo3;
handles.ind = ind;
handles.auto_param = X;
set(handles.pushbutton8, 'Enable', 'on')

% Update handles structure
guidata(handles.figure1, handles);


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



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.*');
if ~isequal(file, 0) && length(file) > 1
    handles.eloc = readlocs([path file]);   % subject's electrode locations
    sens_fn = [path file];
    handles.elocfn = [path file];
    handles.sensorpath = path;
    handles.sensorname = file(1:length(file)-4);
    
    p = handles.MeshFolder; % save the files in mesh folder
    lof = length(p);
    if p(lof) ~= '/';   p(lof+1) = '/'; end;
    save([p 'ori_sen_loc'], 'sens_fn'); % save the location of original sensors
    %varargout{1} = handles.elocfn;
end
set(handles.text13, 'String', handles.elocfn);
handles.init_param = [0 0 0 0 0 0];
handles.auto_param = [0 0 0 0 0 0];
% Update handles structure
guidata(handles.figure1, handles);


if isfield(handles,'MeshFolder') & ~isempty(handles.elocfn)
    set(handles.pushbuttonManualCoreg, 'Enable', 'on')
end



% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p = handles.MeshFolder; % save the files in mesh folder

lof = length(p);
if p(lof) ~= '/'
    p(lof+1) = '/';
end

ind = handles.ind;
fit_sensors = handles.init_sensors;


% save with session and subject name
if isfield(handles,'arg_session') & isfield(handles,'arg_subject')
    f = [handles.arg_subject '_' handles.arg_session];
else
    f = 'temp';
end
ssave.fn = handles.elocfn;
ssave.eloc = handles.eloc;
ssave.pnt = fit_sensors;
ssave.ind = ind;
ssave.param.init = handles.init_param;
ssave.param.auto = handles.auto_param;
save([p f '.sensors'], '-STRUCT', 'ssave')
set(handles.textStatus,'String','Initial registration is saved.'); pause(0.5);

%save([p f '_sensorindex'], 'ind', '-ascii'); % save the index, to use in IP
%save([p f '_headsensors.sens'], 'fit_sensors', '-ascii');

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

p = handles.MeshFolder; % save the files in mesh folder

lof = length(p);
if p(lof) ~= '/'
    p(lof+1) = '/';
end

ind = handles.ind;
fit_sensors = handles.comp_sensors;

% save with session and subject name
if isfield(handles,'arg_session') & isfield(handles,'arg_subject')
    f = [handles.arg_subject '_' handles.arg_session];
else
    f = 'temp';
end
ssave.fn = handles.elocfn;
ssave.eloc = handles.eloc;
ssave.pnt = fit_sensors;
ssave.ind = ind;
ssave.param.init = handles.init_param;
ssave.param.auto = handles.auto_param;
save([p f '.sensors'], '-STRUCT', 'ssave')
set(handles.textStatus,'String','Automatic registration is saved.'); pause(0.5);
%save([p f '_sensorindex'], 'ind', '-ascii'); % save the index, to use in IP
%save([p f '_headsensors.sens'], 'fit_sensors', '-ascii');
