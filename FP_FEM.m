function varargout = FP_FEM(varargin)
% FP_FEM M-file for FP_FEM.fig
%      FP_FEM, by itself, creates a new FP_FEM or raises the existing
%      singleton*.
%
%      H = FP_FEM returns the handle to a new FP_FEM or the handle to
%      the existing singleton*.
%
%      FP_FEM('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FP_FEM.M with the given input arguments.
%
%      FP_FEM('Property','Value',...) creates a new FP_FEM or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Forward_Problem_Solution_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to FP_FEM_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% Author: Zeynep Akalin Acar, SCCN, 2010

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


% Edit the above text to modify the response to help FP_FEM

% Last Modified by GUIDE v2.5 09-Dec-2010 16:34:45

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @FP_FEM_OpeningFcn, ...
                   'gui_OutputFcn',  @FP_FEM_OutputFcn, ...
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


% --- Executes just before FP_FEM is made visible.
function FP_FEM_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to FP_FEM (see VARARGIN)


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


handles.session_changed = 0;


if isfield(handles,'OutputFolder')
    
    % change dir
    path = handles.OutputFolder;
    lof = length(path);
    if path(lof) ~= '/';        path(lof+1) = '/';    end;
    cd(path) % change directory
end    
%else % if no output folder is specified
    % XXX check!
    
    % Choose default command line output for FP_FEM
    handles.bemmesh = [];
    handles.session = [];
    handles.session_changed = 0;
    handles.session = [];
%end

handles.output = hObject;

update_display(handles);

if isfield(handles, 'arg_subject') && ~isempty(handles.arg_subject)
    set(handles.editModelName,'String',handles.arg_subject);
    set_session_changed(handles);
end

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes FP_FEM wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = FP_FEM_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function editMeshName_Callback(hObject, eventdata, handles)
% hObject    handle to editMeshName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editMeshName as text
%        str2double(get(hObject,'String')) returns contents of editMeshName as a double


% --- Executes during object creation, after setting all properties.
function editMeshName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editMeshName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function set_session_changed(handles)
% call when any change is made to the session edit boxes
% if (handles.session_changed ~= 1)
%     set(handles.sessionProgressText,'String','Value Changed!');
%     set(handles.pushbuttonCreateModel,'Enable','on');
% end
% handles.session_changed = 1;

set(handles.pushbuttonCreateModel,'Enable','off');
name = get(handles.editModelName,'String');
cond1 = str2num(get(handles.editScalpCond, 'String'));
cond2 = str2num(get(handles.editSkullCond, 'String'));
cond3 = str2num(get(handles.editCSFCond, 'String'));
cond4 = str2num(get(handles.editBrainCond, 'String'));

if isempty(get(handles.editMeshName,'String'))
    set(handles.sessionProgressText,'String', 'Please load a mesh.');
elseif isempty(name)
    set(handles.sessionProgressText,'String', 'Please enter session name.');
elseif ~isfield(handles, 'sensors')
    set(handles.sessionProgressText,'String', 'Please load sensors.');
elseif (isempty(cond1) || isempty(cond2) || isempty(cond4))
    set(handles.sessionProgressText,'String', 'Please Enter Conductivities');
else
    set(handles.sessionProgressText,'String', 'Ready to create Model');
    set(handles.pushbuttonCreateModel,'Enable','on');
end

guidata(handles.figure1, handles);


function editModelName_Callback(hObject, eventdata, handles)
% hObject    handle to editModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editModelName as text
%        str2double(get(hObject,'String')) returns contents of editModelName as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editModelName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editScalpCond_Callback(hObject, eventdata, handles)
% hObject    handle to editScalpCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editScalpCond as text
%        str2double(get(hObject,'String')) returns contents of editScalpCond as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editScalpCond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editScalpCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editCSFCond_Callback(hObject, eventdata, handles)
% hObject    handle to editCSFCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editCSFCond as text
%        str2double(get(hObject,'String')) returns contents of editCSFCond as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editCSFCond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editCSFCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editBrainCond_Callback(hObject, eventdata, handles)
% hObject    handle to editBrainCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editBrainCond as text
%        str2double(get(hObject,'String')) returns contents of editBrainCond as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editBrainCond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editBrainCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSkullCond_Callback(hObject, eventdata, handles)
% hObject    handle to editSkullCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSkullCond as text
%        str2double(get(hObject,'String')) returns contents of editSkullCond as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editSkullCond_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSkullCond (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in pushbuttonCreateModel.
function pushbuttonCreateModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCreateModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (isempty(handles.mesh_name))
    errordlg('Please Load a Mesh first');
    return
end
 if (handles.bemmesh.num_boundaries < 3)
     errordlg('Mesh must have at least 3 layers');
     return
 end
name = get(handles.editModelName,'String');
if (isempty(name))
    errordlg('Please Enter Session Name');
    return
end
 if (handles.bemmesh.num_boundaries < 3)
     errordlg('Mesh must have at least 3 layers');
     return
 end

 if handles.bemmesh.num_boundaries == 3
     set(handles.editCSFCond, 'String' ,'');
 end
%     
cond1 = str2num(get(handles.editScalpCond, 'String'));
cond2 = str2num(get(handles.editSkullCond, 'String'));
cond3 = str2num(get(handles.editCSFCond, 'String'));
cond4 = str2num(get(handles.editBrainCond, 'String'));

 if (handles.bemmesh.num_boundaries == 3 && ~isempty(cond3))
     errordlg('Mesh has no CSF layer');
     return
 end
cond = [cond1 cond2 cond3 cond4];

if (isempty(cond1) || isempty(cond2) || isempty(cond4))
    errordlg('Please Enter Scalp, Skull and Brain Conductivities');
    return
end

% if (handles.mesh.num_boundaries == 4 && isempty(cond3))
%     errordlg('Please Enter CSF Conductivity');
%     return
% end
mesh_name = handles.mesh_name;
vol2 = metufem_set_mesh(mesh_name);

%handles.model = bem_create_model(name, handles.mesh, cond, mod);
set(handles.sessionProgressText,'String','Generating FEM matrix...'); pause(0.5);

if isfield(handles,'MeshFolder')
    of = handles.MeshFolder; % Output Folder
else
    of=pwd;
end

if isempty(of)
    of = pwd;
end

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
    handles.MeshFolder = of;
end

% vol = metufem_set_mesh(handles.mesh_name);
% vol.cond = cond;
% vol.type = 'metufem';
% save([handles.arg_subject '_vol.mat'], 'vol');

sens.pnt = handles.sensors;
sens.type = 'eeg';
sens = metufem_calcrf(vol2, sens, of, cond);
set(handles.sessionProgressText,'String','FEM Session Created');

handles.session.name = name;
handles.session.cond = cond;
handles.session.sens = sens;
handles.session.type = 'fem';
handles.session.vol = metufem_set_mesh([handles.mesh_path handles.mesh_name]);

% save session
msave.session = handles.session;
msave.mesh_name = handles.mesh_name;
msave.mesh_path = handles.mesh_path;
save([handles.session.name '.session'], '-STRUCT', 'msave')

vol = mesh2volstr(handles.bemmesh);
vol.cond = cond;
vol.type = 'metufem';
save([handles.arg_subject '_vol.mat'], 'vol');

set(handles.pushbuttonCreateModel,'Enable','off');
handles.session_changed = 0;
guidata(handles.figure1, handles);

function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.textComputeLFM,'String','Computing...');
pause(0.5);

mesh_name = handles.mesh_name;
metufem('setup', mesh_name, 'sensors.dat', '');
metufem('setrf', handles.session.sens.rf);

sourcespace = handles.dipoles;
LFM = metufem('pot', sourcespace');

%[handles.potentials, handles.session] = bem_solve_lfm_eeg(handles.session, handles.dipoles);
set(handles.textComputeLFM,'String','LFM Computed');
guidata(hObject, handles);
update_display(handles);

f = handles.session.name;
save([f '_LFM.mat'],'LFM');
clear LFM;


% --- Executes on button press in pushbuttonShowMesh.
function pushbuttonShowMesh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
coordt = handles.mesh.coord;
ma = mean(coordt);
coordt = coordt - ones(length(coordt),1) * ma; 
h = eeglab_plotmesh(handles.mesh.elem, coordt);
set(gcf, 'Name', 'Figure: Mesh', 'NumberTitle', 'off', 'Color', [0.925 0.957 1]);


function editNumberofLayers_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofLayers as text
%        str2double(get(hObject,'String')) returns contents of editNumberofLayers as a double


% --- Executes during object creation, after setting all properties.
function editNumberofLayers_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofLayers (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNumberofNodes_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofNodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofNodes as text
%        str2double(get(hObject,'String')) returns contents of editNumberofNodes as a double


% --- Executes during object creation, after setting all properties.
function editNumberofNodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofNodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNumberofElements_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofElements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofElements as text
%        str2double(get(hObject,'String')) returns contents of editNumberofElements as a double


% --- Executes during object creation, after setting all properties.
function editNumberofElements_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofElements (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editNumberofNE_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofNE as text
%        str2double(get(hObject,'String')) returns contents of editNumberofNE as a double


% --- Executes during object creation, after setting all properties.
function editNumberofNE_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofNE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonLoadSensors.
function pushbuttonLoadSensors_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadSensors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load from coordinates
    [file, path] = uigetfile('*.sens; *.sensors');
    if ~isequal(file, 0) && length(file) > 6
        if file(length(file)-3:length(file)) == 'sens'
            % remove .sens extension
            filen = [path file(1:length(file)-5)];
            handles.sensors = load([filen '.sens'], '-ascii');
%            fileind = [path file(1:length(file)-16)]; % clear headsensors.sens
%            handles.sensorindex = load([fileind 'sensorindex'], '-ascii');
            handles.sensorindex = load([filen '.sensorindex'], '-ascii');
            % Update handles structure
            guidata(hObject, handles);
            set(handles.editNumberofSensors,'String',length(handles.sensors));
        else
            sens = load([path file],'-mat');
            handles.sensors = sens.pnt;
            handles.sensorindex = sens.ind;
            % Update handles structure
            guidata(hObject, handles);
            set(handles.editNumberofSensors,'String',length(handles.sensors));
        end
    end
%set_session_changed(handles);    
%update_display(handles);
set_session_changed(handles);

%

function editNumberofSensors_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofSensors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofSensors as text
%        str2double(get(hObject,'String')) returns contents of editNumberofSensors as a double


% --- Executes during object creation, after setting all properties.
function editNumberofSensors_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofSensors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Convert sensMat to Electrode coordinates.
function C = sensMatToCoord(S, Coord)
% S      Sens Mat
% Coord  Mesh Coordinates
idx = unique(S(:, 1));
C = zeros(length(idx), 3);

for i = 1:length(idx)
    si = find(S(:,1) == idx(i));
    C(i,:) = S(si, 3)' * Coord(S(si,2),:);
end


% --- Executes on button press in pushbuttonLoadDipoles.
function pushbuttonLoadDipoles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadDipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[file, path] = uigetfile('*.dip');
if ~isequal(file, 0) && length(file) > 5
    handles.dipoles_name = file(1:length(file)-4);
    file = [path file(1:length(file)-4)];
    handles.dipoles = load([file '.dip'], '-ascii');
    % Update handles structure
    guidata(hObject, handles);
    set(handles.editNumberofdipoles,'String',size(handles.dipoles,1));
end
%update_display(handles);
%set_session_changed(handles);
set(handles.uipanelDipLoad,'visible','on')
set(handles.textComputeLFM,'String',' ');

% --------------------------------------------------------------------
function File_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to File_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Untitled_2_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% --------------------------------------------------------------------
function update_display(handles)
if (isempty(handles.bemmesh))
    set(handles.editMeshName,'String',[]);
    set(handles.editNumberofLayers,'String',[]);
    set(handles.editNumberofNodes,'String',[]);
%    set(handles.editNumberofElements,'String',[]);
    set(handles.editNumberofNE,'String',[]);
    set(handles.pushbuttonShowMesh,'Enable','off');
else
    set(handles.editMeshName,'String',handles.mesh_name);
    set(handles.editNumberofLayers,'String',handles.bemmesh.num_boundaries);
    set(handles.editNumberofNodes,'String',handles.mesh.num_nodes);
    %set(handles.editNumberofElements,'String',handles.mesh.num_elements);
    set(handles.editNumberofNE,'String',handles.mesh.num_node_elem);
    set(handles.pushbuttonShowMesh,'Enable','on');
end

if (isempty(handles.session))
    set(handles.editModelName,'String',[]);
    %set(handles.editScalpCond,'String',[]);
    %set(handles.editSkullCond,'String',[]);
    %set(handles.editCSFCond,'String',[]);
    %set(handles.editBrainCond,'String',[]);    
else
    set(handles.editModelName,'String',handles.session.name);
    set(handles.editScalpCond,'String',handles.session.cond(1));
    set(handles.editSkullCond,'String',handles.session.cond(2));
    if length(handles.session.cond) == 3
        set(handles.editCSFCond,'String',[]);
        set(handles.editBrainCond,'String',handles.session.cond(3));
    elseif length(handles.session.cond) == 4
        set(handles.editCSFCond,'String',handles.session.cond(3));
        set(handles.editBrainCond,'String',handles.session.cond(4));
    end
    set(handles.sessionProgressText,'String','FEM Session Loaded');
end
set_session_changed(handles);

if isfield(handles,'sensors')
    set(handles.editNumberofSensors,'String',length(handles.sensors));
end

if isfield(handles.bemmesh,'num_boundaries')
    if handles.bemmesh.num_boundaries == 4
        set(handles.uipanelCSF, 'visible', 'on')
    elseif handles.bemmesh.num_boundaries == 3
        set(handles.uipanelCSF, 'visible', 'off')
    end
end
    
% --------------------------------------------------------------------
function Load_Mesh_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Mesh_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.msh');
if ~isequal(file, 0) && length(file) > 5
    handles = load_mesh(handles, file, path);
end

% Update handles structure
guidata(hObject, handles);
update_display(handles);

function handles = load_mesh(handles, file, path)
    handles.mesh_name = file;
    handles.mesh_path = path;
    
    fid = fopen([path file]);
    nonodes = fgetl(fid);
    fclose(fid);
    handles.mesh.num_nodes = nonodes;
    handles.mesh.num_node_elem = 4; % Tetgen linear meshes
    
    % remove .msh extension and load the BEM mesh for layer information
    bemmeshfile = [path file(1:length(file)-6)];
    handles.bemmesh = bem_load_mesh(bemmeshfile);
    handles.session = [];

function handles = load_session(handles, file)
    msave = load(file, '-MAT');

    if ~isfield(msave.session, 'type') ||  ~strcmp(msave.session.type, 'fem')
        errordlg('Not a FEM session','Input Error');
        return;
    end
    
    % load the mesh
    handles = load_mesh(handles, msave.mesh_name, msave.mesh_path);

    % set the session
    handles.session = msave.session;
            
    % re-create sensors.dat
    sens = msave.session.sens;
    sens_name = 'sensors.dat';
    num_sens = size(sens.pnt, 1);
    sensors(:,2:4) = sens.pnt;
    sensors(:,1) = (1:num_sens)';

    fid = fopen(sens_name, 'w');
    fprintf(fid,'%d\n', num_sens);
    fprintf(fid,'%d %5.15f %5.15f %5.15f\r\n', sensors');
    fclose(fid);
    
    handles.sensors = sens.pnt;


% --------------------------------------------------------------------
function Load_Model_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Model_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.session');
if ~isequal(file, 0) && length(file) > 7
    file = [path file];
    handles = load_session(handles, file);
end

% Update handles structure
guidata(hObject, handles);
update_display(handles);


% --------------------------------------------------------------------


function editNumberofdipoles_Callback(hObject, eventdata, handles)
% hObject    handle to editNumberofdipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumberofdipoles as text
%        str2double(get(hObject,'String')) returns contents of editNumberofdipoles as a double


% --- Executes during object creation, after setting all properties.
function editNumberofdipoles_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumberofdipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --------------------------------------------------------------------
%function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%get(handles.uipanel5);



% --------------------------------------------------------------------
function uipanel6_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function uipanel6_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



% --------------------------------------------------------------------
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)









