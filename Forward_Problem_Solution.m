function varargout = Forward_Problem_Solution(varargin)
% FORWARD_PROBLEM_SOLUTION M-file for Forward_Problem_Solution.fig
%      FORWARD_PROBLEM_SOLUTION, by itself, creates a new FORWARD_PROBLEM_SOLUTION or raises the existing
%      singleton*.
%
%      H = FORWARD_PROBLEM_SOLUTION returns the handle to a new FORWARD_PROBLEM_SOLUTION or the handle to
%      the existing singleton*.
%
%      FORWARD_PROBLEM_SOLUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in FORWARD_PROBLEM_SOLUTION.M with the given input arguments.
%
%      FORWARD_PROBLEM_SOLUTION('Property','Value',...) creates a new FORWARD_PROBLEM_SOLUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Forward_Problem_Solution_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Forward_Problem_Solution_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Forward_Problem_Solution

% Last Modified by GUIDE v2.5 20-Mar-2008 12:14:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Forward_Problem_Solution_OpeningFcn, ...
                   'gui_OutputFcn',  @Forward_Problem_Solution_OutputFcn, ...
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


% --- Executes just before Forward_Problem_Solution is made visible.
function Forward_Problem_Solution_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Forward_Problem_Solution (see VARARGIN)


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
handles.model_changed = 0;


if isfield(handles,'OutputFolder')
    
    % change dir
    path = handles.OutputFolder;
    lof = length(path);
    if path(lof) ~= '/';        path(lof+1) = '/';    end;
    cd(path) % change directory
    % try to load mesh, model or session
    if ~isempty(handles.arg_subject) && ~isempty(handles.arg_session)
        % try to load session
        file = [handles.arg_session '.session'];
        try
            ssave = load(file, '-MAT');
            handles = load_model(hObject,handles, [path ssave.model_name, '.model']);
            handles.session = bem_create_session(ssave.name, handles.model, ssave.Smatrix);
            handles.session = bem_load_transfer_matrix(handles.session, 'tmte');
            handles.sensors = sensMatToCoord(handles.session.Smatrix, handles.mesh.coord);
        catch
            try
                % try to load model
                file = [handles.arg_subject '.model'];
                handles = load_model(hObject,handles, file);
                handles.session = [];
            catch
                try
                    % try to load mesh
                    file = handles.arg_subject;
                    handles.mesh = bem_load_mesh(file);
                    handles.model = [];
                    handles.session = [];
                catch
                    handles.mesh = [];
                    handles.model = [];
                    handles.model_changed = 0;
                    handles.session = [];
                    handles.session_changed = 0;
                    rethrow(lasterror);
                end
            end
        end
    elseif ~isempty(handles.arg_subject)
        try
            % try to load model
            file = [handles.arg_subject '.model'];
            handles = load_model(hObject,handles, file);
            handles.session = [];
        catch
            try
                % try to load mesh
                file = handles.arg_subject;
                handles.mesh = bem_load_mesh(file);
                handles.model = [];
                handles.session = [];
            catch
                handles.mesh = [];
                handles.model = [];
                handles.model_changed = 0;
                handles.session = [];
                handles.session_changed = 0;
            end
        end
    else % ~isfield(handles,'arg_subject') & ~isfield(handles,'arg_session')
        handles.mesh = [];
        handles.model = [];
        handles.model_changed = 0;
        handles.session = [];
        handles.session_changed = 0;
    end     %if isfield(handles, 'arg_subject') & isfield(handles,'arg_session')
    
    % try loading source space
    try
        file = 'sourcespace.dip';
%        handles.dipoles_name = file(1:length(file)-4);
        handles.dipoles.pos = load([file], '-ascii');
        handles.dipoles.sym = 0;
        set(handles.editNumberofdipoles,'String',size(handles.dipoles.pos,1));
    end
else % if no output folder is specified
    % XXX check!
    
    % Choose default command line output for Forward_Problem_Solution
    handles.mesh = [];
    handles.model = [];
    handles.model_changed = 0;
    handles.session = [];
    handles.session_changed = 0;
end

handles.output = hObject;

update_display(handles);

if isfield(handles, 'arg_subject') && ~isempty(handles.arg_subject)
    set(handles.editModelName,'String',handles.arg_subject);
    set_model_changed(handles);
end
if isfield(handles,'arg_session') && ~isempty(handles.arg_session)
    set(handles.editSessionName,'String',handles.arg_session);
    set_session_changed(handles);
end

% Update handles structure
guidata(hObject, handles);



% UIWAIT makes Forward_Problem_Solution wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Forward_Problem_Solution_OutputFcn(hObject, eventdata, handles) 
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


function set_model_changed(handles)
% call when any change is made to the model edit boxes
if (handles.model_changed ~= 1)
    set(handles.modelProgressText,'String','Value Changed!');
    set(handles.pushbuttonCreateModel,'Enable','on');
end
handles.model_changed = 1;
guidata(handles.figure1, handles);


function editModelName_Callback(hObject, eventdata, handles)
% hObject    handle to editModelName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editModelName as text
%        str2double(get(hObject,'String')) returns contents of editModelName as a double
set_model_changed(handles);

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

function set_session_changed(handles)
% call when any change is made to the model edit boxes
if (handles.session_changed ~= 1)
    set(handles.sessionProgressText, 'String', 'Value Changed!');
end
if (isempty(handles.session) && ~isempty(get(handles.editSessionName,'String'))...
        && isfield(handles,'sensors') && ~isempty(handles.model))
    set(handles.pushbuttonGenerateTM, 'Enable', 'on');
else
    set(handles.pushbuttonGenerateTM, 'Enable', 'off');
end
 
handles.session_changed = 1;
guidata(handles.figure1, handles);


function editSessionName_Callback(hObject, eventdata, handles)
% hObject    handle to editSessionName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSessionName as text
%        str2double(get(hObject,'String')) returns contents of editSessionName as a double
set_session_changed(handles);

% --- Executes during object creation, after setting all properties.
function editSessionName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSessionName (see GCBO)
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
set_model_changed(handles);

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
set_model_changed(handles);

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
set_model_changed(handles);

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
set_model_changed(handles);

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


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
set_model_changed(handles);

% --- Executes on button press in pushbuttonCreateModel.
function pushbuttonCreateModel_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonCreateModel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if (isempty(handles.mesh))
    errordlg('Please Load a Mesh first');
    return
end
if (handles.mesh.num_boundaries < 3)
    errordlg('Mesh must have at least 3 layers');
    return
end
name = get(handles.editModelName(),'String');
if (isempty(name))
    errordlg('Please Enter Model Name');
    return
end
if (handles.mesh.num_boundaries < 3)
    errordlg('Mesh must have at least 3 layers');
    return
end

if handles.mesh.num_boundaries == 3
    set(handles.editCSFCond, 'String' ,'');
end
    
cond1 = str2num(get(handles.editScalpCond, 'String'));
cond2 = str2num(get(handles.editSkullCond, 'String'));
cond3 = str2num(get(handles.editCSFCond, 'String'));
cond4 = str2num(get(handles.editBrainCond, 'String'));

if (handles.mesh.num_boundaries == 3 && ~isempty(cond3))
    errordlg('Mesh has no CSF layer');
    return
end
cond = [cond1 cond2 cond3 cond4];

if (isempty(cond1) || isempty(cond2) || isempty(cond4))
    errordlg('Please Enter Scalp, Skull and Brain Conductivities');
    return
end

if (handles.mesh.num_boundaries == 4 && isempty(cond3))
    errordlg('Please Enter CSF Conductivity');
    return
end

if (get(handles.checkbox1, 'Value') == 1)
    mod = 3;
else
	mod = -1;
end
handles.model = bem_create_model(name, handles.mesh, cond, mod);
set(handles.modelProgressText,'String','Generating matrices...');
pause(0.5);
bem_generate_eeg_matrices(handles.model);
set(handles.modelProgressText,'String','BEM Model Created');

% save model
msave.name = handles.model.name;
msave.mesh_name = handles.model.mesh.name;
msave.cond = handles.model.cond;
msave.mod = handles.model.mod;
save([handles.model.name '.model'], '-STRUCT', 'msave')

vol = mesh2volstr(handles.model.name);
vol.cond = handles.model.cond;
vol.type = 'metubem';
save([handles.arg_subject '_vol.mat'], 'vol');


set(handles.pushbuttonCreateModel,'Enable','off');
handles.model_changed = 0;
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


% --- Executes on button press in pushbuttonGenerateTM.
function pushbuttonGenerateTM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonGenerateTM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if (isempty(handles.mesh))
    errordlg('Please load a mesh first');
    return
end
if (isempty(handles.model))
    errordlg('Please load or generate a model first');
    return
end
name = get(handles.editSessionName,'String');
if (isempty(name))
    errordlg('Please Enter Session Name');
    return
end

if (isempty(handles.sensors))
    errordlg('Please Load a Sensor file');
    return
end

a = dir(['ori_sen_loc.mat']);
if size(a,1) > 0
    load ori_sen_loc
    handles.elocfn = sens_fn;
    eloc = readlocs(sens_fn);
end

sens.pnt = handles.sensors;
nsens = size(handles.sensors,1);
ind = handles.sensorindex;
% 
if isfield(handles, 'sensorindex') & isfield(handles,'elocfn')
    for ii=1:nsens
        sens.label{ii} = eloc(ind(ii)).labels;
    end
end


if isvector(handles.sensors)
    Smatrix = bem_smatrix_from_nodes(handles.mesh, handles.sensors);
elseif size(handles.sensors,2) == 3
    Smatrix = bem_smatrix_from_coordinates(handles.mesh, handles.sensors);
end

handles.session = bem_create_session(name, handles.model, Smatrix);
set(handles.sessionProgressText,'String','Generating matrices...');

handles.session = bem_generate_eeg_transfer_matrix(handles.session);
set(handles.sessionProgressText,'String','Session Created');

% save session
ssave.name = handles.session.name;
ssave.model_name = handles.session.model.name;
ssave.Smatrix = Smatrix;
ssave.sens = sens;
%ssave.tmte = handles.session.tmte;
save([handles.session.name '.session'], '-STRUCT', 'ssave');

set(handles.pushbuttonGenerateTM,'Enable','off');
handles.session_changed = 0;
guidata(handles.figure1, handles);

update_display(handles);

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
[handles.potentials, handles.session] = bem_solve_lfm_eeg(handles.session, handles.dipoles.pos);
set(handles.textComputeLFM,'String','LFM Computed');
guidata(hObject, handles);
update_display(handles);

f = handles.session.name;
if handles.dipoles.sym == 0
    LFM = handles.potentials;
    save([f '_LFM.mat'],'LFM');
    clear LFM;
else % symmetric dipoles
    Np = size(handles.dipoles.pos,1);
    LFM1 = handles.potentials(:,1:Np/2);
    LFM2 = handles.potentials(:,Np/2+1:Np);
    LFM = LFM1+ LFM2;
    save([f '_sLFM.mat'],'LFM');
    clear LFM*;
end


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

% XXX control node list
if get(handles.radiobuttonNodeList, 'Value') == 1
    [file, path] = uigetfile('*.dat');
    if ~isequal(file, 0) && length(file) > 5
        % remove .dat extension
        file = [path file(1:length(file)-4)];
        handles.sensors = load([file '.dat'], '-ascii');
        % Update handles structure
        guidata(hObject, handles);
        set(handles.editNumberofSensors,'String',length(handles.sensors));
    end
% XXX control coordinate list   1. size =(Nx3) 2. on the mesh?
else % load from coordinates
    [file, path] = uigetfile('*.sens; *.sensors');
    if ~isequal(file, 0) && length(file) > 6
        if file(length(file)-3:length(file)) == 'sens'
            % remove .sens extension
            filen = [path file(1:length(file)-5)];
            handles.sensors = load([filen '.sens'], '-ascii');
            fileind = [path file(1:length(file)-16)]; % clear headsensors.sens
            handles.sensorindex = load([fileind 'sensorindex'], '-ascii');
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
end
set_session_changed(handles);    
%update_display(handles);
%set_model_changed(handles);
set(handles.uipanelSensLoad,'visible','on')
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

% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles,'sensors')
    handles.sensors = sensMatToCoord(handles.session.Smatrix, handles.mesh.coord);
end

if size(handles.sensors, 2) == 3
    eloc = handles.sensors;
else
    eloc = handles.mesh.coord(handles.sensors,:);
end

dipnum = str2num(get(handles.editDipoleNumber, 'String'));
szdip = size(handles.potentials,2);
if dipnum < 1
    dipnum == 1
elseif dipnum > szdip
    dipnum == szdip;
end
dipnum = round(dipnum);
set(handles.editDipoleNumber, 'String', dipnum);

figure;
utilbem_headplot(handles.potentials(:,dipnum), handles.mesh, eloc);
set(gcf, 'Name', 'Figure: Potential distribution', 'NumberTitle', 'off', 'Color', [0.925 0.957 1]);

% --- Executes on button press in pushbuttonLoadDipoles.
function pushbuttonLoadDipoles_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadDipoles (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[file, path] = uigetfile('*.dip;*.sdip');
if ~isequal(file, 0) && length(file) > 5
    if file(length(file)-3)=='s' %sdip
        handles.dipoles_name = file(1:length(file)-5);
        file = [path file(1:length(file)-5)];
        handles.dipoles.pos = load([file '.sdip'], '-ascii');
        handles.dipoles.sym = 1;
        
    else % dip
        handles.dipoles_name = file(1:length(file)-4);
        file = [path file(1:length(file)-4)];
        handles.dipoles.pos = load([file '.dip'], '-ascii');
        handles.dipoles.sym = 0;
    end
    % Update handles structure
    guidata(hObject, handles);
    set(handles.editNumberofdipoles,'String',size(handles.dipoles.pos,1));
end
%update_display(handles);
%set_model_changed(handles);
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
if (isempty(handles.mesh))
    set(handles.editMeshName,'String',[]);
    set(handles.editNumberofLayers,'String',[]);
    set(handles.editNumberofNodes,'String',[]);
    set(handles.editNumberofElements,'String',[]);
    set(handles.editNumberofNE,'String',[]);
    set(handles.pushbuttonShowMesh,'Enable','off');
else
    set(handles.editMeshName,'String',handles.mesh.name);
    set(handles.editNumberofLayers,'String',handles.mesh.num_boundaries);
    set(handles.editNumberofNodes,'String',handles.mesh.num_nodes);
    set(handles.editNumberofElements,'String',handles.mesh.num_elements);
    set(handles.editNumberofNE,'String',handles.mesh.num_node_elem);
    set(handles.pushbuttonShowMesh,'Enable','on');
end

if (isempty(handles.model))
    set(handles.editModelName,'String',[]);
    %set(handles.editScalpCond,'String',[]);
    %set(handles.editSkullCond,'String',[]);
    %set(handles.editCSFCond,'String',[]);
    %set(handles.editBrainCond,'String',[]);
    set(handles.checkbox1,'Value',1);
    set(handles.modelProgressText,'String','No Model');
    set(handles.pushbuttonCreateModel,'Enable','off');
    handles.model_changed = 0;
    guidata(handles.figure1, handles);
else
    set(handles.editModelName,'String',handles.model.name);
    set(handles.editScalpCond,'String',handles.model.cond(1));
    set(handles.editSkullCond,'String',handles.model.cond(2));
    if length(handles.model.cond) == 3
        set(handles.editCSFCond,'String',[]);
        set(handles.editBrainCond,'String',handles.model.cond(3));
    elseif length(handles.model.cond) == 4
        set(handles.editCSFCond,'String',handles.model.cond(3));
        set(handles.editBrainCond,'String',handles.model.cond(4));
    end
    set(handles.modelProgressText,'String','BEM Model Loaded');
end

if (isempty(handles.session))
    set(handles.editSessionName,'String',[]);
    set(handles.pushbuttonGenerateTM,'Enable','off');
    set(handles.sessionProgressText,'String','No Session');
    set(handles.pushbutton5, 'Enable', 'off');
else
    set(handles.editSessionName,'String',handles.session.name);
    set(handles.editNumberofSensors,'String',length(unique(handles.session.Smatrix(:,1))));
    set(handles.sessionProgressText,'String','Session Loaded');
    if ~isfield(handles.session,'tmte')
        set(handles.pushbuttonGenerateTM,'Enable','on');
    end
    set(handles.pushbutton5, 'Enable', 'on');
end
if ~isfield(handles,'potentials')
    set(handles.editDipoleNumber,'Enable','off');
    set(handles.pushbutton8,'Enable','off')
else
    set(handles.editDipoleNumber,'Enable','on');
    set(handles.pushbutton8,'Enable','on')
end
if (isfield(handles,'sensors'))
    set(handles.uipanelSensLoad, 'visible', 'on')
else
    set(handles.uipanelSensLoad, 'visible', 'off')
end
if isfield(handles.mesh,'num_boundaries')
    if handles.mesh.num_boundaries == 4
        set(handles.uipanelCSF, 'visible', 'on')
    elseif handles.mesh.num_boundaries == 3
        set(handles.uipanelCSF, 'visible', 'off')
    end
end

% --------------------------------------------------------------------
function Load_Mesh_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Mesh_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[file, path] = uigetfile('*.bei');
if ~isequal(file, 0) && length(file) > 5
    % remove .bei extension
    file = [path file(1:length(file)-4)];
    handles.mesh = bem_load_mesh(file);
    handles.model = [];
    handles.session = [];
    % Update handles structure
    guidata(hObject, handles);
end
update_display(handles);


function handles = load_model(hObject,handles, file)
    msave = load(file, '-MAT');
    handles.mesh = bem_load_mesh(msave.mesh_name);
    handles.model = bem_create_model(msave.name, handles.mesh, msave.cond, msave.mod);
    handles.session = [];
    % Update handles structure
    guidata(hObject, handles);  % comment out


% --------------------------------------------------------------------
function Load_Model_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Model_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.model');
if ~isequal(file, 0) && length(file) > 7
    file = [path file];
    handles = load_model(hObject,handles, file);
end
update_display(handles);


% --------------------------------------------------------------------
function Load_Session_Menu_Callback(hObject, eventdata, handles)
% hObject    handle to Load_Session_Menu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.session');
if ~isequal(file, 0) && length(file) > 9
    file = [path file];
    ssave = load(file, '-MAT');
    handles = load_model(hObject,handles, [path ssave.model_name, '.model']);
    handles.session = bem_create_session(ssave.name, handles.model, ssave.Smatrix);
    handles.session = bem_load_transfer_matrix(handles.session, 'tmte');
    handles.sensors = sensMatToCoord(handles.session.Smatrix, handles.mesh.coord);
    % Update handles structure
    guidata(hObject, handles);
end
update_display(handles);


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
%if get(handles.radiobuttonNodeList, 'Value') == 1
%    handles.LoadfromNodeList = 1;
%elseif get(handles.radiobuttonCoordinates, 'Value') == 1
%    handles.LoadfromNodeList = 0;
%end
%update_display(handles);   





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




% --- Executes on button press in radiobuttonNodeList.
function radiobuttonNodeList_Callback(hObject, eventdata, handles)
% hObject    handle to radiobuttonNodeList (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radiobuttonNodeList


% --------------------------------------------------------------------
function uipanel5_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)





function editDipoleNumber_Callback(hObject, eventdata, handles)
% hObject    handle to editDipoleNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editDipoleNumber as text
%        str2double(get(hObject,'String')) returns contents of editDipoleNumber as a double


% --- Executes during object creation, after setting all properties.
function editDipoleNumber_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editDipoleNumber (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on button press in pushbuttonShowSensors.
function pushbuttonShowSensors_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonShowSensors (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isfield(handles,'sensors')
    handles.sensors = sensMatToCoord(handles.session.Smatrix, handles.mesh.coord);
end

coordt = handles.mesh.coord;
ma = mean(coordt);
coordt = coordt - ones(length(coordt),1)*ma; 
%figure;
h = eeglab_plotmesh(handles.mesh.elem, coordt);
set(gcf, 'Name', 'Figure: Sensors', 'NumberTitle', 'off', 'Color', [0.925 0.957 1]);
hold

if size(handles.sensors, 2) == 3
    eloc = handles.sensors;
else
    eloc = handles.mesh.coord(handles.sensors,:);
end
eloc = eloc - ones(length(eloc),1)*ma;
plot3(eloc(:,1),eloc(:,2),eloc(:,3),'ko','LineWidth',2,'MarkerEdgeColor',...
    'k','MarkerFaceColor','r','MarkerSize',6);

