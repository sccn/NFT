function varargout = Mesh_generation(varargin)
% MESH_GENERATION M-file for Mesh_generation.fig
%      MESH_GENERATION, by itself, creates a new MESH_GENERATION or raises the existing
%      singleton*.
%
%      H = MESH_GENERATION returns the handle to a new MESH_GENERATION or the handle to
%      the existing singleton*.
%
%      MESH_GENERATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MESH_GENERATION.M with the given input arguments.
%
%      MESH_GENERATION('Property','Value',...) creates a new MESH_GENERATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Mesh_generation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Mesh_generation_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Mesh_generation

% Last Modified by GUIDE v2.5 27-Jun-2011 08:40:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Mesh_generation_OpeningFcn, ...
                   'gui_OutputFcn',  @Mesh_generation_OutputFcn, ...
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


% --- Executes just before Mesh_generation is made visible.
function Mesh_generation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Mesh_generation (see VARARGIN)

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
    set(handles.textof, 'String', handles.MeshFolder);
end
if isfield(handles,'arg_subject')
    set(handles.editMeshName,'String',handles.arg_subject);
end
if isfield(handles, 'MeshFolder') & isfield(handles,'arg_subject')
    handles.data.filename = handles.arg_subject;
    of = handles.MeshFolder;
    lof = length(of);
    if of(lof) ~= filesep
        of(lof+1) = filesep;
    end
    handles.data.filepath = of;

    file = [of handles.data.filename '_segments'];
    handles.param.segmentation = file;

    load(file)
    handles.Sca = Segm.scalpmask;
    handles.Bra = Segm.brainmask;
    handles.Csf = Segm.innerskullmask;
    handles.Skull = Segm.outerskullmask;
    set(handles.textVolumes,'String',file);
end
if isfield(handles,'MeshFolder') & isfield(handles,'Sca')
    set(handles.pushbuttonMesh, 'Enable', 'on')
end

% Update handles structure
%guidata(handles.figure1, handles);

    
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Mesh_generation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Mesh_generation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --------------------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
% hObject    handle to Untitled_1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function Load_Callback(hObject, eventdata, handles)
% hObject    handle to Load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.mat');
if ~isequal(file, 0) && length(file) > 5
    handles.data.filename = file(1:length(file)-4);
    handles.data.filepath = path;
    % remove .hdr extension
    file = [path file];
    load(file)
    handles.Sca = Segm.scalpmask;
    handles.Bra = Segm.brainmask;
    handles.Csf = Segm.innerskullmask;
    handles.Skull = Segm.outerskullmask;
    handles.param.segmentation = [path file];
else
    error('Error loading file...')
end

% Update handles structure
guidata(handles.figure1, handles);



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


% --- Executes on button press in pushbuttonMesh.
function pushbuttonMesh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%tic
mesh_name = get(handles.editMeshName,'String');
if isempty(mesh_name)
    error('You must enter a mesh name');
end
of = handles.MeshFolder; % Output Folder

handles.param.mesh_folder = of;

if isempty(of)
    error('You must enter an output folder...')
end

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
    handles.MeshFolder = of;
end

    transform = size(handles.Sca) / 2;
    save([of 'transform'], 'transform', '-ascii'); 


set(handles.textStatus,'String','Saving volumes in raw...'); pause(0.5)
Mesh_writeraw(handles.Sca, [of 'Scalp']);
Mesh_writeraw(handles.Bra, [of 'Brain']);
Mesh_writeraw(handles.Csf, [of 'Csf']);
Mesh_writeraw(handles.Skull, [of 'Skull']);
set(handles.textStatus,'String','Volumes saved...'); pause(0.5)
[K,L,M] = size(handles.Sca);

% load mesh configuration for path names
conf = nft_get_config;

tis_type = cell(1,4);
tis_type{1} = 'Scalp';
tis_type{2} = 'Brain';
tis_type{3} = 'Csf';
tis_type{4} = 'Skull';

set(handles.textStatus,'String','Triangulating volumes...'); pause(0.5)
hh = waitbar(0,'Triangulating volumes...');
for tis = 1:4
    waitbar(tis/4);
    tt = char(tis_type(tis));
    set(handles.textStatus,'String',['Triangulating ' tt ' volume...']); pause(0.5)
    a= sprintf('"%s" -t 10 -dr1 "%s%s.raw" %d %d %d -f "%s%s.asc" -ot', conf.asc, of,tt, K, L, M,of,tt);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
end
close(hh);

% generate a file for StepSc.txt for coarsening and smoothing
f=fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'correct 5\n');
fprintf(f, 'smooth 1\n');
%fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
%fprintf(f, 'split intersect\n'); % XXX yeni
%fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

% generate a file for StepSc2.txt for final improvement
f=fopen(sprintf('%sStepSc2.txt',of), 'w');
fprintf(f, 'correct 2\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'fill holes\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'fill holes\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

NumberNodes = round(str2num(get(handles.editNumNodes, 'String')));
if get(handles.radiobuttonLinear, 'Value') == 1
    Quad = 0;
else
    Quad = 1;
end

csi = 300000; i = 1;
while NumberNodes < csi
    i = i+1;
    csi(i) = round(csi(i-1) / 1.5);
end
csi(i) = NumberNodes;

nsteps = length(csi); % number of coarsening steps
ntis = 4;  % number of tissues (4)

%%%%%%%%%%%%%%%%%% Coarsening and correcting...
hh = waitbar(0,'Coarsening and correcting...');
for tis = 1:ntis
    tt = char(tis_type(tis));
    set(handles.textStatus,'String',['Coarsening and correcting ' tt ' surface...']); pause(0.5)
    
    a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.asc"',conf.showmesh,of,of,tt);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
    
     for iter = 1:nsteps
         waitbar(((tis-1)*nsteps + iter) / (nsteps*ntis));
         a = sprintf('"%s" -c 0.5 -m 5 -o "%s%s.smf" -t %d "%sScS.smf"', conf.qslim, of, tt, csi(iter), of);
         [status, result] = system(a);
         if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
 
         a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
         [status, result] = system(a);
         if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
         if Quad & iter == 6
             copyfile([of 'ScS.smf'], [of tt 'f.smf'])
         end

     end
     a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
     [status, result] = system(a);
     if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
     
     movefile([of 'ScS.smf'], [of tt '.smf']); % XXX yeni
end
close(hh)


nl = str2num(get(handles.edit_nl, 'String'));
handles.param.no_layers = nl;
handles.param.linear = 1;
set(handles.textStatus,'String','Final mesh correction...'); pause(0.5)

% Final correction with Matlab functions
mesh_final_correction(of, nl);

%if Quad == 0
    % Local mesh refinement
    if (get(handles.checkboxLMR, 'Value') == 1)
        ratio_lmr = str2num(get(handles.editRatioLMR, 'String'));
        set(handles.textStatus,'String','Local mesh refinement...'); pause(0.5)
        mesh_local_refinement(of, nl, ratio_lmr);
        handles.param.lmr = ratio_lmr;
    end
%end

set(handles.textStatus,'String','Generating quadratic meshes...'); pause(0.5)

% generate quadratic meshes
if Quad == 1
    handles.param.linear = 0;
    for tis = 1:ntis
        tt = char(tis_type(tis));
        a = sprintf('"%s" "%s%s.smf" "%s%sf.smf" %s%sq', conf.quad, of,tt,of,tt,of,tt);
        [status, result] = system(a);
        if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
    end
end
       

% delete unnecessary files (.raw, .asc, Scs.smf, and StepSc.txt)
%delete([of 'ScS.smf']);
delete([of 'StepSc.txt']);
delete([of 'StepSc2.txt']);
for tis = 1:ntis
    tt = char(tis_type(tis));
    delete([of tt '.raw']);
    delete([of tt '.asc']);
end


set(handles.textStatus,'String','Mesh Generation done! Saving...'); pause(0.5)
mesh_name = get(handles.editMeshName,'String');
mesh_read_write(of, mesh_name, nl, Quad); % generate nl-layer head model
handles.param.mesh_name = mesh_name;
parameters = handles.param;
save([of mesh_name '_mesh'], '-STRUCT', 'parameters')

set(handles.textStatus,'String','Mesh saved!'); pause(0.5)


%toc
% --- Executes on button press in pushbuttonFolder.
function pushbuttonFolder_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
%global NFM_MESH_FOLDER

handles.MeshFolder = uigetdir;
%NFM_MESH_FOLDER = handles.MeshFolder;


set(handles.textof, 'String', handles.MeshFolder);
% Update handles structure
guidata(handles.figure1, handles);

if isfield(handles,'Sca')
    if ~isempty(handles.MeshFolder) & ~isempty(handles.Sca)
        set(handles.pushbuttonMesh, 'Enable', 'on')
    end
end



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




% --- Executes on button press in pushbuttonLoadVol.
function pushbuttonLoadVol_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonLoadVol (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.mat');
if ~isequal(file, 0) && length(file) > 5
    handles.data.filename = file(1:length(file)-4);
    handles.data.filepath = path;
    % remove .hdr extension
    file = [path file];
    load(file)
    handles.Sca = Segm.scalpmask;
    handles.Bra = Segm.brainmask;
    handles.Csf = Segm.innerskullmask;
    handles.Skull = Segm.outerskullmask;
    set(handles.textVolumes,'String',file);
else
    error('Error loading file...')
end

if isfield(handles,'MeshFolder')
    if ~isempty(handles.MeshFolder) & ~isempty(handles.Sca)
        set(handles.pushbuttonMesh, 'Enable', 'on')
    end
end

% Update handles structure
guidata(handles.figure1, handles);




% --- Executes on button press in checkboxLMR.
function checkboxLMR_Callback(hObject, eventdata, handles)
% hObject    handle to checkboxLMR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkboxLMR



function editRatioLMR_Callback(hObject, eventdata, handles)
% hObject    handle to editRatioLMR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editRatioLMR as text
%        str2double(get(hObject,'String')) returns contents of editRatioLMR as a double


% --- Executes during object creation, after setting all properties.
function editRatioLMR_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editRatioLMR (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





function edit_nl_Callback(hObject, eventdata, handles)
% hObject    handle to edit_nl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_nl as text
%        str2double(get(hObject,'String')) returns contents of edit_nl as a double


% --- Executes during object creation, after setting all properties.
function edit_nl_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_nl (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function Mesh_writeraw(file, fn);
A = single(file);
norm = max(max(max(A)));
A = A * 20 / norm;
filename = [fn '.raw'];
f=fopen(filename, 'w+');
fwrite(f, A, 'uint8');
fclose(f);
clear A norm f 



function editNumNodes_Callback(hObject, eventdata, handles)
% hObject    handle to editNumNodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editNumNodes as text
%        str2double(get(hObject,'String')) returns contents of editNumNodes as a double


% --- Executes during object creation, after setting all properties.
function editNumNodes_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editNumNodes (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbuttonFEM.
function pushbuttonFEM_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFEM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




% load mesh configuration for path names
conf = nft_get_config;

of = handles.MeshFolder; % Output Folder
lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end
mesh_name = get(handles.editMeshName,'String');

fn = [of mesh_name '.bei'];

if exist(fn) == 0
    error('BEM mesh not found!')
end

set(handles.textStatus,'String','Generating FEM mesh...'); pause(0.5)

mesh = bem_load_mesh([of mesh_name]);
R = mesh_find_regions(mesh);
Coord(:,2:4) = mesh.coord;
Elem(:,2:4) = mesh.elem;
Coord(:,1) = [1:length(Coord)]';
Elem(:,1) = [1:length(Elem)]';

fn = [of mesh_name '.smesh'];
WriteSMESH(fn, Coord, Elem, R);

% call tetgen
a = sprintf('"%s" -pq1.4a5A "%s"', conf.tetgen, fn);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

% convert into metu-fem mesh
nl = str2num(get(handles.edit_nl, 'String'));
if nl == 4
    cnd_str = '1=C1 2=C2 3=C3 4=C4';
elseif nl == 3
    cnd_str = '1=C1 2=C2 3=C3';
end
fn = [of mesh_name '.1'];
a = sprintf('"%s" "%s" %s', conf.tetgen2msh, fn, cnd_str);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

set(handles.textStatus,'String','Linear FEM mesh generated!'); pause(0.5)



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% load mesh configuration for path names
conf = nft_get_config;

of = handles.MeshFolder; % Output Folder
lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end
mesh_name = get(handles.editMeshName,'String');

fn = [of mesh_name '.bei'];

if exist(fn) == 0
    error('BEM mesh not found!')
end

set(handles.textStatus,'String','Generating FEM mesh...'); pause(0.5)

mesh = bem_load_mesh([of mesh_name]);
R = mesh_find_regions(mesh);
Coord(:,2:4) = mesh.coord;
Elem(:,2:4) = mesh.elem;
Coord(:,1) = [1:length(Coord)]';
Elem(:,1) = [1:length(Elem)]';

fn = [of mesh_name '.smesh'];
WriteSMESH(fn, Coord, Elem, R);

% call tetgen - generate linear tetrahedral mesh
% a = sprintf('"%s" -pq1.4a12A "%s"', conf.tetgen, fn);
a = sprintf('"%s" -pq1.4a120A "%s"', conf.tetgen, fn);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

% convert into metu-fem mesh
nl = str2num(get(handles.edit_nl, 'String'));
if nl == 4
    cnd_str = '1=C1 2=C2 3=C3 4=C4';
elseif nl == 3
    cnd_str = '1=C1 2=C2 3=C3';
end
fn = [of mesh_name '.1'];
a = sprintf('"%s" "%s" %s', conf.tetgen2msh, fn, cnd_str);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

% convert into quadratic mesh
fn_q = [of mesh_name '.q'];
a = sprintf('"%s" -o "%s".msh "%s".msh', conf.lin2quad, fn_q, fn);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

set(handles.textStatus,'String','Quadratic FEM mesh generated!'); pause(0.5)



function WriteSMESH(name,Coord,Elem,Regions)
% Saves the mesh in .SMESH format for use with Tetgen
% The format is described in:
% http://tetgen.berlios.de/fformats.html
% The .smesh format is slightly simpler than more general .poly
% The Regions parameter must specify a point inside each region
% TetGen will mark output tetrahedra using these region markers

nnp=size(Coord,1);
nel=size(Elem,1);
if ~isempty(Regions)
	nreg = size(Regions,1);
	Reg = zeros(nreg, 6);
	Reg(:,2:4) = Regions;
	Reg(:,1) = 1:nreg;
	Reg(:,5) = 1:nreg;
	Reg(:,6) = -1;
else
	nreg = 0;
	Reg = [];
end

% make sure Node indices are correct
Coord(:,1) = 1:nnp;

fid=fopen(name, 'w');

fprintf(fid, '# Part 1 - node list\n');
fprintf(fid, '# node count, 3 dim, no attr, no boundary \n');
fprintf(fid, '%d 3 0 0\n', nnp);
fprintf(fid, '# Node index, node coordinates\n');
fprintf(fid, '%f %f %f %f\n',Coord');

fprintf(fid, '# Part 2 - facet list\n');
fprintf(fid, '# facet count, no boundary marker\n');
fprintf(fid, '%d\n',nel);
fprintf(fid, '# facets\n');
fprintf(fid, '3 %d %d %d\n', Elem(:,2:4)');

fprintf(fid, '# Part 3 - hole list\n');
fprintf(fid, '0            # no hole\n');

fprintf(fid, '# Part 4 - region list\n');
fprintf(fid, '%d\n', nreg);

if nreg > 0
	fprintf(fid, '%d %f %f %f %d %d\n', Reg');
end

fclose(fid);
