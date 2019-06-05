function varargout = Source_space_generation(varargin)
% SOURCE_SPACE_GENERATION M-file for Source_space_generation.fig
%      SOURCE_SPACE_GENERATION, by itself, creates a new SOURCE_SPACE_GENERATION or raises the existing
%      singleton*.
%
%      H = SOURCE_SPACE_GENERATION returns the handle to a new SOURCE_SPACE_GENERATION or the handle to
%      the existing singleton*.
%
%      SOURCE_SPACE_GENERATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SOURCE_SPACE_GENERATION.M with the given input arguments.
%
%      SOURCE_SPACE_GENERATION('Property','Value',...) creates a new SOURCE_SPACE_GENERATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Source_space_generation_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Source_space_generation_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Source_space_generation

% Last Modified by GUIDE v2.5 13-Jun-2011 10:48:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Source_space_generation_OpeningFcn, ...
                   'gui_OutputFcn',  @Source_space_generation_OutputFcn, ...
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


% --- Executes just before Source_space_generation is made visible.
function Source_space_generation_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Source_space_generation (see VARARGIN)

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


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Source_space_generation wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Source_space_generation_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



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


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sp = str2num(get(handles.edit1,'String'));
th = str2num(get(handles.edit2,'String'));

set(handles.text3,'String','Calculating...'); pause(0.5);

[C1,E1] = mesh_readsmf([handles.MeshFolder '/Brain.smf'],0,0,0,1); % subject's scalp mesh
% C1 and E1 are brain mesh
[rw, ss, dim, inm] = Create_regular_source_space(C1, E1, sp, th);

Ns = size(rw, 1);
so = zeros(3*Ns, 6);
so(1:Ns, 1:3) = rw;
so(1:Ns, 4) = 1;
so(1+Ns:2*Ns, 1:3) = rw;
so(1+Ns:2*Ns, 5) = 1;
so(1+Ns*2:3*Ns, 1:3) = rw;
so(1+Ns*2:3*Ns, 6) = 1;

p = handles.MeshFolder; % save the files in mesh folder
%f = handles.sensorname; % save satrting with the name of sensor file

lof = length(p);
if p(lof) ~= '/'
    p(lof+1) = '/';
end

if isfield(handles,'arg_subject')
    f = [handles.arg_subject];
else
    f = 'temp';
end


% save source space
save([p f '_sourcespace.dip'], 'so', '-ascii'); 

set(handles.text3,'String','Source space saved!');

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MeshFolder = uigetdir;
set(handles.textof, 'String', handles.MeshFolder);
% Update handles structure
guidata(handles.figure1, handles);




% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

sp = str2num(get(handles.edit1,'String'));
th = str2num(get(handles.edit2,'String'));

set(handles.text5,'String','Calculating...'); pause(0.5);

[C1,E1] = mesh_readsmf([handles.MeshFolder '/Brain.smf'],0,0,0,1); % subject's scalp mesh
% C1 and E1 are brain mesh
[so, ss] = Create_symmetric_source_space(C1, E1, sp, th);

p = handles.MeshFolder; % save the files in mesh folder
%f = handles.sensorname; % save satrting with the name of sensor file

lof = length(p);
if p(lof) ~= '/'
    p(lof+1) = '/';
end

if isfield(handles,'arg_subject')
    f = [handles.arg_subject];
else
    f = 'temp';
end


% save source space
save([p f '_sourcespace.sdip'], 'so', '-ascii'); 

set(handles.text5,'String','Source space saved!');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [so, ss, dim, inm] = Create_regular_source_space(Coord, Elem, spacing, thr);

ma = max(Coord(:,2:4));
mi = min(Coord(:,2:4));
md = round((ma - mi)/spacing);

n=0;
for i=1:md(1)+1
    a(i) = mi(1)+spacing*(i-1)-spacing/2;
    for j=1:md(2)+1
        b(j) = mi(2)+spacing*(j-1)-spacing/2;
        for k=1:md(3)+1
            c(k) = mi(3)+spacing*(k-1)-spacing/2;
            n=n+1;
            so(n,1:3)=[a(i) b(j) c(k)];
        end
    end
end

[dim, inm] = utilmesh_check_source_space(so, Coord, Elem);

k = find(inm == 1);   % dipoles inside the mesh
l = find(dim < thr);  % dipoles closer to the mesh less than thr
m = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr
so = so(m,:);

ne=size(so,1);    
ss=zeros(ne*3,7);
ss(1:ne,2:4)=so;
ss(ne+1:2*ne,2:4)=so;
ss(2*ne+1:3*ne,2:4)=so;
ss(1:ne,5)=1;
ss(ne+1:2*ne,6)=1;
ss(2*ne+1:3*ne,7)=1;
ss(:,1)=[1:3*ne]';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [so, ss] = Create_symmetric_source_space(Coord, Elem, spacing, thr);

ma = max(Coord(:,2:4));
mi = min(Coord(:,2:4));
md = round((ma - mi)/spacing);

n=0;
for i=1:md(1)+1
    a(i) = mi(1)+spacing*(i-1)-spacing/2;
    for j=1:md(2)+1
        b(j) = mi(2)+spacing*(j-1)-spacing/2;
        for k=1:md(3)+1
            c(k) = mi(3)+spacing*(k-1)-spacing/2;
            n=n+1;
            ss(n,1:3)=[a(i) b(j) c(k)];
        end
    end
end

% symmetry according to y axis
smy = (max(Coord(:,3))+min(Coord(:,3)))/2;

ksmy = find(ss(:,2) > smy);
ss1 = ss(ksmy,:);
ss2 = ss1;
ss2(:,2) = -ss2(:,2);
ss2(:,2) = ss2(:,2)+smy*2;

[dim1, inm1] = utilmesh_check_source_space(ss1, Coord, Elem);
[dim2, inm2] = utilmesh_check_source_space(ss2, Coord, Elem);


k = find(inm1 == 1);   % dipoles inside the mesh
l = find(dim1 < thr);  % dipoles closer to the mesh less than thr
m1 = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr

k = find(inm2 == 1);   % dipoles inside the mesh
l = find(dim2 < thr);  % dipoles closer to the mesh less than thr
m2 = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr

mi = intersect(m1,m2);

ss1 = ss1(mi,:);
ss2 = ss2(mi,:);

% first Ns y directed dipoles
Ns = size(ss1, 1);
so1 = zeros(3*Ns, 6);
so1(1:Ns, 1:3) = ss1;
so1(1:Ns, 4) = 1;
so1(1+Ns:2*Ns, 1:3) = ss1;
so1(1+Ns:2*Ns, 5) = 1;
so1(1+Ns*2:3*Ns, 1:3) = ss1;
so1(1+Ns*2:3*Ns, 6) = 1;

% second Ns dipoes are the symmetric ones
so2 = zeros(3*Ns, 6);
so2(1:Ns, 1:3) = ss2;
so2(1:Ns, 4) = 1;
so2(1+Ns:2*Ns, 1:3) = ss2;
so2(1+Ns:2*Ns, 5) = -1;
so2(1+Ns*2:3*Ns, 1:3) = ss2;
so2(1+Ns*2:3*Ns, 6) = 1;

so = [so1;so2]; 
