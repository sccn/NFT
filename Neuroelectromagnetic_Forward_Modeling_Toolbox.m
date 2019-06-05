function varargout = Neuroelectromagnetic_Forward_Modeling_Toolbox(varargin)
% NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX M-file for Neuroelectromagnetic_Forward_Modeling_Toolbox.fig
%      NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX, by itself, creates a new NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX or raises the existing
%      singleton*.
%
%      H = NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX returns the handle to a new NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX or the handle to
%      the existing singleton*.
%
%      NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX.M with the given input arguments.
%
%      NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX('Property','Value',...) creates a new NEUROELECTROMAGNETIC_FORWARD_MODELING_TOOLBOX or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Neuroelectromagnetic_Forward_Modeling_Toolbox_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Neuroelectromagnetic_Forward_Modeling_Toolbox_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Neuroelectromagnetic_Forward_Modeling_Toolbox

% Last Modified by GUIDE v2.5 06-Jan-2011 14:33:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Neuroelectromagnetic_Forward_Modeling_Toolbox_OpeningFcn, ...
                   'gui_OutputFcn',  @Neuroelectromagnetic_Forward_Modeling_Toolbox_OutputFcn, ...
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


% --- Executes just before Neuroelectromagnetic_Forward_Modeling_Toolbox is made visible.
function Neuroelectromagnetic_Forward_Modeling_Toolbox_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Neuroelectromagnetic_Forward_Modeling_Toolbox (see VARARGIN)

% Parse arguments and set handles as necessary
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'EEGstruct')
        i = i + 1;
        handles.EEG = varargin{i};
    end
end

% Choose default command line output for Neuroelectromagnetic_Forward_Modeling_Toolbox
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Neuroelectromagnetic_Forward_Modeling_Toolbox wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = Neuroelectromagnetic_Forward_Modeling_Toolbox_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in pushbuttonSegm.
function pushbuttonSegm_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonSegm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    Segmentation('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Segmentation('subject', subj_name, 'session', ses_name);
end


% --- Executes on button press in pushbuttonMesh.
function pushbuttonMesh_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonMesh (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    Mesh_generation('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Mesh_generation('subject', subj_name, 'session', ses_name);
end

% --- Executes on button press in pushbuttonRegis.
function pushbuttonRegis_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonRegis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    handles.elocfn = Coregistration('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    handles.elocfn = Coregistration('subject', subj_name, 'session', ses_name);
end



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    Warping_mesh('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Warping_mesh('subject', subj_name, 'session', ses_name);
end


% --- Executes on button press in pushbuttonFP.
function pushbuttonFP_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    Forward_Problem_Solution('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Forward_Problem_Solution('subject', subj_name, 'session', ses_name);
end




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    Source_space_generation('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Source_space_generation('subject', subj_name, 'session', ses_name);
end

% --- Executes on button press in pushbutton_ip.
function pushbutton_ip_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_ip (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');

if isfield(handles, 'EEG')==0
    EEG = pop_loadset;
    handles.EEG = EEG;
end

if isfield(handles,'SubjectFolder')
    Inverse_Problem_Solution('EEGstruct',handles.EEG,'subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Inverse_Problem_Solution('EEGstruct',handles.EEG,'subject', subj_name, 'session', ses_name);
end

% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.SubjectFolder = uigetdir;
set(handles.textToolboxDir, 'String', handles.SubjectFolder);
% Update handles structure
guidata(handles.figure1, handles);



function editSubjName_Callback(hObject, eventdata, handles)
% hObject    handle to editSubjName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSubjName as text
%        str2double(get(hObject,'String')) returns contents of editSubjName as a double


% --- Executes during object creation, after setting all properties.
function editSubjName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSubjName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function editSesName_Callback(hObject, eventdata, handles)
% hObject    handle to editSesName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of editSesName as text
%        str2double(get(hObject,'String')) returns contents of editSesName as a double


% --- Executes during object creation, after setting all properties.
function editSesName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to editSesName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');
if isfield(handles,'SubjectFolder')
    FP_FEM('subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name)
else
    FP_FEM('subject', subj_name, 'session', ses_name)
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

subj_name = get(handles.editSubjName,'String');
ses_name = get(handles.editSesName,'String');

if isfield(handles, 'EEG')==0
    EEG = pop_loadset;
    handles.EEG = EEG;
end

if isfield(handles,'SubjectFolder')
    Distributed_Source_Localization('EEGstruct',handles.EEG,'subjectdir', handles.SubjectFolder, 'subject', subj_name, 'session', ses_name);
else
    Distributed_Source_Localization('EEGstruct',handles.EEG,'subject', subj_name, 'session', ses_name);
end

