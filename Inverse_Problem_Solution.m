function varargout = Inverse_Problem_Solution(varargin)
% INVERSE_PROBLEM_SOLUTION M-file for Inverse_Problem_Solution.fig
%      INVERSE_PROBLEM_SOLUTION, by itself, creates a new INVERSE_PROBLEM_SOLUTION or raises the existing
%      singleton*.
%
%      H = INVERSE_PROBLEM_SOLUTION returns the handle to a new INVERSE_PROBLEM_SOLUTION or the handle to
%      the existing singleton*.
%
%      INVERSE_PROBLEM_SOLUTION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INVERSE_PROBLEM_SOLUTION.M with the given input arguments.
%
%      INVERSE_PROBLEM_SOLUTION('Property','Value',...) creates a new INVERSE_PROBLEM_SOLUTION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Inverse_Problem_Solution_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Inverse_Problem_Solution_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Inverse_Problem_Solution

% Last Modified by GUIDE v2.5 07-Jun-2011 16:54:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Inverse_Problem_Solution_OpeningFcn, ...
                   'gui_OutputFcn',  @Inverse_Problem_Solution_OutputFcn, ...
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


% --- Executes just before Inverse_Problem_Solution is made visible.
function Inverse_Problem_Solution_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Inverse_Problem_Solution (see VARARGIN)

% Parse arguments and set handles as necessary
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'EEGstruct')
        i = i + 1;
        handles.EEG = varargin{i};
    elseif strcmp(varargin{i}, 'subjectdir')
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


% UIWAIT makes Inverse_Problem_Solution wait for user response (see UIRESUME)
% uiwait(handles.figure1);

sensor_file = [handles.arg_subject '_' handles.arg_session '.sensors'];
a = dir(sensor_file);
if size(a,1) > 0
    se = load(sensor_file,'-mat');
    handles.eloc = se.eloc;
    set(handles.text2, 'String',se.fn);
end

a = dir(['ori_sen_loc.mat']);
if size(a,1) > 0
    load ori_sen_loc
    handles.elocfn = sens_fn;
    handles.eloc = readlocs(sens_fn);
    set(handles.text2, 'String', handles.elocfn);
end

% Choose default command line output for Inverse_Problem_Solution
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = Inverse_Problem_Solution_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function edit_compind_Callback(hObject, eventdata, handles)
% hObject    handle to edit_compind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_compind as text
%        str2double(get(hObject,'String')) returns contents of edit_compind as a double


% --- Executes during object creation, after setting all properties.
function edit_compind_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_compind (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton_dipfit.
function pushbutton_dipfit_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_dipfit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%if get(handles.checkbox1, 'Value') == 1
%    a = dir([handles.arg_subject '_' handles.arg_session '_warping.mat']);
%    if size(a,1) > 0
%        handles.mri_based = 0;
%    else
%        handles.mri_based = 1;
%    end

    %  mri-based realistic
%    if handles.mri_based == 1
%        constr.reduce = [1 2 3];
%        constr.expand = [1 2 3 1 2 3];
%        constr.mirror = [1 1 1 1 -1 1];
%    else
%        constr.reduce = [1 2 3];
%        constr.expand = [1 2 3 1 2 3];
%        constr.mirror = [1 1 1 1 -1 1];
%    end
%else
%    constr = [];
%end
constr = [];

comp_index = str2num(get(handles.edit_compind,'String'));

v = evalin('base','EEG');
if isfield(v.etc,'nft')
    dip1 = v.etc.nft.model;
else
    dip1=[];
end

% check if mr-based realistic or warped mni
%a = dir(['Scalp.smf']);
a = dir([handles.arg_subject '_' handles.arg_session '_warping.mat']);
if size(a,1) == 0
    handles.mri_based = 1;
    [dip, session] = ip_dipolefitting(handles.EEG, handles.eloc, handles.arg_subject, handles.arg_session, comp_index, constr,[]);
else
    handles.mri_based = 0;
    fw = [handles.arg_subject '_' handles.arg_session '_warping'];
    load(fw)
    [dip, session] = ip_dipolefitting(handles.EEG, handles.eloc, handles.arg_subject, handles.arg_session, comp_index, constr, warping_param.back);
end

%save Dipole_soln dip

handles.EEG.etc.nft.session = session;

if ~isempty(dip1)
    for i=1:length(comp_index)
        dip1(comp_index(i)) = dip(comp_index(i));
    end
else
    dip1=dip;
end

% Need two steps to set a structure field
assignin('base','NFT_temp', dip1);
evalin('base', 'EEG.etc.nft.model=NFT_temp; clear NFT_temp;');

handles.dipoles_str = dip1;

guidata(handles.figure1, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.*');
if ~isequal(file, 0) && length(file) > 1
    handles.eloc = readlocs([path file]);   % subject's electrode locations
    sens_fn = [path file];
    handles.elocfn = [path file];
    handles.sensorpath = path;
    handles.sensorname = file(1:length(file)-4);
end
set(handles.text2, 'String', handles.elocfn);
% Update handles structure
guidata(handles.figure1, handles);



% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
a = dir([handles.arg_subject '_' handles.arg_session '_warping.mat']);
if size(a,1) > 0
    handles.mri_based = 0;
else
    handles.mri_based = 1;
end

% mri-based realistic
if handles.mri_based == 1
    mri_file = [handles.arg_subject '_mri'];
else
    % warped mni
    eeglab_folder = dirname(which('eeglab'));
    mri_file = [eeglab_folder '/plugins/dipfit2.2/standard_BEM/standard_mri.mat'];
end

if ~isfield(handles,'dipoles_str')
    v = evalin('base','EEG');
    if isfield(v.etc,'nft')
        handles.dipoles_str = v.etc.nft.model;
    else
        error('Please run dipole fitting!')
    end
end

dip = handles.dipoles_str;
for i = 1:size(dip,2);    dip(i).momxyz = dip(i).momxyz(:)'; end
eeglab_dipplot(dip,'mri',mri_file,'projimg', 'off', 'projlines', 'off', 'axistight', 'on', 'cornermri','on','normlen','on');


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1



% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = handles.arg_subject;
ses_name = handles.arg_session;
if isfield(handles,'MeshFolder')
    Forward_Problem_Solution('subjectdir', handles.MeshFolder, 'subject', subj_name, 'session', ses_name);
else
    Forward_Problem_Solution('subject', subj_name, 'session', ses_name);
end



% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
subj_name = handles.arg_subject;
ses_name = handles.arg_session;
if isfield(handles,'MeshFolder')
    FP_FEM('subjectdir', handles.MeshFolder, 'subject', subj_name, 'session', ses_name)
else
    FP_FEM('subject', subj_name, 'session', ses_name)
end
