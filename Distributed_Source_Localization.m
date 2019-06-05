function varargout = Distributed_Source_Localization(varargin)
% DISTRIBUTED_SOURCE_LOCALIZATION MATLAB code for Distributed_Source_Localization.fig
%      DISTRIBUTED_SOURCE_LOCALIZATION, by itself, creates a new DISTRIBUTED_SOURCE_LOCALIZATION or raises the existing
%      singleton*.
%
%      H = DISTRIBUTED_SOURCE_LOCALIZATION returns the handle to a new DISTRIBUTED_SOURCE_LOCALIZATION or the handle to
%      the existing singleton*.
%
%      DISTRIBUTED_SOURCE_LOCALIZATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISTRIBUTED_SOURCE_LOCALIZATION.M with the given input arguments.
%
%      DISTRIBUTED_SOURCE_LOCALIZATION('Property','Value',...) creates a new DISTRIBUTED_SOURCE_LOCALIZATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Distributed_Source_Localization_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Distributed_Source_Localization_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
% Author: Zeynep Akalin Acar, SCCN, 2016

% Copyright (C) 2017 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
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


% Edit the above text to modify the response to help Distributed_Source_Localization

% Last Modified by GUIDE v2.5 17-Nov-2016 12:42:20

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Distributed_Source_Localization_OpeningFcn, ...
                   'gui_OutputFcn',  @Distributed_Source_Localization_OutputFcn, ...
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


% --- Executes just before Distributed_Source_Localization is made visible.
function Distributed_Source_Localization_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Distributed_Source_Localization (see VARARGIN)

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

mri_file = [handles.arg_subject '_segments'];
a = dir(mri_file);
if size(a,1) > 0
    load(mri_file,'-mat');
    handles.mri = [Segm.parameters.MRpath Segm.parameters.MRfile]
    set(handles.text4, 'String',handles.mri);
end


% Choose default command line output for Distributed_Source_Localization
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes Distributed_Source_Localization wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Distributed_Source_Localization_OutputFcn(hObject, eventdata, handles) 
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


% --- Executes on selection change in popupmenu3.
function popupmenu3_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu3 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu3

% Source localization

comp_index = str2num(get(handles.edit1,'String'));


% --- Executes during object creation, after setting all properties.
function popupmenu3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Load MRI

[file, path] = uigetfile('*.img');
if ~isequal(file, 0) && length(file) > 5
    handles.parameters.MRfile = file;
    handles.parameters.MRpath = path;
    handles.mri = [path file];
    set(handles.text4, 'String',handles.mri);
end
guidata(handles.figure1, handles);

% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Start Freesurfer

set(handles.text5, 'String','Running Freesurfer for cortical segmentation...'); pause(0.5)
%set(handles.text5, 'String','Running Freesurfer completed!'); pause(0.5)

disp('Running Freesurfer...')
a = sprintf('recon-all -subject FS -sd "%s" -i "%s" -all', handles.MeshFolder, handles.mri);
[status, result] = system(a);
if status ~= 0; error('FreeSurfer:system','Failed to execute: %s',result); end
set(handles.text5, 'String','Running Freesurfer completed!'); pause(0.5)
disp('Freesurfer completed!')


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


NumberNodes = round(str2num(get(handles.edit3, 'String')));

of = handles.MeshFolder;
ofFS = [handles.MeshFolder '/FS/surf/'];
cd(ofFS)
[vl, fl] = freesurfer_read_surf('lh.pial');
[vr, fr] = freesurfer_read_surf('rh.pial');

cd(of)
sname = [handles.arg_subject];

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
    handles.MeshFolder = of;
end

disp('Co-registering the FS mesh with the NFT mesh...')
vertices = [vr; vl];
faces = [fr; fl+max(max(fr))];
[nf, nv] = reducepatch(faces, vertices, NumberNodes*2);
[d, nv] = warping_distafterwarping([0 0 0 0 0 0], nv, nv);
nv = nv + 128;

[Cb, Eb] = ReadSMF('Brain.smf', 0, 0, 0, 1);
plotmesh(Eb(:,2:4), Cb(:,2:4)); view(0,90); hold
plot3(nv(:,1),nv(:,2),nv(:,3),'r.'); pause(1);

ma = mean(nv);
nv1 = nv - ones(length(nv),1) * ma;
rat = 0.95;
nv1 = rat * nv1;
nv1 = nv1 + ones(length(nv1),1) * ma;

nvt = [nv; nv1]; nft = [nf; nf+max(max(nf))];
Vfs = mesh2vol3(nvt, nft);
load([sname '_segments'])
Segm = new_segm(Segm, Vfs);
nsname = [sname 'FS'];
save([nsname '_segments.mat'],'Segm');

disp('Generating a new mesh...')
nft_mesh_generation(nsname, of, 4, 'Segm', Segm, 'lin_femmesh', 1)

ma = mean(nv);
F = nv - ones(length(nv),1) * ma;
rat = 0.95;
F = rat * F;
F = F + ones(length(F),1) * ma;

[Coord,Elem] = ReadSMF('Brain.smf',0,0,0,1); % Brain model
disp('Correcting the source space...');
[so2, k1,k2] = correct_source_space(F, Coord, Elem);

% load mesh configuration for path names
conf = nft_get_config;

%save ss_cor k2 k1 so2
disp('Running procmesh for correction...')
WriteSMF2('FSss.smf', so2, nf); % run showmesh once and save XXX

f=fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

a = sprintf('"%s" -c "%sStepSc.txt" FSss.smf', conf.showmesh2, of);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

[Css, Ess] = ReadSMF('ScS.smf',0,0,0,1); % FS sourcespace

save FSss_cor Css Ess   % sourcespace
disp('Calculating the area of nodes...');
[Ae, An] = area_of_nodes(Css,Ess);
save Node_area Ae An

% find node normals for Css Ess mesh
disp('Calculating node normals...')
Nn = NodeNormals(Css,Ess,1);
ss = [Css(:,2:4) Nn];

save([nsname '_ss.dip'],'ss','-ascii');


disp('Checking scalp sensor locations...')
sensor_file = [handles.arg_subject '_' handles.arg_session '.sensors'];
a = dir(sensor_file);
if size(a,1) > 0
    se = load(sensor_file,'-mat');
    handles.eloc = se.eloc;
    % check sensor locations
    [Cs,Es] = ReadSMF('Scalp.smf',0,0,0,1); % Brain model
    [F2, dmi] = warping_distmeshafterwarping([0 0 0 0 0 0], se.pnt, Cs, Es);
    se.pnt = F2;
    save(sensor_file, '-STRUCT', 'se')
else
    disp('Please co-register electrode locations with the head model')
end

disp('Cortical source space is saved!')

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


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1



% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
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
sname = [handles.arg_subject];
subj_name = [sname 'FS'];
ses_name = [handles.arg_session];
if isfield(handles,'MeshFolder')
    Forward_Problem_Solution('subjectdir', handles.MeshFolder, 'subject', subj_name, 'session', ses_name);
else
    Forward_Problem_Solution('subject', subj_name, 'session', ses_name);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
sname = [handles.arg_subject];
subj_name = [sname 'FS'];
ses_name = [handles.arg_session];
if isfield(handles,'MeshFolder')
    FP_FEM('subjectdir', handles.MeshFolder, 'subject', subj_name, 'session', ses_name)
else
    FP_FEM('subject', subj_name, 'session', ses_name)
end


function V = mesh2vol3(C, E)

vsize = floor(max(C)) + 10;
V = zeros(vsize, 'int8');
Elist = cell(vsize(3),1);

ne = size(E,1);

% put elements into Z buckets
for e = 1:ne
    Z = floor(C(E(e,:),3)) + 1;
    Z0 = min(Z);
    Z1 = max(Z) + 1;
    for z = Z0:Z1
        Elist{z} = [Elist{z} e];
    end
end

% Draw element-plane intersections for each z-plane
for z = 1:vsize(3)
    Ez = Elist{z};
    Zp = z - 0.5;
    for Zp = z-1:0.1:z
        for e = Ez
            % intersect element e with plane z = Zp
            ed = [ C(E(e,1),:) C(E(e,2),:)
                   C(E(e,2),:) C(E(e,3),:)
                   C(E(e,3),:) C(E(e,1),:)];
            pl = [];
            for i = 1:3
                z0 = ed(i,3);
                z1 = ed(i,6);
                if ((z0 < Zp && z1 > Zp) || (z0 > Zp && z1 < Zp))
                    p0 = ed(i,1:3);
                    p1 = ed(i,4:6);
                    dp = p1 - p0;
                    pi = p0 + dp * (Zp - z0)/dp(3);
                    pl = [pl; pi];
                end
                if (ed(i,3) == Zp)
                    pl = [pl; ed(i,1:3)];
                end
                if (ed(i,6) == Zp)
                    pl = [pl; ed(i,4:6)];
                end
            end
            np = size(pl,1);
            if (np > 1)
                for i = 1:np - 1
                    v1 = floor(pl(i, 1:2)) + 1;
                    v2 = floor(pl(i+1, 1:2)) + 1;
                    % draw line on z
                    dv = v2 - v1;
                    dvmax = max(abs(dv));
                    if (v1(1) > 0 && v1(2) > 0)
                        V(v1(1),v1(2),z) = 1;
                    end
                    for j = 1:dvmax
                        v = round(v1 + (dv * j /dvmax));
                        if (v(1) > 0 && v(2) > 0)
                            V(v(1),v(2),z) = 1;
                        end
                    end
                end
            end
        end
    end
end
[K,L,M]=size(V);
for i = 1:size(V,3)
    V(:,:,i) = imfill(V(:,:,i), 'holes');
end
for i = 1:size(V,2)
    V(:,i,:) = imfill(reshape(V(:,i,:),K,M), 'holes');
end
for i = 1:size(V,1)
    V(i,:,:) = imfill(reshape(V(i,:,:),L,M), 'holes');
end


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% generate patches


selection = get(handles.popupmenu1, 'Value');
% selection = 1 -> select patch size
% selection = 2 -> 10, 6, 3 mm 
% selection = 3 -> 10 mm 
% selection = 4 -> 6 mm
% selection = 5 -> 3 mm

curr_dir = pwd;
of = handles.MeshFolder;
cd(of)

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
rand('state', 0);                         %comment this statement if you want to produce random mesh every time

load FSss_cor
mesh = geodesic_new_mesh(Css(:,2:4),Ess(:,2:4));
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

if selection == 1
    disp('Please select patch size...');
elseif selection == 2
    %gaussian patches
    disp('Calculating patches with 10 mm...')
    ss_g10 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 10); 
    save ss_g10 ss_g10
    disp('Calculating patches with 10 mm is done!')
    disp('Calculating patches with 6 mm..')
    ss_g6 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 6);
    save ss_g6 ss_g6
    disp('Calculating patches with 6 mm is done!')
    disp('Calculating patches with 3 mm...')
    ss_g3 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 3);
    save ss_g3 ss_g3
    disp('Calculating patches with 3 mm is done!')
elseif selection == 3
    disp('Calculating patches with 10 mm...')
    ss_g10 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 10); 
    save ss_g10 ss_g10
    disp('Calculating patches with 10 mm is done!')
elseif selection ==4
    disp('Calculating patches with 6 mm..')
    ss_g6 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 6);
    save ss_g6 ss_g6
    disp('Calculating patches with 6 mm is done!')
elseif selection ==5    
    disp('Calculating patches with 3 mm...')
    ss_g3 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 3);
    save ss_g3 ss_g3
    disp('Calculating patches with 3 mm is done!')
end
    
cd(curr_dir)


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


selection = ((get(handles.popupmenu3, 'Value')));
% selection = 1 -> select source localization method
% selection = 2 -> SBL
% selection = 3 -> SCS 

comp_index = str2num(get(handles.edit1,'String'));

v = evalin('base','EEG');


% start source localization
curr_dir = pwd;
of = handles.MeshFolder;
cd(of)

EEG = handles.EEG;

Phi_EEG = EEG.icawinv(:,comp_index);
ndip = length(comp_index);

% Solve distributed source loc for all models
sensor_file = [handles.arg_subject '_' handles.arg_session '.sensors'];

sens = load(sensor_file, '-mat');
[ind_fp, ind_eeg] = find_indexes(EEG, sens.ind, sens.eloc);
Phi_EEG = Phi_EEG(ind_eeg,:);

LFM_name = [handles.arg_session '_LFM']; 
load(LFM_name)
LFM2 = LFM(ind_fp,:);
load ss_g10    
load FSss_cor    % sourcespace
load Node_area
max_scs_iter = 25;

if selection == 1
    disp('Please select the source locaaalization method...');
elseif selection == 2
    % SBL
    load ss_g6    
    load ss_g3    
    disp('SBL source localization started...')
    for ij = 1:size(Phi_EEG, 2)
        Vdata = Phi_EEG(:,ij);
        Vdata = Vdata - mean(Vdata);

        Jb = source_loc_SBL_gaus_function(Vdata, ss_g3, ss_g6, ss_g10, LFM2, 4, 0.001, 0);
        sourceJ(:,ij) = Jb;
                
        pot = LFM2 * Jb; pot = pot - mean(pot);
        diff = Vdata - pot; 
        fvalJ(ij) = sum(diff(:).^2) / sum(Vdata(:).^2);

    end
    handles.cort_source = sourceJ;
    handles.cort_fval = fvalJ;
    save cortex_source_sbl sourceJ fvalJ
    disp('SBL source localization finished...')
    
elseif selection == 3
    % SCS
    disp('SCS source localization started...')
    for ij = 1:size(Phi_EEG, 2)
        Vdata = Phi_EEG(:,ij);
        Vdata = Vdata - mean(Vdata);
        [Js1, Jit, fvx] = inverse_cov_sparse_average_noise_whole18_sparse_patchz2(LFM2, Vdata, ss_g10, max_scs_iter, 0);
        
        % find the most compact source in iterations
        for comi = 1:max_scs_iter+1
            pot = Jit(:,comi); pot = pot';
            compact0iter(comi) = calc_compactness(pot, An, Css(:,2:4), 1);
        end
        [maxcom, maxcomi] = max(compact0iter(5:max_scs_iter+1));
        maxcomiter = maxcomi + 4;

        sourceJ(:,ij) = Jit(:,maxcomiter);
                
        pot = LFM2 * Js1; pot = pot - mean(pot);
        diff = Vdata - pot; 
        fvalJ(ij) = sum(diff(:).^2) / sum(Vdata(:).^2);
    end
    handles.cort_source = sourceJ;
    handles.cort_fval = fvalJ;
    save cortex_source_scs sourceJ fvalJ
    disp('SCS source localization finished...')
end

cd(curr_dir)

guidata(hObject, handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

comp_plot = str2num(get(handles.edit4,'String'));
comp_plot = comp_plot(1);

conf = nft_get_config;

curr_dir = pwd;
of = handles.MeshFolder;
cd(of)

pot = handles.cort_source(:,comp_plot);
save pot pot -ascii

f = fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'nfield load %spot',of);
%fprintf(f, 'quit\n');
fclose(f);

a = sprintf('"%s" -c "%sStepSc.txt" FSss.smf', conf.showmesh3, of);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

cd(curr_dir)


function [Coord, Elem]=ReadSMF(name,x,y,z,sc)
fid=fopen(name, 'r');
nnp=0; nel=0;
line=1;
while line~=-1
   line=fgets(fid);
	if line(1)=='v';
   	nnp=nnp+1;
	   [A,count,ERRMSG,NEXTINDEX] = sscanf(line,'%c %f %f %f',4);
   	Coord(nnp,1)=nnp;
	   Coord(nnp,2:4)=A(2:4)';
	elseif (line(1)=='t')|(line(1)=='f');
   	nel=nel+1;
	   [A,count,ERRMSG,NEXTINDEX] = sscanf(line,'%c %d %d %d',4);
	  	Elem(nel,1)=nel;
   	Elem(nel,2:4)=A(2:4)';
   end
end
fclose(fid);

Coord(:,4)=Coord(:,4)/sc+z;
Coord(:,2)=Coord(:,2)+x;
Coord(:,3)=Coord(:,3)+y;

% Vfs = mesh2vol2(v5, E1(:,2:4)); % freesurfer brain surface
function Segm = new_segm(Segm,Vfs);
Vbr0 = Segm.brainmask;
Vc0 = Segm.innerskullmask;
Vsk0 = Segm.outerskullmask;
Vsc0 = Segm.scalpmask;

s(1,:) = size(Vfs);
s(2,:) = size(Vbr0);
sz = max(s);

Vfs0 = zeros(sz, 'int8');
Vfs0(1:s(1,1), 1:s(1,2), 1:s(1,3)) = Vfs;

se1 = strel('ball',3, 3);
se = strel('ball', 3, 3, 0);
%se = strel('ball', 3, 3);

% new brain volume not intersecting with FS brain
Vfs2 = imclose3D(Vfs0, se1);
Vfs2 = imdilate3D(int8(Vfs2), se);
Vfs2 = Vfs2 - min(min(min(Vfs2)));
Vbr2 = Vfs2 | Vbr0; 
clear Vfs2

% new csf volume not intersecting with brain
Vbr3 = imdilate3D(int8(Vbr2), se);
Vbr3 = Vbr3 - min(min(min(Vbr3)));
Vc2 = Vc0 | Vbr3;
clear Vbr3

% new skull volume not intersecting with csf
Vc3 = imdilate3D(int8(Vc2), se);
Vc3 = Vc3 - min(min(min(Vc3)));
Vsk2 = Vsk0 | Vc3;
clear Vc3

% new skull volume not intersecting with csf
Vsk3 = imdilate3D(int8(Vsk2), se);
Vsk3 = Vsk3 - min(min(min(Vsk3)));
Vsc2 = Vsc0 | Vsk3;
clear Vsk3

% farklara bak
A = int8(Vbr0) - int8(Vbr2); 
clear Vbr0
fbr = sum(sum(sum(abs(A))));
A = int8(Vc0) - int8(Vc2);
clear Vc0
fc = sum(sum(sum(abs(A))));
A = int8(Vsk0) - int8(Vsk2);
clear Vsk0
fsk = sum(sum(sum(abs(A))));
A = int8(Vsc0) - int8(Vsc2);
clear Vsc0
fsc = sum(sum(sum(abs(A))));
[fbr fc fsk fsc]
% set the new volumes in Segm structure
Segm.innerskullmask = Vc2;
Segm.outerskullmask = Vsk2;
Segm.brainmask = Vbr2;
Segm.scalpmask = Vsc2;
clear Vsc2 Vbr2 Vsk2 Vc2
% Vfs = mesh2vol2(v5, E1(:,2:4)); % freesurfer brain surface




function [so2, k1,k2] = correct_source_space(so, C1, E1)

% so : sourcespace C1, E1 : mesh
% corrects source space
[N, M] = element_normals(C1, E1);
so2 = so;
[dim, inm] = CheckSourceSpace(so(:,1:3), C1, E1, 2);
k1 = find(inm == 0); %k1x = find(inm==1);
%if length(k1x) < length(k1)
%    k1 = k1x;
%end
no_intnodes = length(k1)

k2 = find(dim < 2); % find the nodes closer than 1mm
k2 = setdiff(k2, k1);
no_closenodes = length(k2)

if length(no_intnodes)>0
    
    for iter=1:2
    for i = 1:no_intnodes
        p1 = so2(k1(i),1:3);
        [dm, Pm, el, in] = DistMeshPoint2(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        % pn = Pm+2*nor;
        pn = Pm-2*N(el,:);
        [dm1, Pm1, el1, in1] = DistMeshPoint2(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k1(i)
        end
        so2(k1(i),1:3) = pn;
    end
    end
end

if length(no_closenodes)>0
    for i = 1:no_closenodes
        p1 = so2(k2(i),1:3);
        [dm, Pm, el, in] = DistMeshPoint2(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        pn = Pm-nor*2;
        [dm1, Pm1, el1, in1] = DistMeshPoint2(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k2(i)
        end
        so2(k2(i),1:3) = pn;
    end
end


function [dim, inm] = CheckSourceSpace(so,C,E,thr)
% so is the source space
% C,E is the linear mesh
% dim gives the vector of distances
% inm gives the vector of inside sources (logical)

%hh = waitbar(0,'computing...');
M = size(so,1);
for i = 1 : M
 %   waitbar(i/M)
 if round(i/1000) == i/1000
     i
 end
    [dm, Pm, el, in] = DistMeshPoint2(so(i,:), C, E);
    dim(i) = dm;
    inm(i) = in;
end
%close(hh);
k = find(inm == 0);   % dipoles outside the mesh
l = find(dim < thr);    % dipoles closer to the mesh less than 1mm
m = setdiff(k, l);   % dipoles outside the mesh, closer to the mesh less than 1mm 
    



function [dm,Pm,el,in]=DistMeshPoint2(P,Coord,Elem);
% looks for if P is inside the mesh Coord, Elem or not
% Pm is the point on the mesh
% dm is the distance
% el is the element of Pm
% in = inside (bool)
% works for LINEAR MESH

% Coord: mesh vertex coordinates (first colum shows the indices)
% Elem: mesh connectivity (first colum shows the indices)


nnp=size(Coord,1);

r=5;
N=0;
Xo=Coord(:,2:4)-ones(nnp,1)*P;
Rad=sum(Xo.*Xo,2);
while length(N)<3
    r=r+5;
    N=find(Rad<r^2);
end

% find the neighbour elements of the closest nodes
E = ElementsOfTheNodes(Coord,Elem,N);


% find the intersection of the line PP1 with the elements of E
Pint=[]; dis=[];
for i=1:length(E)
   Pa=Coord(Elem(E(i),2),2:4);
   Pb=Coord(Elem(E(i),3),2:4);
   Pc=Coord(Elem(E(i),4),2:4);
   [D,Pp]=DistTrianglePoint2(Pa,Pb,Pc,P);
   dis(i)=D;
   Pint(i,:)=Pp;
end
[j,k]=min(abs(dis));
Pm=Pint(k,:);
dm=abs(dis(k));
el=E(k);

N2 = Elem(el,2:4);
el2 = ElementsOfTheNodes(Coord,Elem,N2);

% check if P is inside or outside

% find the normal vector of the element
v1=Coord(Elem(el2,3),2:4)-Coord(Elem(el2,2),2:4);
v2=Coord(Elem(el2,4),2:4)-Coord(Elem(el2,3),2:4);
n=cross(v1, v2); n2=mean(n);
Norm = n2/norm(n2);

% if the vector PPm.Norm > 0 inside, < 0 outside
if dot(Pm-P, Norm) > 0 
    in = 1; % P is inside the mesh
else
    in = 0;
end



function [D,Pp]=DistTrianglePoint2(Pa,Pb,Pc,Px)
% finds the minimum distance of a point Px to triangle Pa, Pb, Pc
% difference from DistTrianglePoint 
% doesn't look if the projection of the point is in the triangle or on the edge 
% of the triangle otherwise MinD is 1000

% find the minimum distance of the point with the 
% plane which is formed by the triangle
% find the normal of the plane
eps=1e-4;
v1=Pa-Pb;
v2=Pc-Pb;
n=[v1(2)*v2(3)-v2(2)*v1(3) v2(1)*v1(3)-v1(1)*v2(3) v1(1)*v2(2)-v1(2)*v2(1)];
n=n/norm(n);
d=-n(1)*Pa(1)-n(2)*Pa(2)-n(3)*Pa(3);
% the plane equation is n(1)*x+n(2)*y+n(3)*z+d=0
% distance of the point to the plane is D
D=(n(1)*Px(1)+n(2)*Px(2)+n(3)*Px(3)+d)/sqrt(n(1)^2+n(2)^2+n(3)^2);
% point on the plane
Pp=Px-D*n;
% check if the point is on the triangle
% Determine whether or not the intersection point is bounded by pa,pb,pc 
Pa1=Pa-Pp;
normPa1=norm(Pa1);
if normPa1>eps
   % normalize the unit vectors
   Pa1=Pa1/normPa1;  
end
Pa2 = Pb - Pp;
normPa2=norm(Pa2);
if normPa2>eps
   Pa2=Pa2/normPa2; 
end
Pa3 = Pc - Pp;
normPa3=norm(Pa3);
if normPa3>eps
   Pa3=Pa3/normPa3;
end
%the angles are 
a1 = acos(Pa1(1)*Pa2(1) + Pa1(2)*Pa2(2) + Pa1(3)*Pa2(3));
a2 = acos(Pa2(1)*Pa3(1) + Pa2(2)*Pa3(2) + Pa2(3)*Pa3(3));
a3 = acos(Pa3(1)*Pa1(1) + Pa3(2)*Pa1(2) + Pa3(3)*Pa1(3));

total = a1+a2+a3;
% if total is 2*pi then the point is in the triangle or on the edges
if (abs(total - 2*pi) < eps)
   MinD=abs(D);
else
    % find the distance of Px to the edges
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pb);
   Dix(1)=Di; Ppi(1,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pc);
   Dix(2)=Di; Ppi(2,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pc,Pb);
   Dix(3)=Di; Ppi(3,:)=Ppl;
   [u,v]=min(Dix);
   MinD=u;
   Pp=Ppi(v,:);
end
D=MinD;

function [D,Pp]=DistPointLineSegment(P,P1,P2);
% Finds distance between the point P and the line segment P1P2

v1=P2-P1;
n1 = v1/norm(v1);

v2 = P-P1;

% projection (component) of v2 along v1 (n1)
vp = dot(n1, v2);

% make sure it is in the range
if (vp < 0)
    vp = 0;
elseif (vp > norm(v1))
    vp = norm(v1);
end


% closest point on edge
Pp = P1 + vp * n1;

D = norm(P-Pp);


function E=ElementsOfTheNodes(Coord,Elem,A);
% A is the list of the nodes
% E is the list of the elements that have nodes in A
E=[];
nop=size(Elem,2);
if nop==4
   for i=1:length(A)
   	n1=find(Elem(:,2)==A(i));
	   n2=find(Elem(:,3)==A(i));
   	n3=find(Elem(:,4)==A(i));
	   n4=union(n1,n2);
   	n5=union(n3,n4);
	   E=union(E,n5);
      clear n1 n2 n3 n4 n5
   end
elseif nop==7
   for i=1:length(A)
   	n1=find(Elem(:,2)==A(i));
	   n2=find(Elem(:,3)==A(i));
      n3=find(Elem(:,4)==A(i));
      n4=find(Elem(:,5)==A(i));
   	n5=find(Elem(:,6)==A(i));
   	n6=find(Elem(:,7)==A(i));
	   n7=union(n1,n2);
      n8=union(n7,n3);
      n9=union(n8,n4);
      n10=union(n9,n5);
      n11=union(n10,n6);
	   E=union(E,n11);
      clear n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11
   end
end


function [N, M] = element_normals(Coord, Elem)

if size(Elem,2)==7

n1 = Coord(Elem(:,2),2:4);
n2 = Coord(Elem(:,4),2:4);
n3 = Coord(Elem(:,6),2:4);
elseif size(Elem,2)==4
n1 = Coord(Elem(:,2),2:4);
n2 = Coord(Elem(:,3),2:4);
n3 = Coord(Elem(:,4),2:4);
end    

M = (n1 + n2 + n3) / 3;

n1 = n1 - n3;
n2 = n2 - n3;

N = cross(n1, n2);


for i=1:size(Elem,1)
    N(i,:)=N(i,:)/norm(N(i,:));
end



function WriteSMF2(name,Coord,Elem);
nnp=size(Coord,1); nel=size(Elem,1);
fid=fopen(name, 'w');
fprintf(fid,'v %f %f %f \r\n',Coord');
fprintf(fid,'t %d %d %d \r\n',Elem');
fclose(fid);


function [Ae, An] = area_of_nodes(C,E);
% Ae : area of elements
% An : area of nodes

Ne = length(E);

for i = 1:Ne
    X = C(E(i,2:4),2:4);
  
    AB = X(1,:)-X(2,:);
    AC = X(1,:)-X(3,:);
    Ae(i) = 1/2 * norm(cross(AB,AC));
end

Nn = length(C);
for i = 1:Nn
    e1 = ElementsOfTheNodes(C,E,i);
    An(i) = mean(Ae(e1));
end
   


function Nn = NodeNormals(Coord,Elem,ccw);
% if ccw=1, then the mesh is CCW
% if ccw=0, mesh is CW

% calculates normals of nodes for a linear mesh

En = ElemNormal(Elem, Coord);
nr = length(Coord(:,1));

hh = waitbar(0,'computing node normals...');
for i = 1:nr
    waitbar(i/nr)

    %ns = findNeigNodes(Coord, Elem, i);
    ns = ElementsOfTheNodes(Coord,Elem,i);
    if length(ns)<3
        nn = findNeigNodes(Coord,Elem,i);
        ns = ElementsOfTheNodes(Coord,Elem,nn);
        if length(ns)<3
            K = Coord(:,2:4)-ones(length(Coord),1)*Coord(i,2:4);
            M = sqrt(sum(K.*K,2));
            [Y,I] = sort(M);
            j=1;
            while length(ns)<3
                ns = ElementsOfTheNodes(Coord,Elem,I(j));
                j=j+1;
            end
        end
    end
    Nn(i,:) = mean(En(ns,:),1);
    Nn(i,:) = Nn(i,:)/norm(Nn(i,:));
end
close(hh)
if ccw==0
    Nn = -Nn;
end

% correct if there are NaNs in the calculation of normals

w = isnan(Nn);
[I,J] = find(w==1);
if length(I)>0
    for i = 1:length(I)
        iw = [];
        nn = I(i);
        while length(iw) == 0
            nn = findNeigNodes(Coord,Elem,nn);
            Nni = Nn(nn,:);
            w =~ isnan(Nni);
            [iw, jw] = find(w==1);
        end
        Nn(I(i),:) = Nni(iw(1),:);
    end
end


function Norm=ElemNormal(Elem,Coord);

ne = length(Elem(:,1));
nr = length(Coord(:,1));

for k = 1:ne
   v1 = Coord(Elem(k,3),2:4)-Coord(Elem(k,2),2:4);
   v2 = Coord(Elem(k,4),2:4)-Coord(Elem(k,3),2:4);
   n = cross(v1, v2);
   Norm(k,1:3) = n/norm(n);
end

function N=findNeigNodes(Coord,Elem,n);
% find the neighbour nodes of the nth node of the mesh

nop=size(Elem,2);
% find the neighbouring elements
if nop==4
    E=[];
    for i=1:length(n)
        E2=find(Elem(:,2)==n(i) | Elem(:,3)==n(i) | Elem(:,4)==n(i));
        E = union(E,E2);
    end

    f=Elem(E,2:4);
    N=unique(f(:));
    
elseif nop==7
    E=[];
    for i=1:length(n)
        E2=find(Elem(:,2)==n(i) | Elem(:,3)==n(i) | Elem(:,4)==n(i) | Elem(:,5)==n(i) | Elem(:,6)==n(i) | Elem(:,7)==n(i));
        E = union(E,E2);
    end

    f=Elem(E,2:7);
    N=unique(f(:));
end
    


function [ind_fp, ind_eeg] = find_indexes(EEG, elp_index, eloc);
% realistic data icin

A = EEG.icawinv;

% neglect the FID electrodes
y = elp_index;
N = length(y);
for i = 1:N
    elocn(i).labels = eloc(y(i)).labels;
    elocn(i).X = eloc(y(i)).X;
    elocn(i).Y = eloc(y(i)).Y;
    elocn(i).Z = eloc(y(i)).Z;
 %   elocn(i).type = eloc(y(i)).type;
end

Nel2 = length(elocn);
Neeg = length(EEG.chanlocs);

Mel2 = zeros(1, Nel2);
Meeg = zeros(1, Neeg);

% Mel2(i) = index of EEG.chanlocs that correspond to electrode i of eloc2
% Meeg(i) = index of elocn that correspond to electrode i of EEG.chanlocs
clear Mel2 Meeg
for i = 1:Nel2
    for j = 1:Neeg
        if strcmp(elocn(i).labels, EEG.chanlocs(j).labels)
            Mel2(i) = j;
            Meeg(j) = i;
            continue;
        end
    end
end

ind_fp = find(Mel2>0); % index for the FP outputs (LFM, TM, session)
ind_eeg = find(Meeg>0); % index for the EEG structure (ICs)

function [J,Jit, fval,stdd_log_a,J_s,dispact_s,prob_ts] = inverse_cov_sparse_average_noise_whole18_sparse_patchz2(F,P,ss_MNI_gaussion,max_it,flag)

%J=inverse_cov_sparse_average(F,P,voxel_position)
%F is the lead fied matrix, P is the observed scalp potential and
% ss_MNI_gaussion for 6mm or 10mm
%max_it = 30, 20
% flag = 1
%voxel_position is the locatoin of the dipoles 
%Edited by Cheng Cao 2011
% A compact function is added 
% the covariance matrix is updated 
%Parallel computation is used 
% modified based on version 7, keep the hidden elements
%dealt with the noise issue, considering the DC shift of the noise
%Using two-point stepsize gradient
%Use log(std) to achieve better performance
%form version 10
%set the initial nsr_level according to the eigen value of M
%Under develovelpment
%Add smooth matrix Mar21,2012smooth_control
%change pinv to inv MAr 22

% initialize
[rt,ty] = size(F);
fval = zeros(1,max_it+1);
Jit = zeros(ty,max_it+1);

stop = 1;
step_size = 0.01; %0.01;
minium_nsr = 0.1; %0.1
n_control_para = 0.00001;
p_std_cof = 0;                              
nsr_coefi = 0.00005; %0.001;%0.01
%nsr_coefi = 0.1; %0.001;%0.01

[number_electrode,number_voxel] = size(F);
smooth_control = 0.1;
J_s = [];
dispact_s = [];
J = ones(number_voxel,1);
prob_ts = [];
MIN_GAMMA = 1e-16; 

P = P - mean(P);
P = P(1:number_electrode-1);
P_norm = norm(P);
P = P / P_norm;

F = F - ones(number_electrode,1) * mean(F);
F = F(1:number_electrode-1,:);
F_norms_sqr = (sum(F.^2))';

for i = 1:number_voxel
    F(:,i) = F(:,i)*((F_norms_sqr(i)).^-0.5);
end

F_e = zeros(number_electrode-1,number_electrode-1);
F_e = F_e-1 / number_electrode;


for i = 1:number_electrode-1
    F_e(i,i) = F_e(i,i)+1;
end

F_e = F_e / norm(F_e(:,1));

M = zeros(number_electrode-1,number_electrode-1);
N = zeros(number_voxel,1);
cov_column = zeros(number_voxel,1);
voxel_position_a = zeros(number_voxel,1);
pre_decompact = inf;

stdd_log_a = zeros(number_voxel,1);
nsr_level = zeros(number_electrode-1,1);
stdd_log_a = stdd_log_a+p_std_cof*log(F_norms_sqr)-0.1*log(min(F_norms_sqr));
n_it = 1;
err = zeros(1,number_electrode-1);

g_k = zeros(number_voxel+number_electrode-1,1);
g_k_old = g_k;
stdd_log_old = stdd_log_a;
nsr_level_old = nsr_level;
prob = 0;
   
J_index = 1:number_voxel;
number_a_voxel = number_voxel;
F_a = F;
n_itt = 1;

ss_matrix = ss_MNI_gaussion;
for i=1:number_voxel
    ss_MNI_gaussion(:,i) = ss_matrix(:,i) / sqrt(sum(ss_matrix(:,i).^2));
end
ss_MNI_gaussion = ss_MNI_gaussion';


while stop && n_it <= max_it
    
    row_temp = zeros(number_electrode-1,number_voxel);
    ss_diag_diag = sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
    
    row_temp = F_a * ss_diag_diag * ss_MNI_gaussion;
    M = row_temp * row_temp';
    a1 = isnan(M);
    if sum(sum(a1))>1
        break
    end
    [~,D_m] = eig(M);
    if n_it == 1
        lam_max = max(diag(D_m));
        nsr_level_init = nsr_coefi*mean(diag(D_m));%Changed on Mar 19 2012
        scale = minium_nsr/nsr_level_init; 
      
        M = scale * M;
        row_temp = row_temp*sqrt(scale);
        stdd_log_a = stdd_log_a+0.5*log(scale);
        nsr_level = zeros(number_electrode-1,1)+log(minium_nsr)/2;
    end
    snr = mean(diag(D_m)) / exp(2*mean(nsr_level));
    M = M + F_e * diag(exp(2*nsr_level)) * F_e';
    [~,D_m]=eig(M);
    c_M = cond(M);
    min_eigv = min(diag(D_m));
    prob = P' * inv(M) * P;
    prob_t = (number_electrode-1) * log(P'*inv(M)*P) - log(det(inv(M))) - n_control_para * (mean(stdd_log_a) + mean(nsr_level));%%add the noise control;
    prob_ts = [prob_ts,prob_t];
    inv_M = inv(M);
   
    ss_diag_diag=sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
   
    %%%%%%%%%%%%% calculate the current %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag == 1
        N = F_a' * inv(M) * P;
        J_tmp = zeros(number_voxel,1);
        J_tmp = (ss_diag_diag*ss_MNI_gaussion)' * N;
        J = ss_diag_diag * ss_MNI_gaussion * J_tmp;
        rv = std(F*J-P)
        J = P_norm*J./sqrt(F_norms_sqr);

        if n_it > 1
            trange_J = norm(J-J_s(:,end)) / norm(J)
        end
        J_s = [J_s, J];
    end

    %%%%%%%%%%% starts to calcualte the gradient %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    g_k = [];
    trace_inM = diag(inv_M);
    cov_init = (ss_MNI_gaussion);
    temp_PM=P'*inv_M;

    temp_PM_F=temp_PM*F;
    temp_cov_rowtemp_temp_PM=cov_init*row_temp'*temp_PM';
    g_1=-(number_electrode-1)*(2*temp_PM_F'.*temp_cov_rowtemp_temp_PM.*exp(stdd_log_a))/(prob);

    %temp_cov_rowtemp=cov_init*row_temp';
    temp_cov_rowtemp=row_temp*cov_init';
    g_2=sum(inv_M*F.*temp_cov_rowtemp);

    g_3=sum(inv_M*temp_cov_rowtemp.*F);
    g_k=[g_k;g_1+(g_3+g_2)'.*exp(stdd_log_a)];
    g_k=g_k+n_control_para/(number_voxel);

    for i=number_a_voxel+1:number_a_voxel+number_electrode-1
        ro_temp=zeros(1,number_electrode-1);
        ro_temp(i-number_a_voxel)=1;
        g_k(i)=-2*(number_electrode-1)*P'*inv_M*F_e*diag(ro_temp)*F_e'*inv_M*P*(exp(2*nsr_level(i-number_a_voxel)))/(prob)+2*trace(inv_M*F_e*diag(ro_temp)*F_e')*(exp(2*nsr_level(i-number_a_voxel)))-n_control_para/(number_electrode);
    end

    %fprintf('finished the calucation of the gradient\n');
    %%%%%%%%%%%%%%%%% transfer the gradient %%%%%%%%%%%%%%%%%%%%%%%%
    g_k_t=[];
    cov_init=ss_MNI_gaussion;
    g_k_t=[g_k_t;cov_init*g_k(1:number_voxel)];
    g_k_t=[g_k_t;g_k(number_a_voxel+1:end)];
    norm(g_k_t);
    surface_norm = ones(number_voxel,1);
    surface_norm = surface_norm/norm(surface_norm);

    g_k_t = g_k_t-mean(g_k_t);
    g_k = g_k'*g_k_t*g_k_t / (norm(g_k_t))^2;

    if n_itt >= 3
        step_size = min(0.5*sqrt(((stdd_log_a-stdd_log_old)'*(stdd_log_a-stdd_log_old)+(norm(nsr_level-nsr_level_old))^2)/((g_k-g_k_old)'*(g_k-g_k_old))),1000);
    end

    stdd_log_old = stdd_log_a;   
    nsr_level_old = nsr_level;

    stdd_log_a = stdd_log_a-step_size*g_k(1:number_voxel);
    nsr_level = nsr_level-step_size*g_k(number_voxel+1:end);

    M = zeros(number_electrode,number_electrode);

    g_k_old = g_k;
    n_it = n_it+1;
    n_itt = n_itt+1;
    %fprintf('%d iter  %d voxels used snr_level= %d current gradient',n_it-1,number_a_voxel,max(exp(nsr_level)),norm(g_k));
    std_max = max([stdd_log_a;nsr_level]);
    if std_max> 5 
        nsr_level = nsr_level-std_max+1;
        stdd_log_a = stdd_log_a-std_max+1;
        stdd_log_old = stdd_log_old-std_max+1;
        nsr_level_old = nsr_level_old-std_max+1;
        %fprintf('sum of std changed\n');
    end
    if norm(g_k) < 0.5   %||isinf(abs(prob_t))
        %stop = 0;
    end

   if norm(g_k) > 1000  % zeynep singular matrix oluyor
        stop = 0;
    end
  
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ss_diag_diag = sparse(1:number_voxel,1:number_voxel,exp(stdd_log_a));
    row_temp = zeros(number_electrode-1,number_voxel);

    cov_init = ss_MNI_gaussion;
    cov_init = ss_diag_diag*cov_init;

    row_temp = F_a * cov_init;
    M = row_temp * row_temp';
    M = M + F_e * diag(exp(2*nsr_level)) * F_e';
    N = F_a' * inv(M) * P;

    J_tmp = zeros(number_voxel,1);
    J_tmp = cov_init' * N;
    J = cov_init*J_tmp;

    err = F*J - P;
    fval(n_it) = sum(err(:).^2) / sum(P(:).^2);
    J = P_norm*J./sqrt(F_norms_sqr);
    Jit(:,n_it) = J;
    %if n_it > 5
    %    if fval(n_it)>fval(n_it-1)
    %        J = Jit(:,n_it-1);
    %        break;
    %    end
    %end

end
J_s = [J_s,J];


function compact = calc_compactness(pot, An, vr, flag, geodis)
% Cr, Er : cortical mesh
% An = area mesh nodes
% pot: cortical source
% vr: voxel_position

%flag = 0; % Cheng's method
%flag = 1; % Zunic's Kst
%flag = 2; % my area method
if nargin < 5
    geodis = [];
end

An  = An(:); An = An';
if flag == 0 
    k = 0.1;
    max_pot = max(abs(pot));
  
    ind_large_value = find(abs(pot) > max_pot/10);
     
    number_main_voxel = length(ind_large_value);
    if  number_main_voxel > 0.05*length(pot)
        compact = inf;
    else
        compact = 0;
        for i = 1:number_main_voxel-1
            for j = i+1:number_main_voxel
                te = exp(k*norm(vr(ind_large_value(j),:)-vr(ind_large_value(i),:)));
                compact = compact + abs(pot(ind_large_value(i))/max_pot) * abs(pot(ind_large_value(j))/max_pot) * te;
            end
        end
        compact = compact / (number_main_voxel*number_main_voxel);
    end
end



if flag == 3 
    k = 0.1;
    max_pot = max(abs(pot));
  
    ind_large_value = find(abs(pot) > max_pot/10);
     
    number_main_voxel = length(ind_large_value);
    if  number_main_voxel > 0.05*length(pot)
        compact = inf;
    else
        compact = 0;
        for i = 1:number_main_voxel-1
            for j = i+1:number_main_voxel
                if i > j
                    dist = geodis(i,j);
                else
                    dist = geodis(j,i);
                end
                if dist == 0
                    dist = 1000;
                end
                if i == j
                    dist = 0;
                end
                te = exp(k * dist);
                compact = compact + abs(pot(ind_large_value(i))/max_pot) * abs(pot(ind_large_value(j))/max_pot) * te;
            end
        end
        compact = compact / (number_main_voxel*number_main_voxel);
    end
end

if flag == 1
    % Compactness measure for 3D shapes
    % J Zunic, K. Hirota, C. Martinez-Ortiz, 2012
    pot = abs(pot);
      if max(pot) == 0
          compact = 0;
          return
      end
 
    pot = pot / max(pot) * 100; % max =100

    max_pot = max(abs(pot));
    bl = 10;
    ind_large_value = find(abs(pot) > max_pot/bl);
    while length(ind_large_value) < 2
        bl = bl*2;
        ind_large_value = find(abs(pot) > max_pot/bl);
    end
    number_main_voxel = length(ind_large_value);
    pot = abs(pot);
    %if  number_main_voxel > 0.05*length(pot)
        %compact = inf;
    %else
        vol_S = sum(An(ind_large_value).*pot(ind_large_value));
        surf_S = sum(An(ind_large_value));
        compact = 36 * pi * vol_S^2 / surf_S^3;
    %end
    compact = compact/max(pdist(vr(ind_large_value,:)));
end



if flag == 2
    pot = abs(pot);
    if max(pot) == 0
          compact = 0;
          return
      end
    pot = pot / max(pot) * 100; % max =100
    ap  = round(prctile(pot,99)); 
    if ap < 1; ap = 1; end
    ni = find(pot > ap);
    clear ai ao
    for i = ap:100
        
        ni = find(pot > i);
        ai(i) = sum(An(ni));
        ao(i) = sum(An(ni) .* pot(ni));
    end
    k = find(ao>std(ao)); jk = max(k)+1;
    compact = (1-ao(jk)/ao(ap)) * 100;
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function [J_esmtime, Jt, fval]  = source_loc_SBL_gaus_function(p,ss_MNI_3_gaus,ss_MNI_6_gaus,ss_MNI_10_gaus,LFM_MNI_sourcespaceNn,max_inter,nsr_lv, flag1)

% if flag1 = 0, EM, if flag1 = 2; McKay

if nargin < 8
    flag1 = 0; % EM
end


LFM = LFM_MNI_sourcespaceNn;
%nsr_lv=0.0001;

ss_10 = sum(ss_MNI_10_gaus,2);
ss_6 = sum(ss_MNI_6_gaus,2);
ss_3 = sum(ss_MNI_3_gaus,2);
ss_10_spnorm = ss_MNI_10_gaus';
ss_6_spnorm = ss_MNI_6_gaus';
ss_3_spnorm = ss_MNI_3_gaus';
ii = length(ss);
for i = 1:ii;
    ss_10_spnorm(:,i) = ss_10_spnorm(:,i) / ss_10(i);
    ss_6_spnorm(:,i) = ss_6_spnorm(:,i) / ss_6(i);
    ss_3_spnorm(:,i) = ss_3_spnorm(:,i) / ss_3(i);
end
ss_10_spnorm = ss_10_spnorm';
ss_6_spnorm = ss_6_spnorm';
ss_3_spnorm = ss_3_spnorm';

New_lfm_10n = LFM * ss_10_spnorm;
New_lfm_6n = LFM * ss_6_spnorm;
New_lfm_3n = LFM * ss_3_spnorm;

%lfm2 = [New_lfm_10 New_lfm_6 New_lfm_3];
lfm3 = [New_lfm_10n New_lfm_6n New_lfm_3n];
ssx = [ss_10_spnorm' ss_6_spnorm' ss_3_spnorm'];

clear New_lfm_* ss_* LFM


[mut3,mu,dmu,kk,gamma,fval] = sparse_learning_ss(lfm3, p, nsr_lv, max_inter, flag1, 0, 0, 1); % EM

itern = size(mut3,1);
vn = mut3(itern,:);
J_esmtime = ssx * vn';
Jt = ssx * mut3'; 
   


function [mut, mu,dmu,k,gamma,fval] = sparse_learning_ss(Phi,T,lambda,iters,flag1,flag2,flag3, nofig)
% *************************************************************************
% 
% *** PURPOSE *** 
% Implements generalized versions of SBL and FOCUSS for learning sparse
% representations from possibly overcomplete dictionaries.
%
%
% *** USAGE ***
% [mu,dmu,k,gamma] = sparse_learning(Phi,T,lambda,iters,flag1,flag2,flag3);
%
%
% *** INPUTS ***
% Phi       = N X M dictionary
% T         = N X L data matrix
% lambda    = scalar trade-off parameter (balances sparsity and data fit)
% iters     = maximum number of iterations
%
% flag1     = 0: fast Mackay-based SBL update rules
% flag1     = 1: fast EM-based SBL update rule
% flag1     = 2: traditional (slow but sometimes better) EM-based SBL update rule
% flag1     = [3 p]: FOCUSS algorithm using the p-valued quasi-norm
%
% flag2     = 0: regular initialization (equivalent to min. norm solution)
% flag2     = gamma0: initialize with gamma = gamma0, (M X 1) vector
%
% flag3     = display flag; 1 = show output, 0 = supress output
%
% *** OUTPUTS ***
% mu        = M X L matrix of weight estimates
% dmu       = delta-mu at convergence
% k         = number of iterations used
% gamma     = M X 1 vector of hyperparameter values
%
%
% *************************************************************************
% Written by:  David Wipf, david.wipf@mrsc.ucsf.edu
% *************************************************************************
    

% *** Control parameters ***
MIN_GAMMA       = 1e-16;  % 1e-4
MIN_DMU         = 1e-12;  
MAX_ITERS       = iters;
DISPLAY_FLAG    = flag3;     % Set to zero for no runtime screen printouts


% *** Initializations ***
[N M] = size(Phi); 
[N L] = size(T);

if (~flag2)         gamma = ones(M,1);    
else                gamma = flag2;  end;   
 
keep_list = [1:M]';
m = length(keep_list);
mu = zeros(M,L);
dmu = -1;
k = 0;

iter=0; % zeynep
fig = 0;
if nargin < 8
   figure
else
    fig=1;
end

% *** Learning loop ***
while (1)
iter=iter+1; % zeynep

    % *** Prune things as hyperparameters go to zero ***
    if (min(gamma) < MIN_GAMMA )
		index = find(gamma > MIN_GAMMA);
		gamma = gamma(index);
		Phi = Phi(:,index);
		keep_list = keep_list(index);
        m = length(gamma);
     
        if (m == 0)   break;  end;
    end;
    
    
    % *** Compute new weights ***
    G = repmat(sqrt(gamma)',N,1);
    PhiG = Phi.*G; 
    [U,S,V] = svd(PhiG,'econ');
    
    [d1,d2] = size(S);
    if (d1 > 1)     diag_S = diag(S);  
    else            diag_S = S(1);      end;
    
    U_scaled = U(:,1:min(N,m)).*repmat((diag_S./(diag_S.^2 + lambda + 1e-16))',N,1);       
    Xi = G'.*(V*U_scaled'); 
        
    mu_old = mu;
    mu = Xi*T; 

    temp = zeros(M,L);
    if (m > 0) temp(keep_list,:) = mu;  end;
    mut(iter,:,:) = temp; % zeynep
    di = ceil(sqrt(iters));
    if fig==0
        subplot(di,di,iter); plot(mu); % zeynep
    end
    
    pot = Phi * mu;
    diff = T - pot;
    err = sum(diff(:).^2) / sum(T(:).^2);
    fval(iter) = err;
    
    % *** Update hyperparameters ***
    gamma_old = gamma;
    mu2_bar = sum(abs(mu).^2,2);
    
    if (flag1(1) == 0)
        % MacKay fixed-point SBL
        R_diag = real( (sum(Xi.'.*Phi)).' );
        te = L*R_diag;
        if min(te) < MIN_GAMMA  % zeynep
            ind_te = find(te<MIN_GAMMA);
            te(ind_te) = MIN_GAMMA;
        end
        gamma = mu2_bar./te;  
        
    elseif (flag1(1) == 1)
        % Fast EM SBL
        R_diag = real( (sum(Xi.'.*Phi)).' );
        gamma = sqrt( gamma.*real(mu2_bar./(L*R_diag)) ); 
        
    elseif (flag1(1) == 2)
        % Traditional EM SBL
        PhiGsqr = PhiG.*G;
        Sigma_w_diag = real( gamma - ( sum(Xi.'.*PhiGsqr) ).' );
        gamma = mu2_bar/L + Sigma_w_diag;
        
    else
        % FOCUSS
        p = flag1(2);
        gamma = (mu2_bar/L).^(1-p/2);
    end;
    
    
    
    % *** Check stopping conditions, etc. ***
  	k = k+1;   
    if (DISPLAY_FLAG) disp(['iters: ',num2str(k),'   num coeffs: ',num2str(m), ...
            '   gamma change: ',num2str(max(abs(gamma - gamma_old))), ...
            '   fval: ',num2str(err)]); end;    
    
    if (k >= MAX_ITERS) break;  end;
    
    % zeynep
    if iter>5
    if (abs(fval(iter-1)-fval(iter)) < 0.0001) break; end; % zeynep 6/11/15
    end
    %
    
	if (size(mu) == size(mu_old))
        dmu = max(max(abs(mu_old - mu)));
        if (dmu < MIN_DMU)  break;  end;
    end;
   
end;


% *** Expand weights, hyperparameters ***
temp = zeros(M,1);
if (m > 0) temp(keep_list,1) = gamma;  end;
gamma = temp;

temp = zeros(M,L);
if (m > 0) temp(keep_list,:) = mu;  end;
mu = temp;
   
return;
