function varargout = Warping_mesh(varargin)
% WARPING_MESH M-file for Warping_mesh.fig
%      WARPING_MESH, by itself, creates a new WARPING_MESH or raises the existing
%      singleton*.
%
%      H = WARPING_MESH returns the handle to a new WARPING_MESH or the handle to
%      the existing singleton*.
%
%      WARPING_MESH('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in WARPING_MESH.M with the given input arguments.
%
%      WARPING_MESH('Property','Value',...) creates a new WARPING_MESH or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Warping_mesh_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Warping_mesh_OpeningFcn via varargin.
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


% Edit the above text to modify the response to help Warping_mesh

% Last Modified by GUIDE v2.5 07-Jan-2011 16:03:50

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @Warping_mesh_OpeningFcn, ...
                   'gui_OutputFcn',  @Warping_mesh_OutputFcn, ...
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


% --- Executes just before Warping_mesh is made visible.
function Warping_mesh_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Warping_mesh (see VARARGIN)

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

if isfield(handles,'MeshFolder')
    set(handles.text5, 'String', handles.MeshFolder);
end

% Choose default command line output for Warping_mesh
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = Warping_mesh_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in warpingpushbutton.
function warpingpushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to warpingpushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

of = handles.MeshFolder;

lof = length(of);
if of(lof) ~= '/'
    of(lof+1) = '/';
end

[Ptm, ind, Cscalp_w,Cskull_w,CCSF_w,Cbrain_w,W,A,e,LMm2,back,Escalp,Eskull,ECSF,Ebrain] = warping_main_function(of,handles.MNImesh.Cscalp, ...
    handles.MNImesh.Escalp, handles.MNImesh.Eskull, handles.MNImesh.ECSF, handles.MNImesh.Ebrain,handles.MNImesh.Cskull, handles.MNImesh.CCSF, handles.MNImesh.Cbrain, ...
    handles.MNILandmarks, handles.MNIFiducials, ...
    handles.fiducials, handles.electrodes, handles.index_kdm);

handles.MNImesh.Escalp = Escalp;
handles.MNImesh.Eskull = Eskull;
handles.MNImesh.ECSF = ECSF;
handles.MNImesh.Ebrain = Ebrain;

handles.warpedMNImesh.Cscalp = Cscalp_w;
handles.warpedMNImesh.Cskull = Cskull_w;
handles.warpedMNImesh.CCSF = CCSF_w;
handles.warpedMNImesh.Cbrain = Cbrain_w;
handles.fitelectrodes = Ptm;
handles.chosenindices = ind;
guidata(handles.figure1, handles);

handles.warping.back=back;
handles.warping.forward.W = W;
handles.warping.forward.A = A;
handles.warping.forward.e = e;
handles.warping.forward.LMm2 = LMm2;

axes(handles.axes2);

% convert to linear mesh if it is quadratic
if size(handles.MNImesh.Escalp,2) == 7
    Elem(:,1:2) = handles.MNImesh.Escalp(:,1:2);
    Elem(:,3) = handles.MNImesh.Escalp(:,4);
    Elem(:,4) = handles.MNImesh.Escalp(:,6);
    np = max(max(Elem(:,2:4)));
    Coord = Cscalp_w(1:np,:);
else
    Elem = handles.MNImesh.Escalp; Coord = Cscalp_w;
end

eeglab_plotmesh(Elem(:,2:4), Coord(:,2:4),[],1); hold; axis on;
plot3(Ptm(:,1), Ptm(:,2), Ptm(:,3), 'b.')
axis([-100 100 -200 100 -100 150])
view(165, 10)

% save the electrode locations, and ind
p = handles.MeshFolder;

lof = length(p);
if p(lof) ~= '/'
    p(lof+1) = '/';
end

% save warped sensors and index
if isfield(handles,'arg_subject')
    fsubj = [handles.arg_subject];
else
    fsubj = 'tempSubj';
end
if isfield(handles,'arg_session')
    fses = [handles.arg_session];
else
    fses = 'tempSes';
end
f = [fsubj '_' fses];
warped_sensors = Ptm;
%save([p f '_headsensors.sens'], 'warped_sensors', '-ascii');
%clear warped_sensors;
%save([p f '_sensorindex'], 'ind', '-ascii'); % save the index, to use in IP

ssave.fn = handles.elocfn;
ssave.eloc = handles.eloc;
ssave.pnt = warped_sensors;
ssave.ind = ind;
save([p f '.sensors'], '-STRUCT', 'ssave')

nl = str2num(get(handles.edit_nl, 'String'));
% save the mesh
[Coord, Elem] = utilbem_add_mesh(Cscalp_w, handles.MNImesh.Escalp, Cskull_w, handles.MNImesh.Eskull);
[Coord, Elem] = utilbem_add_mesh(Coord, Elem, CCSF_w, handles.MNImesh.ECSF);
if nl==4
    [Coord, Elem] = utilbem_add_mesh(Coord, Elem, Cbrain_w, handles.MNImesh.Ebrain);
end
    
save([p fsubj '.bec'], 'Coord', '-ascii');


Info(1,1) = nl;
if nl==3
    Info(1,2) = size(handles.MNImesh.Escalp,1)+size(handles.MNImesh.Eskull,1)+size(handles.MNImesh.ECSF,1);
    Info(1,3) = size(Cscalp_w,1)+size(Cskull_w,1)+size(CCSF_w,1);
elseif nl==4
    Info(1,2) = size(handles.MNImesh.Escalp,1)+size(handles.MNImesh.Eskull,1)+size(handles.MNImesh.ECSF,1)+size(handles.MNImesh.Ebrain,1);
    Info(1,3) = size(Cscalp_w,1)+size(Cskull_w,1)+size(CCSF_w,1)+size(Cbrain_w,1);
    Info(5,1) = 4;
    Info(5,2) = size(handles.MNImesh.Ebrain,1);
    Info(5,3:4) = [4 3];
end
Info(1,4) = size(handles.MNImesh.Escalp,2)-1; % number of nodes per element
Info(2,1) = 1;
Info(2,2) = size(handles.MNImesh.Escalp,1);
Info(2,3:4) = [1 0];
Info(3,1) = 2;
Info(3,2) = size(handles.MNImesh.Eskull,1);
Info(3,3:4) = [2 1];
Info(4,1) = 3;
Info(4,2) = size(handles.MNImesh.ECSF,1);
Info(4,3:4) = [3 2];



fid = fopen([p fsubj '.bei'], 'w');
fprintf(fid, '%d %d %d %d\r\n', Info');
fclose(fid);

if size(Elem,2) == 4
    fid = fopen([p fsubj '.bee'],'w');
    fprintf(fid, '%d %d %d %d\r\n', Elem');
    fclose(fid);
elseif size(Elem,2) ==7
    fid = fopen([p fsubj '.bee'],'w');
    fprintf(fid, '%d %d %d %d %d %d %d\r\n', Elem');
    fclose(fid);
end


% save warping parameters
warping_param = handles.warping;

save([p f '_warping.mat'], 'warping_param'); % save the index, to use in IP

load Warping_so_MNIdata4L

rw = warp_lm(so_MNIdataconv, A, W, LMm2) + so_MNIdataconv;
Ns = size(rw, 1);
so = zeros(3*Ns, 6);
so(1:Ns, 1:3) = rw;
so(1:Ns, 4) = 1;
so(1+Ns:2*Ns, 1:3) = rw;
so(1+Ns:2*Ns, 5) = 1;
so(1+Ns*2:3*Ns, 1:3) = rw;
so(1+Ns*2:3*Ns, 6) = 1;


% save source space
save([p fsubj '_sourcespace.dip'], 'so', '-ascii'); 


set(handles.warpingpushbutton, 'String', 'Mesh Warped!');
set(handles.warpingpushbutton, 'Enable', 'off');

set(handles.pushbuttonFEM, 'Enable', 'on');
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




% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.MeshFolder = uigetdir;
set(handles.text5, 'String', handles.MeshFolder);
% Update handles structure
guidata(handles.figure1, handles);




% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[file, path] = uigetfile('*.*');  % changed 01-06-2011
handles.elocfn = [path file];
if ~isequal(file, 0) && length(file) > 1
    eloc = readlocs([path file]);   % subject's electrode locations
    handles.eloc = eloc;
    if ~strcmp(eloc(1).type,'FID') | ~strcmp(eloc(2).type,'FID') | ~strcmp(eloc(3).type,'FID')
        error('Electrode file does not contain fiducials! Co-registration is done using the fiducials!')
    end
    if ~strcmp(eloc(1).labels,'Nz')
        if ~strcmp(eloc(1).labels,'fidt9')
            warning('Fiducials are assumed to be in this order: [Nz LPA, RPA].')
            h = msgbox('Fiducials are assumed to be in this order: [Nz LPA, RPA]','','warn')
        else   % for .sfp files
            eloc2=eloc;
            eloc2(1)=eloc(2);
            eloc2(2)=eloc(1);
            eloc=eloc2;
        end
        
    end
    sens_fn = [path file];
    for i = 1:length(eloc); elo(i,:) = [eloc(i).X eloc(i).Y eloc(i).Z]; end
    [d, elo] = warping_distafterwarping([0 0 0 0 0 90], elo, elo); % arrange orientation ??? check!
    ne = size(elo,1);
    
    p = handles.MeshFolder; % save the files in mesh folder
    lof = length(p);
    if p(lof) ~= '/';   p(lof+1) = '/'; end;
    save([p 'ori_sen_loc'], 'sens_fn'); % save the location of original sensors
end

handles.input_electrodes = elo;
% set filename removing .elp
handles.filename = file(1:length(file)-4);
handles.filepath = path;

set(handles.warpingpushbutton, 'Enable', 'on');
set(handles.warpingpushbutton, 'String', 'Start Warping');
set(handles.text6,'String',[path file])

updateMeshPlot(handles);

% Update the mesh and plot it
function updateMeshPlot(handles)

nl = str2num(get(handles.edit_nl, 'String'));

load Warping_MNIdata4L

handles.MNImesh.Cscalp = Cscalp;
handles.MNImesh.Escalp = Escalp;
handles.MNImesh.Cskull = Cskull;
handles.MNImesh.Eskull = Eskull;
handles.MNImesh.CCSF = CCSF;
handles.MNImesh.ECSF = ECSF;
handles.MNImesh.Cbrain = Cbrain;
handles.MNImesh.Ebrain = Ebrain;
handles.MNIFiducials = Fm;
handles.MNILandmarks = LMm;

elo = handles.input_electrodes;
a1 = max(elo) - min(elo);
a2 = max(Cscalp(:,2:4)) - min(Cscalp(:,2:4));
rat = mean(a2./a1);
% make the same scale with the mesh
if rat>500
    elo = elo * 1000; 
elseif rat>50
    elo = elo * 100; 
elseif rat>5
    elo = elo * 10;
end
handles.input_electrodes = elo;



cla(handles.axes1);
axes(handles.axes1);
hold off

% convert to linear mesh if it is quadratic
if size(Escalp,2) == 7
    Elem(:,1:2) = Escalp(:,1:2);
    Elem(:,3) = Escalp(:,4);
    Elem(:,4) = Escalp(:,6);
    np = max(max(Elem(:,2:4)));
    Coord = Cscalp(1:np,:);
else
    Elem = Escalp; Coord = Cscalp;
end

eeglab_plotmesh(Elem(:,2:4), Coord(:,2:4),[],1); hold on; axis on;
axis([-100 100 -200 100 -100 150])
view(165, 10)

if isfield(handles, 'input_electrodes')
    elo = handles.input_electrodes;
    [pos, Fd] = initial_registration(elo, elo(1:3,:), Cscalp, Fm);
    
    % find the index of the electrodes that are close to the scalp
    [elox, dm] = warping_distmeshafterwarping([0 0 0 0 0 0], pos, Cscalp, Escalp);
    mdm = median(dm); sdm = std(dm);
    kdm = find((dm < 2*mdm)); % kdm gives the index of the electrodes close to the scalp
    index = [1:length(pos)]; rejected = setdiff(index,kdm)

    
    handles.electrodes = pos;
    handles.fiducials = Fd;
    handles.index_kdm = kdm;
    pos = handles.electrodes;
    plot3(pos(:,1), pos(:,2), pos(:,3), 'b.')
    plot3(pos(rejected,:),pos(rejected,2),pos(rejected,3),'ro')
end

% Update handles structure
guidata(handles.figure1, handles);


% --------------------------------------------------------------------
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to uipanel2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

updateMeshPlot(handles);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [pos, Fd] = initial_registration(elo, F, Cscalp, Fm);

ne = size(elo,1);
[P1e, P2e] = find_new_points_for_reg(elo, F);
[P1m, P2m] = find_new_points_for_reg(Cscalp(:,2:4), Fm);

% find coarse scaling (x) using F2-F3
sx = (F(2,:)-F(3,:))/(Fm(2,:)-Fm(3,:)); sx2=abs(1-sx);
% find coarse scaling (y) using P1-F1
sy = (P1e-F(1,:))/(P1m-Fm(1,:)); sy2=abs(1-sy);
% find coarse scaling (z) using P1-P2
sz = (P1e-P2e)/(P1m-P2m); sz2=abs(1-sz);
% find the one closest to 1, smallest scaling factor
[k,l] = min([sx2 sy2 sz2]);
a = [sx sy sz]; min_sc=a(l);
% after scaling
elo2 = elo;
elo2(:,1) = elo2(:,1) / min_sc;
elo2(:,2) = elo2(:,2) / min_sc;
elo2(:,3) = elo2(:,3) / min_sc;
F2 = elo2(1:3,:);

[P1e2, P2e2] = find_new_points_for_reg(elo2, F2);

% find coarse translation using P1e, P1m
tr = P1e2 - P1m;
elo3 = elo2 - ones(length(elo2),1) * tr;
F3 = F2 - ones(3,1)*tr;

[P1e3, P2e3] = find_new_points_for_reg(elo3, F3);

% find the rotation using fiducials and P2
options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun',1e-6);
Xo = [0 0 0];
X = fminsearch(@(X) funrstPP_rot(X, [F3; P2e3], [Fm;P2m]), Xo, options);
X2=[0 0 0 X];
% find the rotated digitizer locations
[d, elo4] = warping_distafterwarping(X2, elo3, elo3);

[P1x, P2x] = find_new_points_for_reg(elo4, elo4(1:3,:));
%pos = elo4(4:ne,:); % don't take fiducials
pos = elo4; % save with fiducials
Fd = elo4(1:3,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = funrstPP_rot(X, F, Fe);
% F is the point set of the digitizer
% d is the distance between translated and rotated F and Fe
% X is the vector of rotation parameters

alpx = X(1)*pi/180;
alpy = X(2)*pi/180;
alpz = X(3)*pi/180;

x = F(:,1);
y = F(:,2);
z = F(:,3);

% rotation around x-axis
x1 = x;
y1 = y * cos(alpx) - z * sin(alpx);
z1 = y * sin(alpx) + z * cos(alpx);

% rotation around y-axis
x2 = z1 * sin(alpy) + x1 * cos(alpy);
y2 = y1;
z2 = z1 * cos(alpy) - x1 * sin(alpy);

% rotation around z-axis
x3 = x2 * cos(alpz) - y2 * sin(alpz);
y3 = x2 * sin(alpz) + y2 * cos(alpz);
z3 = z2;

N = size(Fe,1);

Ma = Fe - [x3 y3 z3];
d = sum(sqrt(sum(Ma.*Ma,2)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [P1, P2] = find_new_points_for_reg(elo,F);
% elo is the digitizer location (ne x 3)
% F is the fiducials (3 x 3)
%     1st row nasion
%     2nd row LPA
%     3rd row RPA
% P1 is the mean for the ear fiducials
% P2 is the upper point of the line that is perpendicular to the F1-F2-F3 plane
%       that intersects the digitizer locations
 
ne = length(elo); % number of electrodes

F1 = F(1,:); % nasion
F2 = F(2,:); % LPA
F3 = F(3,:); % RPA


% P1 is the mean point of ear fiducials
P1 = (F2 + F3) / 2;

% plane equation for F1-F2-F3
AB = F2 - F1;
AC = F3 - F1;
n = cross(AB,AC); nz=n(3);
% find the line equation perperdicular to F1-F2-F3 plane r(t)
max_d = max(elo(:,3))/nz*2;
min_d = min(elo(:,3))/nz;
incr = (max(elo(:,3))-min(elo(:,3)))/nz/100;
t = min_d:incr:max_d;
r = ones(length(t),1)*P1 + t'*n; %in terms of t

% find the closest electrode point to r
for i=1:ne
    p1 = elo(i,:);
    M = r - ones(length(t),1)*p1;
    M = sqrt(sum(M.*M,2));
    [k,l] = min(M);
    dis(i,1) = k; % minimum distance
    dis(i,2) = l; % index of r
end

[k,l] = min(dis(:,1));

rm = dis(l,2); %the index of closest point on r
P2 = r(rm,:);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rw]=warp_lm(r,A,W,p)
rw = r * A(1:3,1:3) + repmat(A(4,:), size(r,1), 1);
for i = 1 : size(p,1)
    U = sqrt(sum((r - repmat(p(i,:), size(r,1),1)).^2, 2));  
    rw = rw + U * W(i,:);
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

if isfield(handles,'arg_subject')
    mesh_name = [handles.arg_subject];
else
    mesh_name = 'tempSubj';
end


fn = [of mesh_name '.bei'];

if exist(fn) == 0
    error('BEM mesh not found!')
end

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
if status ~= 0; error('Warping_mesh:system','Failed to execute: %s',result); end

% convert into metu-fem mesh
nl = mesh.num_boundaries;
if nl == 4
    cnd_str = '1=C1 2=C2 3=C3 4=C4';
elseif nl == 3
    cnd_str = '1=C1 2=C2 3=C3';
end
fn = [of mesh_name '.1'];
a = sprintf('"%s" "%s" %s', conf.tetgen2msh, fn, cnd_str);
[status, result] = system(a);
if status ~= 0; error('Warping_mesh:system','Failed to execute: %s',result); end


set(handles.pushbuttonFEM, 'String', 'FEM mesh generated!');
set(handles.pushbuttonFEM, 'Enable', 'off');

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
