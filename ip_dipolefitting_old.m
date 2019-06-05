function [dipoles_str, session] = ip_dipolefitting(EEG, eloc, subject_name, session_name, comp_index, constr, warpback)

% Usage:
%   >> dipoles_str = ip_dipolefitting(EEG, sensor_file, subject_name,
%   session_name, comp_index, warpback);
%
% Inputs:
%   EEG - EEGLAB data structure
%   sensor_file - sensor file name
%   subject_name - subject name
%   session_name - session name
%   comp_index - component indices
%   warpback - warping structures
%
% Outputs:
%   dipoles_str - dipole structure as in dipfit

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


% session_name should be entered with its folder
% model_name should be entered with its folder
% if flag1 == 1 mri_based, no warping back

% for old type of sensor structure, this may be removed in future
%eloc = readlocs(sensor_file);

% load sensor index
a = dir([subject_name '_' session_name '.sensors']);
if size(a,1) == 0
    b = [subject_name '_' session_name '_sensorindex'];
    if size(b,1) == 1
        ind_n = [subject_name '_' session_name '_sensorindex'];
        sens_index = load(ind_n);
    else
        error('sensor file is not found!')
    end
else
    se = load([subject_name '_' session_name '.sensors'],'-mat');
    sens_index = se.ind;
end

% load LFM
load([session_name '_LFM']);


[A2, session, LFM2, ind_fp, ind_eeg, elocn] = eloc2eeglab_r(EEG, session_name, LFM, sens_index, eloc);
%elocn2 = elocn(ind_fp); % yeni warp edilmis electrod noktalari
% for the inverse problem solution
ss = load([subject_name '_sourcespace.dip']);

% check session type
if ~isfield(session, 'type') ||  ~strcmp(session.type, 'fem')
    % BEM
    if session.model.mod>0
        if ~isfield(session.model, 'iinv')
            session.model = bem_load_model_matrix(session.model,'iinv');
        end
        if ~isfield(session.model, 'dmt')
            session.model = bem_load_model_matrix(session.model,'dmt');
        end
    end
    %load vol 
    [vol, sens] = session2vol(session);
    sens.label = sens.label(ind_fp);
    sens.pnt = sens.pnt(ind_fp,:);
else
    %FEM
    vol = session.vol;
    sens = session.sens;
    metufem('setup',session.vol.mesh_name,'','')
    metufem('setrf',session.sens.rf)
end

%if isfield(session.model.mesh,'transform')
%    if length(session.model.mesh.transform) == 3
%        ss(:,1:3) = ss(:,1:3) - ones(size(ss,1),1) * session.model.mesh.transform;
%    end
%end
%tr=session.model.mesh.transform;
Ncomp = length(comp_index); % number of components
parfor compi = 1:Ncomp
    comp = comp_index(compi)
    Vdata = A2(:, comp);
    Vdata = Vdata - mean(Vdata);
    [pos_bin, griderror] = Grid_dipole2(LFM2, Vdata, ss);
    [Y1,I1] = sort(griderror);
    pos_bin_grid = pos_bin(I1(1),:);
    dip.pos = pos_bin_grid;
    %for bf = 1:1
        [dip_bin] = dipole_fit(dip, sens, vol, Vdata, 'constr', constr,'maxiter',500);
    %end
    for i=1:size(dip_bin.pos,1)
        mom(i,:) = dip_bin.mom((i-1)*3+1:i*3);
    end
    dip_bin.mom = mom;
    % compute fval
    if vol.type == 'metufem'
        
        for i=1:size(dip_bin.pos,1)
             p3(i,:) = [dip_bin.pos(i,:) dip_bin.mom(i,:)];
        end
        lf = metufem('pot', p3', 'interp');
    
    elseif vol.type == 'metubem'
          
          session = vol.session;
          for i=1:size(dip_bin.pos,1)
              p3(i,:) = [dip_bin.pos(i,:) dip_bin.mom(i,:)];
          end
          [lf, session] = bem_solve_lfm_eeg(session, p3);
          
    end
    lf = sum(lf,2);

    Vdata = Vdata-mean(Vdata);
    lf = lf-mean(lf);
    dif = Vdata(:) - lf(:);
    % relative residual variance
    fval = sum(dif(:).^2) / sum(Vdata(:).^2);
    if nargin > 6 % warp back dipole locations
        % warp back
        for i=1:size(dip_bin.pos,1)
            dip_bin.pos(i,:) = dip_bin.pos(i,:) + warp_lm(dip_bin.pos(i,:), warpback.A, warpback.W, warpback.LMd);
        end
    end
    
    % check this !!! (before or after warping???)
    %if isfield(session.model.mesh,'transform')
    %    if length(session.model.mesh.transform) == 3
    %        dip_bin.pos = dip_bin.pos + ones(size(dip_bin.pos,1),1)*session.model.mesh.transform;
    %    end
    %end
    
    dipoles_str(comp).posxyz = dip_bin.pos;
    dipoles_str(comp).momxyz = dip_bin.mom';
    dipoles_str(comp).rv = fval;
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,griderror] = Grid_dipole2(lf, Vdata, BrSS)
% compute grid error of individual dipole-sets in the LFM
% The LFM consists of sets of three orthagonal dipoles occupying
% the same location. 

% Zeynep Akalin Acar, 2008

nchans = length(Vdata);
ndip = size(lf,2) / 3;
m=size(BrSS,2)-6;
griderror=zeros(ndip,1);
pos=zeros(ndip,3);

s=1;
for d=1:ndip;
    i = [d d+ndip d+ndip+ndip];
    griderror(d) = sum(((eye(nchans)-lf(:,i)*pinv(lf(:,i)))*Vdata).^2);
    pos(d,:) = BrSS(d,1+m:3+m);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pos,griderror] = Grid_dipole2sym(lf, Vdata, BrSS)
% compute grid error of individual dipole-sets in the symmetric LFM
% The symmetric LFM consists of sets of six orthagonal dipoles
% for the same symmetric pair. 

% Zeynep Akalin Acar, 2011

nchans = length(Vdata);
ndip = size(lf,2) / 6;
m=size(BrSS,2)-6;
griderror=zeros(ndip,1);
pos=zeros(ndip,3);

s=1;
for d=1:ndip;
    i = [d d+ndip d+2*ndip d+3*ndip d+4*ndip d+5*ndip];
    griderror(d) = sum(((eye(nchans)-lf(:,i)*pinv(lf(:,i)))*Vdata).^2);
    pos(d,:) = BrSS(d,1+m:3+m);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [rw]=warp_lm(r,A,W,p)
% performs warp transformation with linear 3D RFB see Ermer's Thesis
rw = r * A(1:3,1:3) + repmat(A(4,:), size(r,1), 1);
for i = 1 : size(p,1)
    U = sqrt(sum((r - repmat(p(i,:), size(r,1),1)).^2, 2));  
    rw = rw + U * W(i,:);
end    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A2, session, LFM2, ind_fp, ind_eeg, elocn] = eloc2eeglab_r(EEG, session_name, LFM, elp_index, eloc);
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

if ~isfield(EEG.etc,'nft')
    EEG.etc.nft=[];
end
if ~isfield(EEG.etc.nft,'session')
    % load the session and the transfer matrix
    ssave = load([session_name '.session'], '-MAT');
    if isfield(ssave, 'model_name') 
        % BEM
        model = load_model([ssave.model_name, '.model']);
        session = bem_create_session(ssave.name, model, ssave.Smatrix);
        session = bem_load_transfer_matrix(session, 'tmte');
        % update the transfer matrix of the session wrt ind_fp
        session.tmte = session.tmte(:,ind_fp);
        session.num_electrodes = length(ind_fp);
        if session.model.mod>0
            session.model = bem_load_model_matrix(session.model,'iinv');
            session.model = bem_load_model_matrix(session.model,'dmt');
        end
   else
        session = ssave.session;
        session.sens.rf = session.sens.rf(:,ind_fp);
        session.sens.pnt = session.sens.pnt(ind_fp,:);    
    end
else
    session = EEG.etc.nft.session;
end




% update LFM wrt. ind_fp
LFM2 = LFM(ind_fp,:);

% update EEG.icawinv wrt. ind_eeg

A2 = A(ind_eeg,:);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function model = load_model(file)
msave = load(file, '-MAT');
%msave.mesh_name='pb03'
mesh = bem_load_mesh(msave.mesh_name);
model = bem_create_model(msave.name, mesh, msave.cond, msave.mod);
session = [];