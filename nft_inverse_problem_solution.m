% nft_inverse_problem_solution() - inverse problem solution
%
% Usage:
%   >> dip1 = nft_inverse_problem_solution(subject_name, session_name, of, EEG, comp_index, plotting, elec_file)
%
% Inputs:
%   subject_name - subject name as entered in NFT main window
%   session_name - session name as entered in NFT main window
%   of - output folder where all the outputs are saved
%   EEG - EEG data structure
%   comp_index - component index for dipole source localization
%   plotting - flag for plotting dipoles
%   elec_file - [Optional] electrode file
%
% Optional keywords:
%
%   sensor_name:  sensor name (default: [subject_name '_' session_name '.sensors'])
%   ss_name : sourcespace name (default: [subject_name '_sourcespace.dip'])
%
% Outputs:
%   dip1 - dipole structure with fields posxyz, momxyz, rv 
%
%
% Author: Zeynep Akalin Acar, SCCN, 2012

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

function dip1 = nft_inverse_problem_solution(subject_name, session_name, of, EEG, comp_index, plotting, elec_file,varargin)
current_folder = pwd;

cd(of)
sensor_name = [subject_name '_' session_name '.sensors'];
ss_name = [subject_name '_sourcespace.dip'];


for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};

      if ~isstr(Keyword)
         fprintf('keywords must be strings')
         return
      end

      if strcmp(Keyword,'sensor_name')
         if ~isstr(Value)
            fprintf('sensor_name must be a string');
            return
         else
             sensor_name = Value;
         end
      elseif strcmp(Keyword,'ss_name')
         if ~isstr(Value)
            fprintf('ss_name must be a string');
            return
         else
             ss_name = Value;
         end
      
      end
end

a = dir(sensor_name);
if size(a,1) > 0
    se = load(sensor_name,'-mat');
    eloc = se.eloc;
end

a = dir(['ori_sen_loc.mat']);
if size(a,1) > 0
    load ori_sen_loc
    elocfn = sens_fn;
    eloc = readlocs(sens_fn);
end

if nargin>6
    eloc = readlocs(elec_file);   % subject's electrode locations
    sens_fn = elec_file;
    elocfn = elec_file;
end

constr = [];

v = evalin('base','EEG');
if isfield(v.etc,'nft')
    dip1 = v.etc.nft.model;
else
    dip1 = [];
end

% check if mr-based realistic or warped mni
%a = dir(['Scalp.smf']);
a = dir([subject_name '_' session_name '_warping.mat']);
if size(a,1) == 0
    handles.mri_based = 1;
    [dip, session] = ip_dipolefitting(EEG, eloc, subject_name, session_name, comp_index, constr, [], 'sensor_name', sensor_name,'ss_name', ss_name);
else
    handles.mri_based = 0;
    fw = [subject_name '_' session_name '_warping'];
    load(fw)
    [dip, session] = ip_dipolefitting(EEG, eloc, subject_name, session_name, comp_index, constr, warping_param.back, 'sensor_name', sensor_name,'ss_name', ss_name);
end

EEG.etc.nft.session = session;

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

dipoles_str = dip1;



if plotting
    a = dir([subject_name '_' session_name '_warping.mat'])
    if size(a,1) > 0
        mri_based = 0;
    else
         mri_based = 1;
    end

    % mri-based realistic
    if mri_based == 1
        mri_file = [subject_name '_mri'];
    else
        % warped mni
        eeglab_folder = dirname(which('eeglab'));
        mri_file = [eeglab_folder '/plugins/dipfit/standard_BEM/standard_mri.mat'];
    end
    dip = dipoles_str;
    for i=1:size(dip,2);    dip(i).momxyz = dip(i).momxyz(:)'; end
    eeglab_dipplot(dip,'mri',mri_file,'projimg', 'off', 'projlines', 'off', 'axistight', 'on', 'cornermri','on','normlen','on');
end

cd(current_folder)


