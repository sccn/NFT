% nft_fem_forward_problem_solution() - FEM forward problem solution
%
% Usage:
%   >> nft_fem_forward_problem_solution(subject_name, session_name, of)
%
% Inputs:
%   subject_name - subject name as entered in NFT main window
%   session_name - session name as entered in NFT main window
%   of - output folder where all the outputs are saved
%
% Optional keywords:
%
%   cond : conductivity vector (default = [0.33 0.0132 1.79 0.33])
%   mesh_name : mesh name that will be loaded (default: [subject_name '.1.msh'])
%   sensor_name:  sensor name (default: [subject_name '_' session_name '.sensors'])
%   ss_name : sourcespace name (default: [subject_name '_sourcespace.dip'])
%   LFM_name :  LFM name (default: session_name_LFM)
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


function nft_fem_forward_problem_solution(subject_name, session_name, of, varargin)

% default conductivity values
cond = [0.33 0.0132 1.79 0.33];
% load mesh

mesh_name = [subject_name '.1.msh'];
sensor_name = [subject_name '_' session_name '.sensors'];
ss_name = [subject_name '_sourcespace.dip'];
LFM_name = [session_name '_LFM.mat'];


for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};

      if ~isstr(Keyword)
         fprintf('keywords must be strings')
         return
      end

      if strcmp(Keyword,'cond')
         if isstr(Value)
            fprintf('cond must be vector');
            return
         else
            cond = Value;
         end
      elseif strcmp(Keyword,'mesh_name')
         if ~isstr(Value)
            fprintf('mesh_name must be a string');
            return
         else
             mesh_name = Value;
         end
      elseif strcmp(Keyword,'sensor_name')
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
       elseif strcmp(Keyword,'LFM_name')
         if ~isstr(Value)
            fprintf('LFM_name must be a string');
            return
         else
             LFM_name = Value;
         end
      end
end
current_folder = pwd;
lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end
cd(of)

% load the FEM mesh
[Coordf, Elemf, Sig] = fem_load_mesh (mesh_name);
nl = length(unique(Sig));
sprintf('Number of layers in the model is: %d',nl)

if (nl == 3 && length(cond) == 4)
    cond = [cond(1) cond(2) cond(4)];
end


% set conductivity values
sens = load(sensor_name, '-mat'); % sensor locations
ss = load(ss_name); % sourcespace

vol = metufem_set_mesh(mesh_name);
sens.type = 'eeg';
sens = metufem_calcrf(vol, sens, of, cond);


session.name = session_name;
session.cond = cond;
session.sens = sens;
session.type = 'fem';
session.vol = metufem_set_mesh([of mesh_name]);

% save session
msave.session = session;
msave.mesh_name = mesh_name;
msave.mesh_path = of;
save([session.name '.session'], '-STRUCT', 'msave')


metufem('setup', mesh_name, 'sensors.dat', '');
metufem('setrf', session.sens.rf);

LFM = metufem('pot', ss','interp');


save(LFM_name,'LFM');
clear LFM;
cd(current_folder)



