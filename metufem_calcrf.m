% metufem_calcrf() - Calculates rf vectors.
%
% Usage:
%   >> sens = metufem_calcrf(vol, sens, of, cond)
%
% Inputs:
%   vol - volume 
%   sens - sensor structure
%   of  - output folder
%   cond - conductivities
%   
%
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


function sens = metufem_calcrf(vol, sens, of, cond)
% conductivity values cond=[0.33 0.0042 1.79 0.33]
nl = length(cond);
if ~ft_voltype(vol, 'metufem')
    error('Volume type is not metufem')
end

if ~isfield(vol, 'mesh_name')
    error('FEM mesh name is not set in volume');
end

if ~ft_senstype(sens, 'eeg')
    warning('Only EEG sensors are supported for now!');
end

conf = nft_get_config;

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end

curdir = pwd;
cd(of);

sens_name = 'sensors.dat';
num_sens = size(sens.pnt, 1);
sensors(:,2:4) = sens.pnt;
sensors(:,1) = (1:num_sens)';

fid = fopen(sens_name, 'w');
fprintf(fid,'%d\n', num_sens);
fprintf(fid,'%d %5.15f %5.15f %5.15f\r\n', sensors');
fclose(fid);

c = '';
for i=1:nl
    if ~isempty(c)
        c = [c ','];
    end
    c = sprintf('%s%d=%0.5g', c, i, cond(i));
end

% call 'forward' with vol.mesh_name and sens_name
a = sprintf('"%s" -m "%s%s" -s "%s%s" -c %s -rfpot', conf.metufem, of, vol.mesh_name, of, sens_name, c);
[status, result] = system(a);
if status ~= 0;
    warning('MetuFEM:system', 'Warning:: %s\nOutput:\n%s\n', a, result);
end

% read rf files to matrix
rf = load('rf000.pot');
tmte = zeros(length(rf), num_sens);
tmte(:,1) = rf;

for i = 2:num_sens
    fn = sprintf('rf%03d.pot', i - 1);
    rf = load(fn);
    tmte(:,i) = rf;
end

cd(curdir);

% set sens.rf
sens.rf = tmte;

%metufem('setup', vol.mesh_name, '', '');
%metufem('setrf', sens.rf);
