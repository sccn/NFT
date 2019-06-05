% metufem_calcsens() - Calculates sensitivity vectors.
%
% Usage:
%   >> sens = metufem_calcsens(vol, sens, of, cond)
%
% Inputs:
%   vol - volume 
%   sens - sensor structure
%   of  - output folder
%   cond - conductivities
%   
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


function [pot0] = metufem_calcpot(vol, of, cond, dipoles)
% conductivity values cond=[0.33 0.0042 1.79 0.33]
nl = length(cond);
if ~ft_voltype(vol, 'metufem')
    error('Volume type is not metufem')
end

if ~isfield(vol, 'mesh_name')
    error('FEM mesh name is not set in volume');
end

% XXX for testing  only. Fix before release
conf.metufem = '/home/zeynep/Programs/metu_fem-0.5b/forward';
conf.coordmap = '/home/zeynep/Programs/fem_util/coordmap/coordmap';


lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end

curdir = pwd;
cd(of);

c = '';
for i=1:nl
    if ~isempty(c)
        c = [c ','];
    end
    c = sprintf('%s%d=%0.5g', c, i, cond(i));
end

% convert dipoles
dipole_in = 'dipoles';
dipole_file = 'dipoles.fem';
num_dip = size(dipoles, 1);
if size(dipoles,2) == 6
    % dipole solutions
    nd = size(dipoles,1);
    dip(:,2:7) = dipoles;
    dip(:,1) = (1:num_dip)';

    fid = fopen(dipole_in, 'w');
    fprintf(fid,'%d\n', num_dip);
    fprintf(fid,'%d %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f\r\n', dip');
    fclose(fid);
elseif size(dipoles,2) == 7
    % patch solutions
    nd = length(unique(dipoles(:,7)));
    dip(:,2:8) = dipoles;
    dip(:,1) = (1:num_dip)';

    fid = fopen(dipole_in, 'w');
    fprintf(fid,'%d\n', num_dip);
    fprintf(fid,'%d %5.15f %5.15f %5.15f %5.15f %5.15f %5.15f %d\r\n', dip');
    fclose(fid);
end

a= sprintf('%s -d -m %s%s -i %s%s -o %s%s', conf.coordmap, of, vol.mesh_name, of, dipole_in, of, dipole_file);

[status, result] = system(a);
if status ~= 0;
    warning('coordmap:system', 'Warning:: %s\nOutput:\n%s\n', a, result);
end

% call 'forward' with vol.mesh_name and sens_name
a = sprintf('"%s" -m "%s%s" -d "%s%s" -c "%s"', conf.metufem, of, vol.mesh_name, of, dipole_file, c);

[status, result] = system(a);
if status ~= 0;
    warning('MetuFEM:system', 'Warning:: %s\nOutput:\n%s\n', a, result);
end

% read pot files to matrix
pot = load('x000.pot');
pot0 = zeros(length(pot), num_dip);
pot0(:,1) = pot;

for i = 2:num_dip
    fn = sprintf('x%03d.pot', i - 1);
    pot = load(fn);
    pot0(:,i) = pot;
end
