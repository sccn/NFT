% nft_dsl_source_visualization() - Distributed source localization
% (dsl) source visualization
%
% Usage:
%   >> nft_dsl_source_visualization(sourceJ, comp_plot, of)
%
% Inputs:
%   sourceJ - source localization obtained in
%   nft_dsl_inverse_problem_solution
%   comp_plot - component number to plot
%   of - output folder where all the outputs are saved
%
%
% Author: Zeynep Akalin Acar, SCCN, 2021

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

function nft_dsl_source_visualization(sourceJ, comp_plot, of)


conf = nft_get_config;

curr_dir = pwd;
cd(of)

pot = sourceJ(:,comp_plot);
save pot pot -ascii

f = fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'nfield load %spot',of);
%fprintf(f, 'quit\n');
fclose(f);

a = sprintf('"%s" -c "%sStepSc.txt" FSss.smf', conf.showmesh3, of);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

cd(curr_dir)

