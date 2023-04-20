% nft_get_config() - Returns configuration information program path names,
%       defaults etc. Other BEM functions call this m-file
%       to read their configurations.
%
% Outputs:
%   conf - config structure.
%
% Config structures:
%   bem_matrix_program - name of the program that creates
%                        BEM matrices.
%   asc - name of the adaptive skeleton climbing program.
%   qslim - name of the coarsening program.
%   showmesh - name of the correction and smoothing program.
%
% Author: Zeynep Akalin Acar, SCCN, 2009

% Copyright (C) 2009 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
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


function conf = nft_get_config

% Set toolbox_dir to the installation path of the toolbox or leave as is for autodetect:
%mfiledir = dirname(mfilename('fullpath'));
mfiledir = mfilename('fullpath');
mfiledir = mfiledir(1:max(find(mfiledir == filesep) - 1));
% Strip out the last component which should be 'mfiles'

     %conf.nft_dir = mfiledir(1:max(find(mfiledir == filesep) - 1)); % comment
     %out for eeglab

% Add bin
     %bindir = [conf.nft_dir filesep 'bin' filesep]; % for eeglab
     bindir  = [mfiledir filesep];   % for eeglab
if ismac
    sufx = '.osx';
elseif isunix
    sufx = '';
elseif ispc
    sufx = '.exe';
else
    error('Platform not supported')
end

% Settings for 32 and 64 bit Linux
conf.bem_matrix_program = [bindir 'bem_matrix' sufx];
conf.asc                = [bindir 'asc1' sufx];
conf.qslim              = [bindir 'qslim' sufx];
conf.metufem            = [bindir 'forward' sufx];
conf.quad               = [bindir 'quadmesh' sufx];
conf.tetgen             = [bindir 'tetgen' sufx];
conf.tetgen2msh         = [bindir 'tetgen2msh.sh']; % ? windows?
conf.showmesh           = [bindir 'procmesh' sufx];
conf.lin2quad           = [bindir 'lin2quad' sufx];
conf.showmesh3          = [bindir 'Showmesh' sufx];

conf.metufem = '/home/zeynep/Programs/metu_fem-0.5/forward';
conf.metufem = '/home/zeynep/Programs/metu_fem-0.4/forward';
%conf.metufem = '/home/zeynep/Programs/metu_fem-0.5b/forward';

%conf.metufem = '/home/zeynep/Programs/metu_fem-0.5x/forward';

%conf.metufem = '/home/zeynep/Programs/metu_fem-0.6/forward.sh';
