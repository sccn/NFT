% nft_eeglab_dipplot() - Dipplot
%
% Usage:
%   >> nft_eeglab_dipplot(subject_name, session_name, of, dip)
%
% Inputs:
%   subject_name - subject name as entered in NFT main window
%   session_name - session name as entered in NFT main window
%   dip - dipole structure with fields posxyz, momxyz, rv 
%   of - output folder where all the outputs are saved
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

function nft_eeglab_dipplot(subject_name, session_name, of, dip)

current_folder = pwd;
cd(of)
a = dir([subject_name '_' session_name '_warping.mat']);
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
     mri_file = [eeglab_folder '/plugins/dipfit2.2/standard_BEM/standard_mri.mat'];
 end
 
 for i=1:size(dip,2);    dip(i).momxyz = dip(i).momxyz(:)'; end
 eeglab_dipplot(dip,'mri',mri_file,'projimg', 'off', 'projlines', 'off', 'axistight', 'on', 'cornermri','on','normlen','on');
 cd(current_folder)
