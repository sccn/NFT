% bem_load_session() - Loads the BEM session
%
% Usage:
%   >> session = bem_load_session(file);
%
% Inputs:
%   file - sesion name
%
% Outputs:
%   session - session structure.
%
% Author: Zeynep Akalin Acar, SCCN, 2007

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

function session = bem_load_session(file)

ssave = load(file, '-MAT');
mfile = [ssave.model_name, '.model'];
msave = load(mfile, '-MAT');
mesh = bem_load_mesh(msave.mesh_name);
model = bem_create_model(msave.name, mesh, msave.cond, msave.mod);
session = bem_create_session(ssave.name, model, ssave.Smatrix);    
session = bem_load_transfer_matrix(session, 'tmte');
if isfield(ssave,'sens')
    session.sens = ssave.sens;
end

            