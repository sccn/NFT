% bem_create_model() - Creates a model structure combining a mesh,
%                      conductivity information and BEM parameters.
%
% Usage:
%   >> model = bem_create_model(name, mesh, cond, mod);
%
% Inputs:
%   name - model name, used as a base filename for matrices
%   mesh - mesh structure obtained from BEM_LOAD_MESH
%   cond - conductivity values for mesh tissue classes
%          vector, starting from first tissue type;
%          0 is reserved for air
%   mod - index of the modified boundary, for use with IPA.
%         If mod <= 0, IPA is not used.
%
% Outputs:
%   model - model structure with the following fields.
%
% Model Structure:
%   name - name of the model
%   mesh - mesh structure
%   cond - conductivity vector
%   mod  - modified boundary information
%   node_cond - average conductivity around a node
%
% Optional Fields:
%   ind_mod       - node indices of the modified boundary
%   ind_imesh     - node indices of the inner mesh
%   ind_imesh_mod - node indices of the modified boundary
%                   relative to the inner mesh
%  
% Author: Zeynep Akalin Acar, SCCN, 2007
%
% Notes: The optional fields in the model structure are only computed if
%        IPA is in use (mod > 0)

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

function model = bem_create_model(name, mesh, cond, mod)

if ~isempty(find(isfield(mesh, {'num_boundaries','coord','elem','bnd'}) == 0,1))
    error('BEM:bem_create_model:mesh','%s','Invalid mesh');
end

if ~isscalar(mod) || mod > mesh.num_boundaries
    error('BEM:bem_create_model:mod','%s','Invalid "mod" parameter');
end

if (~isvector(cond) || ~isnumeric(cond) || length(cond) ~= mesh.num_class)
    error('BEM:bem_create_model:cond','%s','Invalid "cond" vector');
end

model.name = name;
model.mesh = mesh;
model.cond = cond(:);
model.mod = mod;

[ncond, ncondE] = utilbem_compute_cond(mesh.coord, mesh.elem, mesh.bnd, cond);
model.node_cond = ncond;
model.node_cond_diff = ncondE;

if mod > 0
    [indMsh, indMod, indMshMod] = utilbem_compute_indices( ...
        mesh.coord, mesh.elem, mesh.bnd, mod);

    model.ind_imesh = indMsh;
    model.ind_mod = indMod;
    model.ind_imesh_mod = indMshMod;
end

% XXX TODO: check if IPA is better or not necessary using cond and mod
