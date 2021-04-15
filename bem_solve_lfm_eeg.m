% bem_solve_lfm_eeg() - Computes the LFM arising from given dipoles
% at the sensor locations defined by the session. 
% 
% Usage:
%   >> [pot, session] = bem_solve_lfm_eeg(session, dipoles);
%
% Inputs:
%   session - session structure defining the model and the sensors
%   dipoles - dipole matrix. Each row defines a dipole as follows:
%            [x y z px py pz]
%
% Outputs:
%   pot - the LFM at the sensors produced by activation of each dipole, 
%   session - updated session structure, in case any matrices were
%             loaded that modified the underlying model.
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

function [pot, session] = bem_solve_lfm_eeg(session, dipoles)

%if isfield(session.model.mesh,'transform')
%    if length(session.model.mesh.transform) == 3
%        dipoles(:,1:3) = dipoles(:,1:3) + ones(size(dipoles,1),1) * session.model.mesh.transform;
%    end
%end

% check session
if ~isempty(find(isfield(session, {'name', 'model'}) == 0,1))
    error('BEM:bem_solve_dipoles_eeg:session','%s','Invalid session');
end
model = session.model;
if ~isempty(find(isfield(model, {'name', 'mesh', 'node_cond', ...
        'cond','mod'}) == 0,1))
    error('BEM:bem_solve_dipoles_eeg:model','%s','Invalid model');
end
mesh = model.mesh;
if ~isempty(find(isfield(mesh, {'name','bnd','coord'}) == 0,1))
    error('BEM:bem_solve_dipoles_eeg:mesh','%s','Invalid mesh');
end

rhs = zeros(mesh.num_nodes, size(dipoles,1));

if model.mod < 1
    parfor k = 1:size(dipoles,1)
        rhs(:,k) = utilbem_pot_unbound(mesh.coord, model.node_cond, dipoles(k,:));
    end
else
    if ~isempty(find(isfield(model, {'ind_mod', 'ind_imesh', 'ind_imesh_mod'}) == 0,1))
        error('BEM:bem_solve_dipoles_eeg:model','%s','Invalid model');
    end
    if ~isfield(model,'iinv')
        %session.model = bem_load_model_matrix(model, 'iinv');
        model = bem_load_model_matrix(model, 'iinv');   % changed 02.29
    end
    if ~isfield(model,'dmt')
        session.model = bem_load_model_matrix(model, 'dmt');
    end
    model = session.model;
    s2 = model.cond(mesh.bnd(model.mod,3)); % outer
    s3 = model.cond(mesh.bnd(model.mod,2)); % inner
    parfor k = 1:size(dipoles,1)
        rhs(:,k) = utilbem_multilayer_rhs(mesh.coord, model.node_cond, model.ind_mod, ...
            model.ind_imesh, model.ind_imesh_mod, model.iinv, model.dmt, s2, s3, dipoles(k,:));
    end
end
if ~isfield(session, 'tmte')
    session = bem_load_transfer_matrix(session, 'tmte');
end
pot = session.tmte' * rhs;
