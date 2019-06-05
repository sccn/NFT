% bem_create_session() - Creates a session structure combining a model,
%       and sensor data. The session structure contains a model for a complete
%       head and electrode 'recording session'. Data recorded using the 
%       same set of sensor locations is considered to be in the same session.
%
% Usage:
%       >> session = bem_create_session(name, model, Smatrix);
%
% Inputs:
%   name    - session name, used as a base filename for matrices
%   model   - model structure obtained from bem_create_model().
%   Smatrix -  matrix that defines EEG electrodes in terms of 
%     the BEM mesh. Each electrode is a weighted sum of the nodes
%     of the element. The weights are determined by the element
%     shape functions. The format of the Smatrix is as follows:
%          [electrode_index node_index weight]
%     The rows of the matrix must be sorted by electrode_index;
%     there can be more than one row with a given electrode
%     index. The Smatrix can be constructed using one of these:
%     bem_smatrix_from_nodes() and bem_smatrix_from_points().
%
% Outputs:
%   session - session structure with the fields defined in the next section
%
% Session Structure:
%   name   - name of the session
%   model  - model structure
%   Smatrix - EEG Sensor information matrix.
%   num_electrodes - number of EEG sensors
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


function session = bem_create_session(name, model, Smatrix)

if ~isfield(model, 'mesh')
    error('BEM:bem_create_session:model','%s','Invalid model');
end

if ~isfield(model.mesh, 'num_nodes')
    error('BEM:bem_create_session:mesh','%s','Invalid mesh');
end

if size(Smatrix,2) ~= 3 || ~isnumeric(Smatrix) 
    error('BEM:bem_create_session:Smatrix','%s','Invalid Smatrix');
end

ne = max(Smatrix(:,1));

if min(Smatrix(:,1)) ~= 1 || length(unique(Smatrix(:,1))) ~= ne
    error('BEM:bem_create_session:Smatrix','%s','Invalid electrode index');
end

if min(Smatrix(:,2)) < 1 || max(Smatrix(:,2)) > model.mesh.num_nodes
   error('BEM:bem_create_session:Smatrix','%s','Invalid node index');
end 

session.name = name;
session.model = model;
session.Smatrix = Smatrix;
session.num_electrodes = ne;


