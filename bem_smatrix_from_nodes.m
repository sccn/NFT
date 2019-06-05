% bem_smatrix_from_nodes() - Generates Smatrix from nodes of a mesh. 
%       Smatrix is the sensor information matrix used in
%       bem_create_session(). See bem_create_session() for more information.
%
% Usage:
%   >> Smatrix = bem_smatrix_from_nodes(mesh, nodes);
%
% Inputs:
%   mesh - mesh structure
%   nodes - vector of node indices.
% 
% Outputs:
%   Smatrix - defines one electrode per node.
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

function Smatrix = bem_smatrix_from_nodes(mesh, nodes)

if ~isempty(find(isfield(mesh, {'num_nodes'}) == 0,1))
    error('BEM:bem_smatrix_from_nodes:mesh','%s','Invalid mesh');
end
if ~isvector(nodes) || min(nodes)<1 || max(nodes)>mesh.num_nodes
    error('BEM:bem_smatrix_from_nodes:nodes','%s','Invalid nodes vector');
end

Smatrix(:,1) = (1:length(nodes))';
Smatrix(:,2) = nodes(:);
Smatrix(:,3) = 1;



    