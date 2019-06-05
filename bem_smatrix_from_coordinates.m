% bem_smatrix_from_coordinates() - Generates Smatrix from sensor coordinates. 
%       Smatrix is the sensor information matrix used in
%       bem_create_session(). See bem_create_session() for more information.
%
% Usage:
%   >> Smatrix = bem_smatrix_from_coordinates(mesh, coords);
%
% Inputs:
%   mesh - mesh structure
%   coords - matrix of coordinate locations (Nx3).
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

function Smatrix = bem_smatrix_from_coordinates(mesh, coords)

if ~isempty(find(isfield(mesh, {'num_nodes'}) == 0,1))
    error('BEM:bem_smatrix_from_coords:mesh','%s','Invalid mesh');
end
if size(coords,2) ~= 3
    error('BEM:bem_smatrix_from_coords:coords','%s','Invalid coord matrix');
end
nne = mesh.num_node_elem;
Coord = zeros(mesh.num_nodes,4);
Elem = zeros(mesh.num_elements,nne+1);
Coord(:,2:4) = mesh.coord;
Elem(:,2:nne+1) = mesh.elem;
Coord(:,1) = [1:mesh.num_nodes]';
Elem(:,1) = [1:mesh.num_elements]';

Smatrix = utilbem_sens_mat(Coord, Elem, coords);




    