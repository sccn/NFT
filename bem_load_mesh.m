% bem_load_mesh() - Loads a BEM mesh.
%
% Usage:
%   >> mesh = bem_load_mesh(name);
%
% Inputs:
%   name - mesh name excluding the extension.
%
% Outputs:
%   mesh - mesh structure.
%
% Mesh Structure:
%   name - mesh name
%   coord - node coordinates
%   elem - elements (connectivity information)
%   num_nodes - number of nodes
%   num_elements - number of elements
%   num_boundaries - number of boundaries
%   num_node_elem - number of nodes per element
%   num_class - number of tissue classes
%   bnd - boundary information array
%        [num_elements inner_tissue_class outer_tissue_class]
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


function mesh = bem_load_mesh(name)

bec = load(sprintf('%s.bec', name));
bee = load(sprintf('%s.bee', name));
bei = load(sprintf('%s.bei', name));

mesh.name = name;

if size(bei,2) ~= 4 || size(bei,1) ~= (bei(1,1) + 1)
    error('BEM:bem_load_mesh:info','Invalid Info file %s.bei',name);
end

mesh.num_nodes = bei(1,3);
mesh.num_elements = bei(1,2);
mesh.num_boundaries = bei(1,1);
mesh.num_node_elem = bei(1,4);

if size(bec,1) ~= mesh.num_nodes || size(bec,2) ~= 4
    error('BEM:bem_load_mesh:coord','Inconsistent Coordinate file %s.bec',name);
end
if size(bee,1) ~= mesh.num_elements || size(bee,2) ~= mesh.num_node_elem+1
    error('BEM:bem_load_mesh:elem','Inconsistent Element file %s.bee',name);
end

mesh.bnd = bei(2:mesh.num_boundaries+1,2:4);
mesh.num_class = max(max(mesh.bnd(:,2:3)));

% Check validity of boundary information
% All elements must be integer and > 0
%if ~isinteger(mesh.bnd) || ~isempty(find(mesh.bnd < 0))
if ~isempty(find(mesh.bnd < 0))
    error('BEM:bem_load_mesh:bnd','%s', 'Invalid tissue class');
end
% Sum of boundary elements must match total number of elements
if sum(mesh.bnd(:,1)) ~= mesh.num_elements
    error('BEM:bem_load_mesh:bnd','%s', 'Boundary size mismatch');
end

mesh.coord = bec(:,2:4);
mesh.elem = bee(:,2:mesh.num_node_elem+1);

a = dir('transform');
if size(a,1) > 0 
    mesh.transform = load('transform');
end

