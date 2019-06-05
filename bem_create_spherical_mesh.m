% bem_create_spherical_mesh() - Create a spherical BEM mesh.
%
% Usage:
%   >> mesh = bem_create_spherical_mesh(name, radii, order);
% 
% Inputs:
%   name - mesh name excluding the extension.
%   radii - vector of radii of the spheres.
%   order - order of the mesh. Either 'l' for linear or 'q' for quadratic meshes.
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


function mesh = bem_create_spherical_mesh(name, radii, order)

radii = sort(radii,'descend');

mesh.name = name;
mesh.num_boundaries = length(radii);

if order == 'l'
    mesh.num_node_elem = 3;
elseif order == 'q'
    mesh.num_node_elem = 6;
else
    error('BEM:bem_create_spherical_mesh:order','Order is not defined %s', order);
end

if mesh.num_node_elem == 3
    [Coord,Elem] = utilbem_generate_meshocta(4);
else
    [Coord,Elem] = utilbem_generate_meshocta(3);
    [Coord,Elem] = utilbem_form_quad_mesh(Coord,Elem);
end

Coord1(:,2:4) = radii(1) * Coord(:,2:4);
Elem1 = Elem;
mesh.bnd(1,:) = [size(Elem1,1) 1 0];
if mesh.num_boundaries > 1
    for layer = 1:mesh.num_boundaries-1
        Coord2(:,2:4) = radii(layer+1) * Coord(:,2:4);
        Elem2 = Elem;
        mesh.bnd(layer+1,:) = [size(Elem2,1) layer+1 layer];
        [Coord1, Elem1] = utilbem_add_mesh(Coord1, Elem1, Coord2, Elem2);
    end
end

mesh.num_nodes = size(Coord1, 1);
mesh.num_elements = size(Elem1,1);
mesh.num_class = mesh.num_boundaries;
mesh.num_elements = size(Elem1, 1);

mesh.coord = Coord1(:,2:4);
mesh.elem = Elem1(:,2:mesh.num_node_elem + 1);

nb = mesh.num_boundaries;
Info = int32(zeros(nb + 1, 4));
Info(2:nb+1,1) = 1:nb;
Info(2:nb+1,2:4) = mesh.bnd;
Info(1,:) = [nb mesh.num_elements mesh.num_nodes mesh.num_node_elem];

save([name '.bec'], 'Coord1', '-ascii');

f = fopen([name '.bei'], 'w');
fprintf(f, '%d %d %d %d\r\n', Info');
fclose(f);

noe=size(Elem1,2);
fid=fopen([name '.bee'],'w');
if noe == 7
    fprintf(fid, '%d %d %d %d %d %d %d\r\n', Elem1');
elseif noe == 4
    fprintf(fid, '%d %d %d %d\r\n', Elem1');
end
fclose(fid);

