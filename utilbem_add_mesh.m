% utilbem_add_mesh() - This function adds two meshes and generate a single
%                      mesh.
%
% Usage:
%   >> [Coord3,Elem3] = utilbem_add_mesh(Coord1, Elem1, Coord2,Elem2);
%
% Inputs:
%    Coord1 - coordinates of the first mesh.
%    Elem1  - element connectivity matrix of the first mesh.
%    Coord2 - coordinates of the second mesh.
%    Elem2  - element connectivity matrix of the second mesh.
%
% Outputs:
%    Coord3 - coordinates of the resulting mesh.
%    Elem3  - element connectivity matrix of the resulting mesh.
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


function [Coord3,Elem3] = utilbem_add_mesh(Coord1, Elem1, Coord2,Elem2);

[nel1, nne] = size(Elem1);
[nnp1, dum] = size(Coord1);
[nel2, nne] = size(Elem2);
[nnp2, dum] = size(Coord2);

Elem3 = zeros(nel1+nel2,nne);
Coord3 = zeros(nnp1+nnp2,dum);

Elem2(:,1) = Elem2(:,1) + nel1;
Elem2(:,2:nne) = Elem2(:,2:nne) + nnp1;
Coord2(:,1) = Coord2(:,1) + nnp1;

Elem3(1:nel1,:) = Elem1;
Coord3(1:nnp1,:) = Coord1;
Elem3(nel1+1:nel1+nel2,:) = Elem2;
Coord3(nnp1+1:nnp1+nnp2,:) = Coord2;
nnp = size(Coord3,1);
Coord3(:,1) = (1:nnp)';