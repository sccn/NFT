% utilmesh_correct() - Performs topological correction for a given mesh.
% Checks intersecting elements, flips elements according to aspect ratio.
%
% Usage:
%   >> [Coord,Elem] = utilmesh_correct(Coord, Elem);
%
% Inputs:
%   Coord - coordinate matrix of input mesh
%   Elem  - connectivity matrix of input mesh
%
% Outputs:
%   Coord - coordinate matrix of output mesh
%   Elem  - connectivity matrix of output mesh
%
% Author: Zeynep Akalin Acar, SCCN, 2008

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

function [Coord,Elem] = utilmesh_correct(Coord, Elem);
ne=1;
while ne > 0
    [Coord, Elem, ne] = utilmesh_flip_elements(Coord,Elem);
end

[Coord, Elem, Check] = utilmesh_modify(Coord, Elem, []);
[Coord, Elem, Check] = utilmesh_modify(Coord, Elem, Check);
[Coord, Elem] = utilmesh_improve(Coord,Elem);
ne = 1;
while ne > 0
    [Coord, Elem, ne] = utilmesh_flip_elements(Coord,Elem);
end

[Coord, Elem, Check] = utilmesh_modify(Coord, Elem, []);
[Coord, Elem, Check] = utilmesh_modify(Coord, Elem, Check);
[Coord, Elem] = utilmesh_improve(Coord,Elem);
ne = 1;
while ne > 0
    [Coord, Elem, ne] = utilmesh_flip_elements(Coord,Elem);
end

