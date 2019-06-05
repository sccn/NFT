% warping_distmeshpoint() - Computes distance between a mesh and a point.
%
% Usage:
%   >> [dm, Pm] = warping_distmeshpoint(P, Coord, Elem)
%
% Inputs:
%   P - Point
%   Coord, Elem - mesh
%
% Outputs:
%   dm - distance
%   Pm - closest point on the mesh
%
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

function [dm, Pm] = warping_distmeshpoint(P, Coord, Elem)
% Pm is the point on the mesh
% dm is the distance
% works for LINEAR MESH

% convert to linear mesh if it is quadratic
if size(Elem,2) == 7
    E(:,1:2) = Elem(:,1:2);
    E(:,3) = Elem(:,4);
    E(:,4) = Elem(:,6);
    np = max(max(E(:,2:4)));
    Coord = Coord(1:np,:);
    Elem = E; clear E;
end

nnp = size(Coord,1);

Xo = Coord(:,2:4) - ones(nnp,1) * P;
[Rad Ir] = sort(sum(Xo.*Xo,2));
N = Ir(1:4);

% find the neighbour elements of the closest nodes
E = mesh_elementsofthenodes(Coord,Elem,N);

% find the intersection of the line PP1 with the elements of E
nE=length(E);
Pint = zeros(nE,3); dis = zeros(nE,1);
for i = 1 : nE
   Pa = Coord(Elem(E(i),2),2:4);
   Pb = Coord(Elem(E(i),3),2:4);
   Pc = Coord(Elem(E(i),4),2:4);
   [D, Pp] = warping_disttrianglepoint(Pa,Pb,Pc,P);
   dis(i) = D;
   Pint(i,:) = Pp;
end
[j, k] = min(abs(dis));
Pm = Pint(k,:);
dm = abs(dis(k));

