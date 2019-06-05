% utilbem_PlotPot() - Plot potential map on the mesh
%
% Usage:
%   >> utilbem_PlotPot(Coord, Elem, Pot)
%
% Inputs:
%   Coord - Coordinate file
%   Elem - Mesh connectivity file
%   Pot - potential distribution
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

function utilbem_PlotPot(Coord, Elem, Pot)
% Coord : coordinate file
% Elem :  element connectivity file
% Pot : vector of potentials at Coord points

[nel,dum]=size(Elem); 

%figure
hold on
view(180,0)
nnode=size(Elem,2)-1;
for i=1:nel
   for n=1:nnode
      a=Elem(i,n+1);
		X(n)=Coord(a,2);
		Y(n)=Coord(a,3);
      Z(n)=Coord(a,4);
      Col(n)=Pot(a);
   end
   fill3(X',Y',Z',Col')
end
rotate3d

%shading('interp')
colormap default
xlabel('x')
ylabel('y')
zlabel('z')

