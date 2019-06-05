% mesh_elementsofthenodes() - finds the elements that are connected to
% given nodes
%
% Usage:
%   >> E = mesh_elementsofthenodes(Coord,Elem,A);
%
% Inputs:
%   Coord, Elem - mesh
%   A - list of the noes
%
% Outputs:
%   E - list of elements
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

function E = ElementsOfTheNodes(Coord,Elem,A);
% A is the list of the nodes
% E is the list of the elements that have nodes in A
E = [];
nop = size(Elem,2);
if nop == 4
   for i = 1 : length(A)
   	   n1 = find(Elem(:,2)==A(i));
	   n2 = find(Elem(:,3)==A(i));
   	   n3 = find(Elem(:,4)==A(i));
	   n4 = union(n1,n2);
   	   n5 = union(n3,n4);
	   E = union(E,n5);
       clear n1 n2 n3 n4 n5
   end
elseif nop == 7
   for i = 1 : length(A)
   	   n1 = find(Elem(:,2)==A(i));
	   n2 = find(Elem(:,3)==A(i));
       n3 = find(Elem(:,4)==A(i));
       n4 = find(Elem(:,5)==A(i));
   	   n5 = find(Elem(:,6)==A(i));
   	   n6 = find(Elem(:,7)==A(i));
	   n7 = union(n1,n2);
       n8 = union(n7,n3);
       n9 = union(n8,n4);
       n10 = union(n9,n5);
       n11 = union(n10,n6);
	   E = union(E,n11);
       clear n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11
   end
end

   


   