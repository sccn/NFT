% utilbem_form_quad_mesh() - This function generates a quadratic mesh from
%                   a linear mesh for spherical surfaces.
%
% Usage:
%   >> [Coord,Elem] = utilbem_form_quad_mesh(Coord,Elem);
%
% Inputs:
%    Coord - Coordinates of the linear mesh.
%    Elem - Element connectivity matrix of the linear mesh.
%
% Outputs:
%    Coord - Coordinates of the quadratic mesh.
%    Elem - Element connectivity matrix of the quadratic mesh.
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

function [Coord,Elem] = utilbem_form_quad_mesh(Coord,Elem);


[nel,dum]=size(Elem);
[nnp,dum]=size(Coord);

newnode=zeros(1,3);
% form the newnode matrix
n=1;
for i=1:nel
   [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,2))));
   if isempty(I)
      newnode(n,1)=Elem(i,2);
      newnode(n,2)=Elem(i,3);
      n=n+1;
   end
   [I,J]=find((newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,4))));
   if isempty(I)
      newnode(n,1)=Elem(i,3);
      newnode(n,2)=Elem(i,4);
      n=n+1;
   end
   [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,4)))|(newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,2))));
   if isempty(I)
      newnode(n,1)=Elem(i,2);
      newnode(n,2)=Elem(i,4);
      n=n+1;
   end
end
nlin=n-1;
newnode(:,3)=(nnp+1:nnp+nlin)';

% add nodes between the two nodes
for i=1:nlin
   l=newnode(i,1);
   k=newnode(i,2);
	X=[Coord(l,2) Coord(k,2)];
	Y=[Coord(l,3) Coord(k,3)];
	Z=[Coord(l,4) Coord(k,4)];
  	R=sqrt((mean(X))^2+(mean(Y))^2+(mean(Z))^2);
   Coord(nnp+i,1)=nnp+i;
   Coord(nnp+i,2)=mean(X)/R;
   Coord(nnp+i,3)=mean(Y)/R;
   Coord(nnp+i,4)=mean(Z)/R;
end

Elemqd=zeros(nel,7);
Elemqd(:,1:2)=Elem(:,1:2);
Elemqd(:,4)=Elem(:,3);
Elemqd(:,6)=Elem(:,4);
% add new nodes to element matrix
for i=1:nel
   [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,2))));
   Elemqd(i,3)=newnode(I,3);
   [I,J]=find((newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,4))));
   Elemqd(i,5)=newnode(I,3);
   [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,4)))|(newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,2))));
   Elemqd(i,7)=newnode(I,3);

end
nnp=nnp+nlin;
nlin=2*nlin;
Elem=Elemqd;


