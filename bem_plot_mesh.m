% bem_plot_mesh() - Plot a BEM mesh.
%
% Usage:
%   >> bem_plot_mesh(mesh, subE);
%
% Inputs:
%    mesh - mesh structure
%    subE - [optional] the indices of the elements that will be plotted in
%           a different color.
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

function bem_plot_mesh(mesh, subE);

if exist('subE')~=1
   subE=[]; 
end

Coord = zeros(mesh.num_nodes,4);
Elem = zeros(mesh.num_elements,mesh.num_node_elem+1);
Coord(:,1) = [1:mesh.num_nodes]';
Elem(:,1) = [1:mesh.num_elements]';
Coord(:,2:4) = mesh.coord;
Elem(:,2:mesh.num_node_elem+1) = mesh.elem;

c=0;
[nel,dum]=size(Elem);
[nnp,dum]=size(Coord);
figure
view(0,0)
hold on
nnode=size(Elem,2)-1;
for i=1:nel
   for n=1:nnode-c
      a=Elem(i,n+1);
      a=find(Coord(:,1)==a);
      X(n)=Coord(a,2); 
      Y(n)=Coord(a,3); 
      Z(n)=Coord(a,4);
   end
   if ~isempty(subE)
      t2=find(subE==Elem(i,1));
   else
      t2=[];
   end

   %   if isempty(subE) | isempty(find(subE==Elem(i,1)))
   if length(t2)>0
      fill3(X,Y,Z,[0.5 0.9 0.9])
   else
      fill3(X,Y,Z,[0.8 0.9 1])
   end
end
rotate3d on
xlabel('x'); ylabel('y'); zlabel('z');

if ~isempty(subE)
    figure; hold
    for i=1:length(subE)
        j=subE(i);
        for n=1:nnode-c
            a=Elem(j,n+1);
            a=find(Coord(:,1)==a);
            X(n)=Coord(a,2); 
            Y(n)=Coord(a,3); 
            Z(n)=Coord(a,4);
        end
        fill3(X,Y,Z,[0.5 0.9 0.9])
    end
end



