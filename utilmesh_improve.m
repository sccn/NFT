% utilmesh_improve() - Improves mesh quality by checking and correcting
%           discontinuities
%
% Usage:
%   >> [Coord,Elem] = utilmesh_improve(Coord,Elem);
%
% Inputs:
%   Coord - coordinate matrix of input mesh
%   Elem  - connectivity matrix of input mesh
%
% Outputs:
%   Coord - coordinate matrix of output mesh
%   Elem  - connectivity matrix of output mesh
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

function [Coord,Elem] = utilmesh_improve(Coord,Elem);
nnp = size(Coord,1);
nel = size(Elem,1);
for i = 1:nel
   % find neighbour elements by edge
   N = findNeigE2(Elem,i);
   N = setdiff(N,i);
   NormN = ElemNormal(Elem(N,:),Coord);
   Normi = ElemNormal(Elem(i,:),Coord);
   Ti = NormN*Normi';
   h = find(Ti < -0.5); 
   % if the angle between ith element and neighbour element is more than
   if length(h) == 3
      % check its norm by flipping each element with ith element
      for j = 1:3
         ji = N(j); % j. neighbour element
         E = [Elem(i,2:4);Elem(ji,2:4)];
         Ef = flipElem(E);
         Elem(i,2:4) = Ef(1,:);
         Elem(ji,2:4) = Ef(2,:);
         NormN = ElemNormal(Elem(N,:),Coord);
         Normi = ElemNormal(Elem(i,:),Coord);
         T = NormN * Normi';
         sumT(j) = sum(T);
         Elem(i,2:4) = E(1,:);
         Elem(ji,2:4) = E(2,:);
      end
      [Y,I] = max(sumT);
      if Y > sum(Ti)
         ji = N(I);
	     E = [Elem(i,2:4);Elem(ji,2:4)];
   	     Ef = flipElem(E);
      	 Elem(i,2:4) = Ef(1,:);
         Elem(ji,2:4) = Ef(2,:);
         EL = EdgeList(Coord, Elem, [i ji]);
         K = find(EL(:,5) ~= 0);
         if ~isempty(K)
%            disp('not improved...'); pause(0.1)
            Elem(i,2:4) = E(1,:);
            Elem(ji,2:4) = E(2,:);
         end
      end
   end
end

Coord(:,1)=[1:length(Coord(:,1))]';
Elem(:,1)=[1:length(Elem(:,1))]';

  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function N=findNeigE2(Elem2,n);
% find the neighbour elements of the nth element of the mesh
% finds the only neighbours if any edge is common

n1=Elem2(n,2);
n2=Elem2(n,3);
n3=Elem2(n,4);

nA=n1; nB=n2;
Elab1=find(((Elem2(:,2)==nA)&(Elem2(:,3)==nB))|((Elem2(:,2)==nB)&(Elem2(:,3)==nA))|...
   ((Elem2(:,3)==nA)&(Elem2(:,4)==nB))|((Elem2(:,3)==nB)&(Elem2(:,4)==nA))|...
   ((Elem2(:,2)==nA)&(Elem2(:,4)==nB))|((Elem2(:,2)==nB)&(Elem2(:,4)==nA)));

nA=n2; nB=n3;
Elab2=find(((Elem2(:,2)==nA)&(Elem2(:,3)==nB))|((Elem2(:,2)==nB)&(Elem2(:,3)==nA))|...
   ((Elem2(:,3)==nA)&(Elem2(:,4)==nB))|((Elem2(:,3)==nB)&(Elem2(:,4)==nA))|...
   ((Elem2(:,2)==nA)&(Elem2(:,4)==nB))|((Elem2(:,2)==nB)&(Elem2(:,4)==nA)));

nA=n3; nB=n1;
Elab3=find(((Elem2(:,2)==nA)&(Elem2(:,3)==nB))|((Elem2(:,2)==nB)&(Elem2(:,3)==nA))|...
   ((Elem2(:,3)==nA)&(Elem2(:,4)==nB))|((Elem2(:,3)==nB)&(Elem2(:,4)==nA))|...
   ((Elem2(:,2)==nA)&(Elem2(:,4)==nB))|((Elem2(:,2)==nB)&(Elem2(:,4)==nA)));

N=union(Elab1,Elab2);
N=union(N,Elab3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Norm=ElemNormal(Elem,Coord);

ne=length(Elem(:,1));
nr=length(Coord(:,1));

for k=1:ne
   v1=Coord(Elem(k,3),2:4)-Coord(Elem(k,2),2:4);
   v2=Coord(Elem(k,4),2:4)-Coord(Elem(k,3),2:4);
   n=cross(v1, v2);
   Norm(k,1:3)=n/norm(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ef=flipElem(E);

K=unique(E(:)); n=0;
for i=1:length(K)
   a=find(E==K(i));
   if length(a)==2
      n=n+1;
      Ko(n)=K(i);
   end
end

% current edges
EL=[E(1,1) E(1,2);E(1,2) E(1,3);E(1,3) E(1,1);E(2,1) E(2,2);E(2,2) E(2,3);E(2,3) E(2,1)];

% common edge is Ke in new elements 
Ke = setdiff(K,Ko);
Ef = [Ke Ko(1);Ke Ko(2)];
a = findEdge(EL,[Ef(1,2) Ef(1,3)]);
if isempty(a)
   Ef(1,:) = fliplr(Ef(1,:));
end
a = findEdge(EL,[Ef(2,2) Ef(2,3)]);
if isempty(a)
   Ef(2,:) = fliplr(Ef(2,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a=findEdge(KE,e1);
% KE : N x 2
% e1 : 1 x 2
% looks for e1 in KE
X=find(KE(:,1)==e1(1));
Y=find(KE(:,2)==e1(2));
a=intersect(X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EL = EdgeList(Coord,Elem,E);
% forms an edge list
% E is the list of the elements that you want to list their edges
% EL is the edge list
% the first two columns are the edges, the other columns are
% the elements that have that edge.
nel = size(Elem,1);
e1 = [Elem(E(1),2) Elem(E(1),3)]; e1 = sort(e1);
e2 = [Elem(E(1),3) Elem(E(1),4)]; e2 = sort(e2);
e3 = [Elem(E(1),4) Elem(E(1),2)]; e3 = sort(e3);
clear els; els = ElemsOfEdge(Elem,e1); els = sort(els); ne = length(els);
EL(1,1:2) = e1; EL(1,3:ne+2) = els'; EL(1,ne+3:ne+4) = 0;
clear els; els = ElemsOfEdge(Elem,e2); els = sort(els); ne = length(els);
EL(2,1:2) = e2; EL(2,3:ne+2) = els';
clear els; els = ElemsOfEdge(Elem,e3); els = sort(els); ne = length(els);
EL(3,1:2) = e3; EL(3,3:ne+2) = els';
n = 3;
for k = 2:length(E)
   i = E(k);
   e1 = [Elem(i,2) Elem(i,3)]; e1 = sort(e1);
   a = findEdge(EL,e1);
   if isempty(a)
      n = n+1;
      clear els; els = ElemsOfEdge(Elem,e1); els = sort(els); ne = length(els);
		EL(n,1:2) = e1; EL(n,3:ne+2) = els';
   end
   e2 = [Elem(i,3) Elem(i,4)]; e2 = sort(e2);
   a = findEdge(EL,e2);
   if isempty(a)
      n = n+1;
      clear els; els = ElemsOfEdge(Elem,e2); els = sort(els); ne = length(els);
		EL(n,1:2) = e2; EL(n,3:ne+2) = els';
   end
   e3 = [Elem(i,4) Elem(i,2)]; e3 = sort(e3);
   a = findEdge(EL,e3);
	if isempty(a)
      n = n+1;
      clear els; els = ElemsOfEdge(Elem,e3); els = sort(els); ne = length(els);
		EL(n,1:2) = e3; EL(n,3:ne+2) = els';
   end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Elab1 = ElemsOfEdge(Elem,e1);
% finds the elements that have the edge e1 from Elem
% e1 : 1x2
% Elem : Nx4
nA = e1(1); nB = e1(2);
Elab1 = find(((Elem(:,2)==nA)&(Elem(:,3)==nB))|((Elem(:,2)==nB)&(Elem(:,3)==nA))|...
((Elem(:,3)==nA)&(Elem(:,4)==nB))|((Elem(:,3)==nB)&(Elem(:,4)==nA))|...
((Elem(:,2)==nA)&(Elem(:,4)==nB))|((Elem(:,2)==nB)&(Elem(:,4)==nA)));

