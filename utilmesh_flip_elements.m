% utilmesh_flip_elements() - Flips elements with it's neighbor element
% according to aspect ratio in a given mesh
%
% Usage:
%   >> [Coord,Elem,number_of_flipped_elems] =
%   utilmesh_flip_elements(Coord,Elem);
%
% Inputs:
%   Coord - coordinate matrix of input mesh
%   Elem  - connectivity matrix of input mesh
%
% Outputs:
%   Coord - coordinate matrix of output mesh
%   Elem  - connectivity matrix of output mesh
%   number_of_flipped_elems - number of flipped elements
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


function [Coord,Elem,number_of_flipped_elems] = utilmesh_flip_elements(Coord,Elem);
% for linear meshes
nel = size(Elem,1);
nop = size(Elem,2);

thr_ang = 140/180*pi;
number_of_flipped_elems = 0;
for i = 1:nel
   A1 = Coord(Elem(i,2),2:4);
   A2 = Coord(Elem(i,3),2:4);
   A3 = Coord(Elem(i,4),2:4);
   a1 = angle2Lines(A2,A1,A3);
   a2 = angle2Lines(A1,A2,A3);
   a3 = angle2Lines(A2,A3,A1);
   ang = [a1 a2 a3];
   k = find(ang>thr_ang);
   if ~isempty(k) % if the element is too thin
      % find the edge in front of the large angle
      if k==1
         n1 = Elem(i,3); n2=Elem(i,4);
      elseif k==2
         n1 = Elem(i,2); n2=Elem(i,4);
      elseif k==3
         n1 = Elem(i,2); n2=Elem(i,3);
      end
      % find the element that have the edge n1-n2
      El = ElemsOfEdge(Elem,[n1 n2]);
      if (length(El) ~= 2)
         error('oops');
      end
      
      Crit1 = sum(ElemCritic(Elem(El,:),Coord));
      E = Elem(El,2:4); 
      Ef = flipElem(E); % flipped elements
      Elem(El,2:4) = Ef;
      
      EList = EdgeList(Coord,Elem,El);
      K = find(EList(:,5)~=0);
      Crit2 = sum(ElemCritic(Elem(El,:),Coord));
      if (~isempty(K))|(Crit2>Crit1)
         % disp('not improved...'); pause(0.1)
         Elem(El,2:4) = E;
      else
         %disp('flipped...');pause(0.1);
         number_of_flipped_elems = number_of_flipped_elems+1;
      end
   end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ang = angle2Lines(P1,P2,P3);
% finds the angle between the lines P1-P2 and P2-P3
P12=P2-P1;
P23=P2-P3;
normP12=norm(P12);
if normP12>eps
   P12=P12/normP12;  % normalize the unit vectors
end

normP23=norm(P23);
if normP23>eps
   P23=P23/normP23;  % normalize the unit vectors
end

ang=acos(P12(1)*P23(1)+P12(2)*P23(2)+P12(3)*P23(3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Elab1 = ElemsOfEdge(Elem,e1);
% finds the elements that have the edge e1 from Elem
% e1 : 1x2
% Elem : Nx4
nA = e1(1); nB = e1(2);
Elab1 = find(((Elem(:,2)==nA)&(Elem(:,3)==nB))|((Elem(:,2)==nB)&(Elem(:,3)==nA))|...
((Elem(:,3)==nA)&(Elem(:,4)==nB))|((Elem(:,3)==nB)&(Elem(:,4)==nA))|...
((Elem(:,2)==nA)&(Elem(:,4)==nB))|((Elem(:,2)==nB)&(Elem(:,4)==nA)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Crit = ElemCritic(Elem,Coord);
ne = length(Elem(:,1));
for k = 1:ne
   
   d(1) = EdgeLength(Elem(k,2), Elem(k,3), Coord);
   d(2) = EdgeLength(Elem(k,3), Elem(k,4), Coord);
   d(3) = EdgeLength(Elem(k,4), Elem(k,2), Coord);
   a = sort(d);
   Crit(k) = a(1) / a(3);
   
   A1 = Coord(Elem(k,2),2:4);
   A2 = Coord(Elem(k,3),2:4);
   A3 = Coord(Elem(k,4),2:4);
   a(1) = angle2Lines(A2,A1,A3);
   a(2) = angle2Lines(A1,A2,A3);
   a(3) = angle2Lines(A2,A3,A1);
   
	Crit(k) = Crit(k)+sum(abs(a-pi/3));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Len = EdgeLength(n1,n2,Coord);

p1=Coord(n1,2:4);
p2=Coord(n2,2:4)-p1;
Len=sqrt(p2*p2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Ef=flipElem(E);
% flip neighbor elements
% E matrix size: 2*3

K=unique(E(:)); n=0;
for i=1:length(K)
   a=find(E==K(i));
   if length(a)==2
      n=n+1;
      Ko(n)=K(i);
   end
end
% Ko is the common edge

% current edges
EL=[E(1,1) E(1,2);E(1,2) E(1,3);E(1,3) E(1,1);E(2,1) E(2,2);E(2,2) E(2,3);E(2,3) E(2,1)];

% Ke common edge
Ke=setdiff(K,Ko);
Ef=[Ke Ko(1);Ke Ko(2)];
a=findEdge(EL,[Ef(1,2) Ef(1,3)]);
if isempty(a)
   Ef(1,:)=fliplr(Ef(1,:));
end
a=findEdge(EL,[Ef(2,2) Ef(2,3)]);
if isempty(a)
   Ef(2,:)=fliplr(Ef(2,:));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function a=findEdge(KE,e1);
% KE : N x 2
% e1 : 1 x 2
X=find(KE(:,1)==e1(1));
Y=find(KE(:,2)==e1(2));
a=intersect(X,Y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EL=EdgeList(Coord,Elem,E);
% forms an edge list
% E is the list of the elements that you want to list their edges
% EL is the edge list
% the first two columns are the edges, the other columns are
% the elements that have that edge.
nel=size(Elem,1);
e1=[Elem(E(1),2) Elem(E(1),3)]; e1=sort(e1);
e2=[Elem(E(1),3) Elem(E(1),4)]; e2=sort(e2);
e3=[Elem(E(1),4) Elem(E(1),2)]; e3=sort(e3);
clear els; els=ElemsOfEdge(Elem,e1); els=sort(els); ne=length(els);
EL(1,1:2)=e1; EL(1,3:ne+2)=els'; EL(1,ne+3:ne+4)=0;
clear els; els=ElemsOfEdge(Elem,e2); els=sort(els); ne=length(els);
EL(2,1:2)=e2; EL(2,3:ne+2)=els';
clear els; els=ElemsOfEdge(Elem,e3); els=sort(els); ne=length(els);
EL(3,1:2)=e3; EL(3,3:ne+2)=els';
n=3;
for k=2:length(E)
   i=E(k);
   e1=[Elem(i,2) Elem(i,3)]; e1=sort(e1);
   a=findEdge(EL,e1);
   if isempty(a)
      n=n+1;
      clear els; els=ElemsOfEdge(Elem,e1); els=sort(els); ne=length(els);
		EL(n,1:2)=e1; EL(n,3:ne+2)=els';
   end
   e2=[Elem(i,3) Elem(i,4)]; e2=sort(e2);
   a=findEdge(EL,e2);
   if isempty(a)
      n=n+1;
      clear els; els=ElemsOfEdge(Elem,e2); els=sort(els); ne=length(els);
		EL(n,1:2)=e2; EL(n,3:ne+2)=els';
   end
   e3=[Elem(i,4) Elem(i,2)]; e3=sort(e3);
   a=findEdge(EL,e3);
	if isempty(a)
      n=n+1;
      clear els; els=ElemsOfEdge(Elem,e3); els=sort(els); ne=length(els);
		EL(n,1:2)=e3; EL(n,3:ne+2)=els';
   end
end




            
      
   



   