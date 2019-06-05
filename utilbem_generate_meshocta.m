% utilbem_generate_meshocta() - This function generates coordinate and
%            connectivity files starting from a octahedron. 
%
% Usage:
%   >> [Coord,Elem] = utilbem_generate_meshocta(degr);
%
% Inputs:
%    degr  - degree of the mesh. If degr=1 -> number of nodes = 18
%                                If degr=2 -> number of nodes = 66
%                                If degr=3 -> number of nodes = 258
%                                If degr=4 -> number of nodes = 1026
%                                If degr=5 -> number of nodes = 4098
%
% Outputs:
%    Coord - Coordinate matrix of the nodes.
%    Elem - Element connectivity matrix.
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


function [Coord,Elem] = utilbem_generate_meshocta(degr)

b=sqrt(2)/2;
Coord=[1 0 0 1;
   2 -b b 0;
   3 -b -b 0;
   4 b -b 0;
   5 b b 0;
   6 0 0 -1];

Elem=[1 1 2 3;
   2 1 5 2;
   3 1 4 5;
   4 1 3 4;
   5 6 3 2;
   6 6 2 5;
   7 6 5 4;
   8 6 4 3];

nel=8;
nnp=6;
nlin=12;
newnode=[1     2     7
	      1     3     8
     		1     4     9
     		1     5    10
     		2     3    11
     		2     5    12
     		2     6    13
     		3     4    14
     		3     6    15
     		4     5    16
     		4     6    17
     		5     6    18];
  
for iter=1:degr
   % change coordinate matrix
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
   
   % change element matrix
   for i=1:nel
      [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,2))));
      nn(1)=newnode(I,3);
      [I,J]=find((newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,3)))|(newnode(:,1)==Elem(i,3)&(newnode(:,2)==Elem(i,4))));
      nn(2)=newnode(I,3);
      [I,J]=find((newnode(:,1)==Elem(i,2)&(newnode(:,2)==Elem(i,4)))|(newnode(:,1)==Elem(i,4)&(newnode(:,2)==Elem(i,2))));
      nn(3)=newnode(I,3);
       
      Elem(nel+1+(i-1)*3,1)=nel+1+(i-1)*3;
      Elem(nel+1+(i-1)*3,2)=nn(1);
      Elem(nel+1+(i-1)*3,3)=Elem(i,3);
      Elem(nel+1+(i-1)*3,4)=nn(2);
      Elem(nel+2+(i-1)*3,1)=nel+2+(i-1)*3;
      Elem(nel+2+(i-1)*3,2)=nn(3);
      Elem(nel+2+(i-1)*3,3)=nn(2);
      Elem(nel+2+(i-1)*3,4)=Elem(i,4);
      Elem(nel+3+(i-1)*3,1)=nel+3+(i-1)*3;
      Elem(nel+3+(i-1)*3,2)=nn(1);
      Elem(nel+3+(i-1)*3,3)=nn(2);
      Elem(nel+3+(i-1)*3,4)=nn(3);
      Elem(i,3)=nn(1);
      Elem(i,4)=nn(3);
      
      A(1+(i-1)*3,1)=nn(1);
      A(1+(i-1)*3,2)=nn(2);
      A(2+(i-1)*3,1)=nn(2);
      A(2+(i-1)*3,2)=nn(3);
      A(3+(i-1)*3,1)=nn(3);
      A(3+(i-1)*3,2)=nn(1);
   end
   nnp=nnp+nlin;
   nlin=2*nlin+3*nel;
   nel=nel*4;
   
   % find newnode matrix
   [sznn1 sznn2]=size(newnode);
   [szA1 szA2]=size(A);
   A(szA1+1:szA1+sznn1,1)=newnode(:,1);
   A(szA1+1:szA1+sznn1,2)=newnode(:,3);
   A(szA1+1+sznn1:szA1+2*sznn1,1)=newnode(:,2);
   A(szA1+1+sznn1:szA1+2*sznn1,2)=newnode(:,3);
   
   clear newnode
   newnode=zeros(nlin,3);
   newnode(1:nlin,1:2)=A;
   newnode(1:nlin,3)=(nnp+1:nnp+nlin)';
end


