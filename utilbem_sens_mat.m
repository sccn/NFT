% utilbem_sens_mat() - Calculate sensor matrix (Smatrix).
%
% Usage:
%   >> S = utilbem_sens_mat(Coord, Elem, meshp);
% 
% Inputs:
%   Coord - Coordinate file
%   Elem - Connectivity file
%   meshp - digitizer points on the mesh (# of electrodes x 3 matrix)
%
% Outputs:
%   S - Smatrix.
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

function S = utilbem_sens_mat(Coord, Elem, meshp);

if size(Elem, 2) == 4 
    S = ubsm_sensMatL(Coord, Elem, meshp);
elseif size(Elem,2) == 7
    S = ubsm_sensMat(Coord,Elem, meshp);
else
    error('BEM:utilbem_sens_mat:mesh','%s','Invalid number of nodes per element');
end

% -------------------------------------------------------------
function d = ubsm_funG2L(x,gx,gy,gz,x0,y0,z0);
% callback function for GlobalToLocal
[xi,yi,zi] = ubsm_LocalToGlobal(x(1),x(2),gx,gy,gz);
d(1) = xi-x0;
d(2) = yi-y0;
d(3) = zi-z0;

% -------------------------------------------------------------
function [zet,eta] = ubsm_GlobalToLocal(xi,yi,zi,x,y,z);
% x,y,z : 6 points of the element
% xi,yi,zi : point
x = x(:); y = y(:); z = z(:);

options = optimset('Display','off');
X = lsqnonlin(@ubsm_funG2L,[0.3 0.3], [-inf -inf], [inf inf], options, x,y,z,xi,yi,zi);
%d=funG2L(x,gx,gy,gz,x0,y0,z0);
zet = X(1);
eta = X(2);

% -------------------------------------------------------------
function [xi,yi,zi] = ubsm_LocalToGlobal(zet,eta,x,y,z);
% x,y,z are the global coordinates
x = x(:); y = y(:); z = z(:);

%   ^ eta
%   |
% 0 o
%   |\
% 1 o o 5
%   |  \
%   o-o-o-> zet
%  2  3  4
%

if (length(x) == 3)
    Shape(1) = eta;
    Shape(2) = 1-zet-eta;
    Shape(3) = zet;
else
    Shape(1) = 2*eta*eta-eta;
    Shape(2) = -4*(eta*eta+zet*eta-eta);
    Shape(3) = 2*(zet*zet+eta*eta)+4*zet*eta-3*(zet+eta)+1;
    Shape(4) = -4*(zet*zet+zet*eta-zet);
    Shape(5) = 2*zet*zet-zet;
    Shape(6) = 4*zet*eta;
end

xi = Shape*x;
yi = Shape*y;
zi = Shape*z;

% -------------------------------------------------------------
function S = ubsm_sensMatL(Coord, Elem, meshp)

% calculation of the 'S' matrix
% Coord, Elem : linear mesh
% meshp : digitizer points on the mesh

H = zeros(length(meshp),3);

hh=waitbar(0,'calculating sensor matrix...');
for i = 1 : length(meshp)
    waitbar(i/length(meshp));
    [dm,Pm,el] = ubsm_DistMeshPoint(meshp(i,:),Coord,Elem);
    x = Coord(Elem(el,2:4),2);
    y = Coord(Elem(el,2:4),3);
    z = Coord(Elem(el,2:4),4);
    [zet,eta] = ubsm_GlobalToLocal(meshp(i,1),meshp(i,2),meshp(i,3),x,y,z);
    H(i,1:3) = [el zet eta];
end
close(hh);

%nds=unique(Elem2(H(:,1),2:7));
nn = size(Coord, 1);
ns = size(H,1);

Smat = zeros(ns, nn);

for i = 1:ns
    el = H(i,1);
    zet = H(i,2);
    eta = H(i,3);
    
    Shape(1) = eta;
    Shape(2) = 1-zet-eta;
    Shape(3) = zet;
    
    for j = 1:3
        Smat(i, Elem(el, j+1)) = Shape(j);
    end
end

[i,j,s] = find(Smat');
S = [j i s];

% -------------------------------------------------------------
function S = ubsm_sensMat(Coord2,Elem2, meshp);
% calculation of the 'S' matrix (sensor matrix)
% Coord2, Elem2 : quadratic mesh
% Coord, Elem : linear mesh
% meshp : digitizer points on the mesh (# of electrodes x 3 matrix)

[Coord, Elem] = ubsm_MakeLinearMesh(Coord2,Elem2);
hh=waitbar(0,'calculating sensor matrix...');
for i = 1:length(meshp)
    waitbar(i/length(meshp));
    [dm,Pm,el] = ubsm_DistMeshPoint(meshp(i,:),Coord,Elem);
    x = Coord2(Elem2(el,2:7),2);
    y = Coord2(Elem2(el,2:7),3);
    z = Coord2(Elem2(el,2:7),4);
    [zet,eta] = ubsm_GlobalToLocal(meshp(i,1),meshp(i,2),meshp(i,3),x,y,z);
    H(i,1:3) = [el zet eta];
end
close(hh);


%nds=unique(Elem2(H(:,1),2:7));
nds=Coord2(:,1);
nn = length(nds);
ns = size(H,1);

%calculate a map from node index to vector index
ndmap = zeros(1,max(nds));
for i = 1:nn
    ndmap(nds(i)) = i;
end

Smat = zeros(ns, nn);

for i = 1:ns
    el = H(i,1);
    zet = H(i,2);
    eta = H(i,3);
    
    Shape(1) = 2*eta*eta-eta;
    Shape(2) = -4*(eta*eta+zet*eta-eta);
    Shape(3) = 2*(zet*zet+eta*eta)+4*zet*eta-3*(zet+eta)+1;
    Shape(4) = -4*(zet*zet+zet*eta-zet);
    Shape(5) = 2*zet*zet-zet;
    Shape(6) = 4*zet*eta;
    
    for j = 1:6
        Smat(i, ndmap(Elem2(el, j+1))) = Shape(j);
    end
end


[i,j,s] = find(Smat');
S = [j i s];

% -------------------------------------------------------------
function [Coordl,Eleml] = ubsm_MakeLinearMesh(Coordq,Elemq)

nel = size(Elemq,1);
nnpq = size(Coordq,1);
Eleml = zeros(nel,4);
Eleml(:,1:2) = Elemq(:,1:2);
Eleml(:,3) = Elemq(:,4);
Eleml(:,4) = Elemq(:,6);

E = Eleml(:,2:4);
Ei = sort(unique(E(:)));
Ej = zeros(nnpq,1);
Ej(Ei) = 1:length(Ei);

Coordl = Coordq(Ei,:);
nnp = size(Coordl,1);

for i=1:nel
    Eleml(i,2:4) = Ej(Eleml(i,2:4));
end

Coordl(:,1) = 1:nnp;

%-------------------------------------------------------------
function [dm,Pm,el,in] = ubsm_DistMeshPoint(P,Coord,Elem);
% looks for if P is inside the mesh Coord, Elem or not
% Pm is the point on the mesh
% dm is the distance
% el is the element of Pm
% in = inside (bool)
% works for LINEAR MESH

nnp = size(Coord,1);

r = 5;
N = 0;
Xo = Coord(:,2:4)-ones(nnp,1)*P;
Rad = sum(Xo.*Xo,2);
while length(N)<3
    r=r+5;
    N=find(Rad<r^2);
end

% find the neighbour elements of the closest nodes
E = ubsm_ElementsOfTheNodes(Coord,Elem,N);

% find the intersection of the line PP1 with the elements of E
Pint = []; dis = [];
for i = 1:length(E)
   Pa = Coord(Elem(E(i),2),2:4);
   Pb = Coord(Elem(E(i),3),2:4);
   Pc = Coord(Elem(E(i),4),2:4);
   [D,Pp] = ubsm_DistTrianglePoint2(Pa,Pb,Pc,P);
   dis(i) = D;
   Pint(i,:) = Pp;
end
[j,k] = min(abs(dis));
Pm = Pint(k,:);
dm = abs(dis(k));
el = E(k);

% check if P is inside or outside

% find the normal vector of the element
v1 = Coord(Elem(el,3),2:4) - Coord(Elem(el,2),2:4);
v2 = Coord(Elem(el,4),2:4) - Coord(Elem(el,3),2:4);
n = cross(v1, v2);
Norm = n/norm(n);

% if the vector PPm.Norm > 0 inside, < 0 outside
if dot(Pm-P, Norm) > 0 
    in = 1; % P is inside the mesh
else
    in = 0;
end

% -------------------------------------------------------------
function E = ubsm_ElementsOfTheNodes(Coord,Elem,A);
% A is the list of the nodes
% E is the list of the elements that have nodes in A
E = [];
nop = size(Elem,2);
if nop == 4
   for i = 1:length(A)
       n1 = find(Elem(:,2)==A(i));
	   n2 = find(Elem(:,3)==A(i));
       n3 = find(Elem(:,4)==A(i));
	   n4 = union(n1,n2);
       n5 = union(n3,n4);
	   E = union(E,n5);
      clear n1 n2 n3 n4 n5
   end
elseif nop == 7
   for i = 1:length(A)
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

% -------------------------------------------------------------
function [D,Pp]=ubsm_DistTrianglePoint2(Pa,Pb,Pc,Px);
% finds the minimum distance of a point Px to triangle Pa, Pb, Pc
% difference from DistTrianglePoint 
% doesn't look if the projection of the point is in the triangle or on the edge 
% of the triangle otherwise MinD is 1000

% find the minimum distance of the point with the 
% plane which is formed by the triangle
% find the normal of the plane
eps=1e-4;
v1=Pa-Pb;
v2=Pc-Pb;
n=[v1(2)*v2(3)-v2(2)*v1(3) v2(1)*v1(3)-v1(1)*v2(3) v1(1)*v2(2)-v1(2)*v2(1)];
n=n/norm(n);
d=-n(1)*Pa(1)-n(2)*Pa(2)-n(3)*Pa(3);
% the plane equation is n(1)*x+n(2)*y+n(3)*z+d=0
% distance of the point to the plane is D
D=(n(1)*Px(1)+n(2)*Px(2)+n(3)*Px(3)+d)/sqrt(n(1)^2+n(2)^2+n(3)^2);
% point on the plane
Pp=Px-D*n;
% check if the point is on the triangle
% Determine whether or not the intersection point is bounded by pa,pb,pc 
Pa1=Pa-Pp;
normPa1=norm(Pa1);
if normPa1>eps
   % normalize the unit vectors
   Pa1=Pa1/normPa1;  
end
Pa2 = Pb - Pp;
normPa2=norm(Pa2);
if normPa2>eps
   Pa2=Pa2/normPa2; 
end
Pa3 = Pc - Pp;
normPa3=norm(Pa3);
if normPa3>eps
   Pa3=Pa3/normPa3;
end
%the angles are 
a1 = acos(Pa1(1)*Pa2(1) + Pa1(2)*Pa2(2) + Pa1(3)*Pa2(3));
a2 = acos(Pa2(1)*Pa3(1) + Pa2(2)*Pa3(2) + Pa2(3)*Pa3(3));
a3 = acos(Pa3(1)*Pa1(1) + Pa3(2)*Pa1(2) + Pa3(3)*Pa1(3));

total = a1+a2+a3;
% if total is 2*pi then the point is in the triangle or on the edges
if (abs(total - 2*pi) < eps)
   MinD=abs(D);
else
    % find the distance of Px to the edges
   [Di,Ppl]=ubsm_DistPointLineSegment(Px,Pa,Pb);
   Dix(1)=Di; Ppi(1,:)=Ppl;
   [Di,Ppl]=ubsm_DistPointLineSegment(Px,Pa,Pc);
   Dix(2)=Di; Ppi(2,:)=Ppl;
   [Di,Ppl]=ubsm_DistPointLineSegment(Px,Pc,Pb);
   Dix(3)=Di; Ppi(3,:)=Ppl;
   [u,v]=min(Dix);
   MinD=u;
   Pp=Ppi(v,:);
end
D=MinD;

% -------------------------------------------------------------
function [D,Pp]=ubsm_DistPointLineSegment(P,P1,P2);
% Finds distance between the point P and the line segment P1P2

v1=P2-P1;
n1 = v1/norm(v1);

v2 = P-P1;

% projection (component) of v2 along v1 (n1)
vp = dot(n1, v2);

% make sure it is in the range
if (vp < 0)
    vp = 0;
elseif (vp > norm(v1))
    vp = norm(v1);
end


% closest point on edge
Pp = P1 + vp * n1;

D = norm(P-Pp);

% -------------------------------------------------------------
