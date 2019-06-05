% Warping_DistTrianglePoint() - Computes the distance between a triangle
% and a point.
%
% Usage:
%   >> [D, Pp] = Warping_DistTrianglePoint(Pa, Pb, Pc, Px);
%
% Inputs:
%   Pa, Pb, Pc - Coordinates of three corners of the element
%   Px - Point
%
% Outputs:
%   D - distance of the point to the plane
%   Pp - point on the plane
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

function [D, Pp] = Warping_DistTrianglePoint(Pa, Pb, Pc, Px);
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
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pb);
   Dix(1)=Di; Ppi(1,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pc);
   Dix(2)=Di; Ppi(2,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pc,Pb);
   Dix(3)=Di; Ppi(3,:)=Ppl;
   [u,v]=min(Dix);
   MinD=u;
   Pp=Ppi(v,:);
end
D=MinD;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [D,Pp]=DistPointLineSegment(P,P1,P2)
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


