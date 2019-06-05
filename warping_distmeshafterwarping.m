% warping_distmeshafterwarping() - Finds the distance between a mesh and translated and
% rotated point set.
%
%
% Usage:
%   >> [F2, dmi] = warping_distmeshafterwarping(X, F, Coord, Elem);
%
% Inputs:
%   F - point set
%   Coord, Elem - mesh
%   X - translation and rotations
%
% Outputs:
%   dmi - distance vector
%   F2 - translated rotated point set.
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

function [F2, dmi] = warping_distmeshafterwarping(X, F, Coord, Elem);
% F : digitizer points
% F2 : are the rotated and translated version of F according to X
% they are moved to the closest points on the mesh.

tx = X(1);
ty = X(2);
tz = X(3);
alpx = X(4) * pi / 180;
alpy = X(5) * pi / 180;
alpz = X(6) * pi / 180;

x2 = F(:,1);
y2 = F(:,2);
z2 = F(:,3);

% rotation around x-axis
x3 = x2;
y3 = y2 * cos(alpx) - z2 * sin(alpx);
z3 = y2 * sin(alpx) + z2 * cos(alpx);

% rotation around y-axis
x4 = z3 * sin(alpy) + x3 * cos(alpy);
y4 = y3;
z4 = z3 * cos(alpy) - x3 * sin(alpy);

% rotation around z-axis
x5 = x4 * cos(alpz) - y4 * sin(alpz);
y5 = x4 * sin(alpz) + y4 * cos(alpz);
z5 = z4;

% translation
x1 = x5 + tx;
y1 = y5 + ty;
z1 = z5 + tz;

F2 = [x1 y1 z1];

hh = waitbar(0,'calculating the distance between digitizer and mesh');
for i = 1 : size(F2,1);
    waitbar(i/length(F2));
    [dm, Pm] = warping_distmeshpoint(F2(i,:),Coord,Elem);
    dmi(i) = dm;    Pmi(i,:) = Pm;
end; 
close(hh);
F2 = Pmi;




