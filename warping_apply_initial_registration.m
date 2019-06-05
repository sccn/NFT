% warping_apply_initial_registration() - Transforms a set of points
%                 using scaling, translation and rotation parameters
%                 computed during initial registration
%
%
% Usage:
%   >> Pnew = warping_apply_initial_registration(sc, T, R, P);
%
% Inputs:
%   sc - scaling (scalar
%   T  - translation [x y z]
%   R  - rotation [x y z] in degrees
%   P  - point set
%
% Outputs:
%   Pnew - point set after scaling, translation and rotation of P
%
%
% Author: Zeynep Akalin Acar, SCCN, 2013

% Copyright (C) 2013 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
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

function Pnew = warping_apply_initial_registration(sc, T, R, P)

alpx = R(1) * pi / 180;
alpy = R(2) * pi / 180;
alpz = R(3) * pi / 180;

x = P(:,1) * sc + T(1);
y = P(:,2) * sc + T(2);
z = P(:,3) * sc + T(3);

% rotation around x-axis
x1 = x;
y1 = y*cos(alpx) - z * sin(alpx);
z1 = y*sin(alpx) + z * cos(alpx);

% rotation around y-axis
x2 = z1 * sin(alpy) + x1 * cos(alpy);
y2 = y1;
z2 = z1 * cos(alpy) - x1 * sin(alpy);

% rotation around z-axis
x3 = x2 * cos(alpz) - y2 * sin(alpz);
y3 = x2 * sin(alpz) + y2 * cos(alpz);
z3 = z2;

Pnew = [x3 y3 z3];

