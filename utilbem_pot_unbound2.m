% utilbem_pot_unbound() - Computes the unbounded potential at the
%     given coordinates arising from activation of a single dipole. The 
%     result is weighted by the average conductivity around the node.
%
% Usage:
%   >> pot = utilbem_pot_unbound2(Coord, dcond, dip)
%
% Inputs:
%    Coord - node coordinate matrix
%    dcond  - derivative of inverse average conductivity for each node (vector)
%    dip   - dipole [x y z px py pz]
%
% Outputs:
%    pot -  node potentials for the homogeneous model
%
% Notes: When IPA is not used the weighted potential is the Right-Hand-Side
%         of the BEM matrix equation .
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

function pot = utilbem_pot_unbound2(Coord, dcond, dip)
% cond: average conductivity for each node (vector)
% dip: dipole [x y z px py pz]
% potential unbounded

% substract dipole coordinate from mesh Coord
R = Coord - ones(size(Coord,1),1)*dip(1:3);
MagR = sum(R.*R,2).^1.5;
pot = 1/4/pi*(R*dip(4:6)').*dcond./MagR;
