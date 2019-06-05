% utilbem_multilayer_rhs() - Computes the Right-Hand-Side of the
%     BEM matrix equation when IPA (Isolated Problem Approach) is
%     used. The RHS vector corrected by IPA reduces numerical errors
%     resulting from the low conductivity skull layer.
%
% Usage:
%   >> RHS = utilbem_multilayer_rhs(Coord, cond, indMod, indMsh, indMshMod, iinv, dmt, s2, s3, dip);
%
% Inputs:
%    Coord     - node coordinate matrix
%    cond      - average conductivity for each node (vector)
%                as computed by utilbem_compute_cond()
%    indMod    - coordinate indices of modified boundary nodes (vector).
%    indMsh    - coordinate indices of inner mesh nodes (vector).
%    indMshMod - coordinate indices of modified boundary relative to inner mesh.
%                The indices are computed by utilbem_compute_indices()
%    iinv      - inverse inner matrix.
%    dmt       - sub coefficient matrix (# of nodes) x (# of mod. boundary nodes).
%    s2        - outer conductivity of the modified layer.
%    s3        - inner conductivity of the modified layer
%    dip       - dipole parameters [x y z px py pz]
%
% Outputs:
%    RHS - right-hand-side of the BEM matrix equation for a given
%    dipole.
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

function RHS = utilbem_multilayer_rhs(Coord, cond, indMod, indMsh, indMshMod, iinv, dmt, s2, s3, dip)

iCoord = Coord(indMsh,:);
icond = cond(indMsh,:);

icond(indMshMod,:) = icond(indMshMod) - (s2/2 * ones(length(indMshMod),1));

beta = s2/s3;
if beta == 1
    RHS = utilbem_pot_unbound(Coord, cond, dip);
    return;
end
mod = beta / (beta-1);
iRHS = utilbem_pot_unbound(iCoord, icond, dip);
ipot = iinv * iRHS;
RHS = mod * dmt * ipot(indMshMod);
mod = beta / (beta + 1);
RHS(indMod,:) = RHS(indMod,:) - mod * ipot(indMshMod);
