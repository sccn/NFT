% utilbem_compute_indices() - When the Isolated Problem Approach (IPA) is used,
%     the BEM equations are modified to reduce the numerical errors due to
%     the low conductivity skull layer.  For this purpose an "inner mesh"
%     is defined consisting of skull and the inner layers. 
%     Computes the node indices corresponding to these layers to be used
%     in BEM computations.
%
% Usage:
%   >> [indMsh, indMod, indMshMod] = utilbem_compute_indices(Coord, Elem, layers, mod);
%
% Inputs:
%    Coord  - node coordinate matrix
%    Elem   - element connectivity matrix
%    layers - mesh boundary information
%    mod    - index of the boundary modified for IPA.
%
% Outputs:
%    indMsh - coordinate indices of inner mesh nodes (vector).
%    indMod - coordinate indices of modified boundary nodes (vector).
%    indMshMod - coordinate indices of modified boundary relative to inner
%                mesh.
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

function [indMsh, indMod, indMshMod] = utilbem_compute_indices(Coord, Elem, layers, mod)

num_layers = size(layers,1);
cmap = zeros(size(Coord,1),1);

ei = 1;
for l = 1:num_layers
    for i = 1:layers(l,1)
        if l == mod
            cmap(Elem(ei,:)) = 1;
        elseif l > mod
            cmap(Elem(ei,:)) = 2;
        end
        ei = ei + 1;
    end
end

indMsh = find(cmap);
indMod = find(cmap == 1);

% We also need the index of boundary nodes
% relative to the inner mesh
[C,indMshMod,iB] = intersect(indMsh, indMod);

if length(C) ~= length(indMod)
    error('UtilBEM:utilbem_compute_indices:mesh','%s', ...
        'Mesh does not cover boundary nodes!');
end
