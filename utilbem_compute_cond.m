% utilbem_compute_cond() - Computes the average and difference
%     conductivity around a node used in BEM computations.
%
% Usage:
%   >> [condN, condE] = utilbem_compute_cond(Coord, Elem, layers, sigma)
%
% Inputs:
%    Coord  - node coordinate matrix
%    Elem   - nlement connectivity matrix
%    layers - mesh boundary information
%    sigma  - vector of element conductivities
%
% Outputs:
%    condN - average conductivity for each node (vector)
%    condE - difference conductivity (inner-outer)
%
% Note: The average conductivity is normally the average of the inner and outer
%    conductivities of the layer that a node belongs to. However, for
%    meshes with intersecting boundaries there may be three or more tissues
%    around a node. This function also handles this general case.
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

function [condN, condE] = utilbem_compute_cond(Coord, Elem, layers, sigma)
num_layers = size(layers,1);

% increse layer values to prevent 0 index and move sigma
layers(:,2:3) = layers(:,2:3) + 1;
sigma = [0 sigma(:)']';

bmax = max(max((layers(:,2:3))));
cmap = zeros(size(Coord,1), bmax);
condE = zeros(size(Coord,1), 1);

ei = 1;
for l = 1:num_layers
    for i = 1:layers(l,1)
        cmap(Elem(ei,:), layers(l,2)) = 1;
        cmap(Elem(ei,:), layers(l,3)) = 1;
        
        %% XXX does not work with intersecting meshes!
        condE(Elem(ei,:)) = sigma(layers(l,2)) - sigma(layers(l,3));
        ei = ei + 1;
    end
end

csum = sum(cmap,2);
condN = (cmap * sigma) ./ csum;
