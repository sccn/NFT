% mesh_find_regions() - Find regions in a BEM mesh for Tetgen
%
% Usage:
%   >> R = mesh_find_regions(mesh)
%
% Inputs:
%   mesh - mesh
%   
% Outputs:
%   R - coordinates of points between the layers
%
%
% Author: Zeynep Akalin Acar, SCCN, 2010

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



function R = mesh_find_regions(mesh)
vol = [];
ncoordp = 0;
nelemp = 0;
for ii = 1:mesh.num_boundaries
    nelem = mesh.bnd(ii,1) + nelemp;
    vol.bnd(ii).tri = mesh.elem(nelemp+1:nelem,:) - ncoordp;
    ncoord = max(max(vol.bnd(ii).tri)) + ncoordp;
    vol.bnd(ii).pnt = mesh.coord(ncoordp+1:ncoord,:);
    ncoordp = ncoord;
    nelemp = nelem;
end

Ci(:,2:4) = vol.bnd(1).pnt;
Ci(:,1) = [1:length(Ci)]';

for ii = 1:mesh.num_boundaries-1
    clear C* E*
    % outer mesh for iith surface
    Co(:,2:4) = vol.bnd(ii).pnt;
    Co(:,1) = [1:length(Co)]';
    Eo(:,2:4) = vol.bnd(ii).tri;
    Eo(:,1) = [1:length(Eo)]';
    
    % inner mesh for iith surface
    Ci(:,2:4) = vol.bnd(ii+1).pnt;
    Ci(:,1) = [1:length(Ci)]';
    Ei(:,2:4) = vol.bnd(ii+1).tri;
    Ei(:,1) = [1:length(Ei)]';
    k1 = find(Co(:,4)==max(Co(:,4)));
    % max z in mesh Co
    p1 = Co(k1(1),2:4);
    K = Ci(:,2:4)-ones(length(Ci),1)*p1;
    M = sqrt(sum(K.*K,2));
    [h,j] = min(M);
    % closest point in mesh Ci
    p2 = Ci(j,2:4); 
    R(ii,:)=mean([p1; p2]);
end
if mesh.num_boundaries == 1
    ii = 0;
end
R(ii+1,:) = mean(Ci(:,2:4));
