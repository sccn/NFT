% mesh_check_intersection() - Detects and corrects mesh intersections
%
% Usage:
%   >> [so2, k1,k2] = mesh_check_intersection(so, C1, E1);
%
% Inputs:
%   so - mesh nodes
%   C1, E1 : mesh
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

function [so2, k1,k2] = mesh_check_intersection(so, C1, E1);
% so : sourcespace C1, E1 : mesh
% corrects source space
so2 = so;
[dim, inm] = utilmesh_check_source_space(so(:,1:3), C1, E1);
k1 = find(inm == 0); k1x = find(inm==1);
%if length(k1x) < length(k1);    k1 = k1x;end
no_intnodes = length(k1);

k2 = find(dim < 1); % find the nodes closer than 1mm
k2 = setdiff(k2, k1);
no_closenodes = length(k2);

if length(no_intnodes)>0
    for i = 1:no_intnodes
        p1 = so(k1(i),1:3);
        [dm, Pm, el, in] = utilmesh_dist_mesh_point(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        pn = Pm+nor;
        [dm1, Pm1, el1, in1] = utilmesh_dist_mesh_point(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k1(i); % intersecting nodes
        end
        so2(k1(i),1:3) = pn;
    end
end

if length(no_closenodes)>0
    for i = 1:no_closenodes
        p1 = so(k2(i),1:3);
        [dm, Pm, el, in] = utilmesh_dist_mesh_point(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        pn = Pm-nor;
        [dm1, Pm1, el1, in1] = utilmesh_dist_mesh_point(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k2(i); % close nodes
        end
        so2(k2(i),1:3) = pn;
    end
end
