% session2vol() - converts NFT session to Fieldtrip vol and sens structure
%
% Usage:
%   >> [vol, sens] = session2vol(session)
%
%   
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

function [vol, sens] = session2vol(session)

% session_name : session name (with the extension)

%session = bem_load_session(session_name);

mesh = session.model.mesh;
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

vol.type = 'metubem';
vol.cond = session.model.cond';
vol.session = session;

if isfield(session, 'sens')
    sens = session.sens;
else
    % compute from Smatrix
    S = session.Smatrix;
    ns = max(S(:,1));
    pnt = zeros(ns, 3);
    label = cell(1, ns);
    for i = 1:ns
        j = find(S(:,1) == i);
        pnt(i,:) = S(j,3)' *  mesh.coord(S(j,2),:);
        label{i} = sprintf('Sens-%d', i);
    end
    sens.label = label;
    sens.pnt = pnt;
end

sens.type = 'eeg';
