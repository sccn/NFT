% mesh_local_refinement() - Refines the meshes in a given folder. Loads
% Scalp.smf, Skull.smf, Csf.smf, Brain.smf meshes in smf format and saves 
% them with the same names. 
%
% Usage:
%   >> mesh_local_refinement(of, nl, ratio_lmr);
%
% Inputs:
%   of - mesh folder
%   nl - number of layers (3 or 4)
%   ratio_lmr - ratio of local edge length to local distance between meshes 
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


function mesh_local_refinement(of, nl, ratio_lmr)

% call from the folder where the meshes are
% nl is the number of layers ,3 or 4

if not(nl==3|nl==4)
    error('number of layers should be entered as 3 or 4');
end

[C1,E1] = mesh_readsmf([of 'Scalp.smf'],0,0,0,1);       % subject's scalp mesh
[C2,E2] = mesh_readsmf([of 'Skull.smf'], 0, 0, 0, 1);   % subject's skull mesh
[C3,E3] = mesh_readsmf([of 'Csf.smf'], 0, 0, 0, 1);     % subject's CSF mesh
[C4,E4] = mesh_readsmf([of 'Brain.smf'], 0, 0, 0, 1);   % subject's brain mesh

% refine and correct
if nl == 4
    disp('refining brain and csf...')
    [C4, E4, C3, E3] = utilmesh_refine_correct(C4, E4, C3, E3,ratio_lmr,of);
    disp('refining csf and skull...')
    [C3, E3, C2, E2] = utilmesh_refine_correct(C3, E3, C2, E2,ratio_lmr,of);
        disp('refining skull and scalp...')
    [C2, E2, C1, E1] = utilmesh_refine_correct(C2, E2, C1, E1,ratio_lmr,of);
else
    [C4, E4, C2, E2] = utilmesh_refine_correct(C4, E4, C2, E2,ratio_lmr,of);
    [C2, E2, C1, E1] = utilmesh_refine_correct(C2, E2, C1, E1,ratio_lmr,of);
end

% write meshes
Mesh_WriteSMF(of, 'Scalp.smf', C1, E1);
Mesh_WriteSMF(of, 'Skull.smf', C2, E2);
Mesh_WriteSMF(of, 'Csf.smf', C3, E3);
Mesh_WriteSMF(of, 'Brain.smf', C4, E4);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh_WriteSMF(of, name, Coord, Elem);
nnp = size(Coord,1); 
nel = size(Elem,1);
fid = fopen([of name], 'w');
fprintf(fid,'v %f %f %f \r\n',Coord(:,2:4)');
fprintf(fid,'t %d %d %d \r\n',Elem(:,2:4)');
fclose(fid);
