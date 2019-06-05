% mesh_read_write() - Reads Scalp, Skull, Csf and Brain meshes in .smf 
% format in a given folder and writes the total head mesh  bec, bee, bei
% format.
%
% Usage:
%   >> mesh_read_write(of, mesh_name, nl);
%
% Inputs:
%   of - mesh folder
%   mesh_name - name of the mesh that will be saved (in bec, bee, andbei
%   format)
%   nl - number of layers (3 or 4)
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

function mesh_read_write(of, mesh_name, nl, Quad);

% call from the folder where the meshes are
% reads smf meshes and converts it to .bee, .bec, .bei format
% nl is the number of layers ,3 or 4

if not(nl==3|nl==4)
    error('number of layers should be entered as 3 or 4');
end

bec_name=[mesh_name '.bec'];
bee_name=[mesh_name '.bee'];
bei_name=[mesh_name '.bei'];

if Quad == 1
    fn = [of 'Scalpq.bec']; C1 = load(fn);
    fn = [of 'Scalpq.bee']; E1 = load(fn);
    fn = [of 'Skullq.bec']; C2 = load(fn);
    fn = [of 'Skullq.bee']; E2 = load(fn);
    fn = [of 'Csfq.bec'];   C3 = load(fn);
    fn = [of 'Csfq.bee'];   E3 = load(fn);
    fn = [of 'Brainq.bec']; C4 = load(fn);
    fn = [of 'Brainq.bee']; E4 = load(fn);
else
    [C1,E1] = mesh_readsmf([of 'Scalp.smf'],0,0,0,1);       % subject's scalp mesh
    [C2,E2] = mesh_readsmf([of 'Skull.smf'], 0, 0, 0, 1); % subject's skull mesh
    [C3,E3] = mesh_readsmf([of 'Csf.smf'], 0, 0, 0, 1);   % subject's CSF mesh
    [C4,E4] = mesh_readsmf([of 'Brain.smf'], 0, 0, 0, 1);   % subject's brain mesh
end

% save in the .bei, .bee, .bec format
[Coord,Elem]=utilbem_add_mesh(C1,E1,C2,E2);
[Coord,Elem]=utilbem_add_mesh(Coord,Elem,C3,E3);
if nl==4
    [Coord,Elem]=utilbem_add_mesh(Coord,Elem,C4,E4);
end

writeMesh(of, bec_name, bee_name, Coord, Elem);
info(1,1) = nl;
info(1,2) = size(Elem,1);
info(1,3) = size(Coord,1);
info(1,4) = size(Elem,2)-1;
info(2,1) = 1;
info(3,1) = 2;
info(4,1) = 3;
info(5,1) = 4;
info(2,2) = size(E1,1);
info(3,2) = size(E2,1);
info(4,2) = size(E3,1);
info(5,2) = size(E4,1);
info(2,3:4) = [1 0];
info(3,3:4) = [2 1];
info(4,3:4) = [3 2];
info(5,3:4) = [4 3];

if nl==3
    info = info(1:4,:);
end

fid=fopen([of bei_name], 'w');
fprintf(fid,'%d %d %d %d\r\n',info');
fclose(fid);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function writeMesh(of, nameC, nameE, Coord, Elem);
% of is the folder
fid=fopen([of nameC],'w');
fprintf(fid,'%5.15f %5.15f %5.15f %5.15f\r\n',Coord');
fclose(fid);

noe=size(Elem,2);
if noe==7
    fid=fopen([of nameE],'w');
    fprintf(fid,'%d %d %d %d %d %d %d\r\n',Elem');
    fclose(fid);
elseif noe==4
    fid=fopen([of nameE],'w');
    fprintf(fid,'%d %d %d %d\r\n',Elem');
    fclose(fid);
end
