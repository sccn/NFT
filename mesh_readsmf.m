% mesh_readsmf() - Reads meshes in smf format.
%
% Usage:
%   >> [Coord, Elem] = mesh_readsmf(name,x,y,z,sc);
%
% Inputs:
%   name - name of the mesh
%   x,y,z - translation
%   sc - scale
%
% Outputs:
%   Coord - Coordinate matrix
%   Elem - Connectivity matrix
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

function [Coord, Elem] = mesh_readsmf(name,x,y,z,sc);
fid=fopen(name, 'r');
nnp=0; nel=0;
line=1;
while line~=-1
   line=fgets(fid);
	if line(1)=='v';
   	nnp=nnp+1;
	   [A,count,ERRMSG,NEXTINDEX] = sscanf(line,'%c %f %f %f',4);
   	Coord(nnp,1)=nnp;
	   Coord(nnp,2:4)=A(2:4)';
	elseif (line(1)=='t')|(line(1)=='f');
   	nel=nel+1;
	   [A,count,ERRMSG,NEXTINDEX] = sscanf(line,'%c %d %d %d',4);
	  	Elem(nel,1)=nel;
   	Elem(nel,2:4)=A(2:4)';
   end
end
fclose(fid);

Coord(:,4)=Coord(:,4)/sc+z;
Coord(:,2)=Coord(:,2)+x;
Coord(:,3)=Coord(:,3)+y;

   