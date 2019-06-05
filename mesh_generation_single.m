function mesh = mesh_generation_single(volume, of, tis_type, NumberElem)

% mesh generation for a single layer
% volume : input volume
% of : output folder
% tis_type: tissue type (scalp, skull, csf, brain)
% NumberElem : number of elements for the resulting mesh

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


lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end

Mesh_writeraw(volume, [of tis_type]);
[K,L,M] = size(volume);

% load mesh configuration for path names
conf = nft_get_config;

disp('Triangulating...')    
tt = char(tis_type);
a= sprintf('"%s" -t 10 -dr1 "%s%s.raw" %d %d %d -f "%s%s.asc" -ot', conf.asc, of,tt, K, L, M,of,tt);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

% generate a file for StepSc.txt for coarsening and smoothing
f=fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'correct 5\n');
fprintf(f, 'smooth 1\n');
%fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
%fprintf(f, 'split intersect\n'); % XXX yeni
%fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

% generate a file for StepSc2.txt for final improvement
f=fopen(sprintf('%sStepSc2.txt',of), 'w');
fprintf(f, 'correct 2\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'fill holes\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'fill holes\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);


csi = 300000; i = 1;
while NumberElem < csi
    i = i+1;
    csi(i) = round(csi(i-1) / 1.5);
end
csi(i) = NumberElem;

nsteps = length(csi); % number of coarsening steps

%%%%%%%%%%%%%%%%%% Coarsening and correcting...
tt = char(tis_type);

disp('Coarsening and correcting...')
a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.asc"',conf.showmesh,of,of,tt);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
   for iter = 1:nsteps
      a = sprintf('"%s" -c 0.5 -m 5 -o "%s%s.smf" -t %d "%sScS.smf"', conf.qslim, of, tt, csi(iter), of);
      [status, result] = system(a);
      if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
      a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
      [status, result] = system(a);
      if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
   end
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
copyfile([of 'ScS.smf'], [of tt '.smf'])

[C1,E1] = mesh_readsmf([of tis_type '.smf'],0,0,0,1);   
mesh.coord = C1(:,2:4);
mesh.elem = E1(:,2:4);

% delete unnecessary files (.raw, .asc, Scs.smf, and StepSc.txt)
delete([of 'ScS.smf']);
delete([of 'StepSc.txt']);
delete([of 'StepSc2.txt']);
tt = char(tis_type);
delete([of tt '.raw']);
delete([of tt '.asc']);



function Mesh_writeraw(file, fn);
A = single(file);
norm = max(max(max(A)));
A = A * 20 / norm;
filename = [fn '.raw'];
f=fopen(filename, 'w+');
fwrite(f, A, 'uint8');
fclose(f);
clear A norm f 


