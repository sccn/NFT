% utilmesh_refine_correct() - Locally refines the meshes where the edge length 
% is larger than distance between the meshes with a given ratio. 
%
% Usage:
%   >> [C1,E1,C2,E2] = utilmesh_refine_correct(C1, E1, C2, E2, ratio_lmr);
%
% Inputs:
%   C1, C2 - coordinate matrices of input meshes
%   E1, E2 - connectivity matrices of input meshes
%   ratio_lmr - ratio of edge length to distance between meshes
%
% Outputs:
%   C1, C2 - coordinate matrices of output meshes
%   E1, E2 - connectivity matrices of output meshes
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

function [C1,E1,C2,E2] = utilmesh_refine_correct(C1, E1, C2, E2, ratio_lmr,of)

[ind_C1, ind_C2] = find_close_regions(C1, E1, C2, E2, ratio_lmr);

% load mesh configuration for path names
conf = nft_get_config;

iter = 0;
diff_length_ind_C1 = 10;
diff_length_ind_C2 = 10;
while (length(ind_C1)>0 | length(ind_C2)>0) & ((diff_length_ind_C1 > 0 | diff_length_ind_C2 > 0) & (iter < 5))
    iter = iter + 1;

    len_ind_C1 = length(ind_C1);
    len_ind_C2 = length(ind_C2);
    
    el1 = ElementsOfTheNodes(C1, E1, ind_C1);
    el2 = ElementsOfTheNodes(C2, E2, ind_C2);
    
  % refine first mesh  
    [C1, E1] = Local_mesh_refine(C1, E1, el1);
    [C2, E2] = Local_mesh_refine(C2, E2, el2);
    
    Mesh_WriteSMF(of, 'temp.smf', C1, E1);
    a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
    [status, result] = system(a);
    if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
    movefile([of 'ScS.smf'], [of  'temp1.smf'])
    [C1,E1] = mesh_readsmf([of 'temp1.smf'],0,0,0,1); 
     
    Mesh_WriteSMF(of, 'temp.smf', C2, E2);
    a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
    [status, result] = system(a);
    if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
    movefile([of 'ScS.smf'], [of  'temp2.smf'])
    [C2,E2] = mesh_readsmf([of 'temp2.smf'],0,0,0,1); 

    [ind_C1, ind_C2] = find_close_regions(C1, E1, C2, E2, ratio_lmr);
    
    diff_length_ind_C1 = len_ind_C1 - length(ind_C1);
    diff_length_ind_C2 = len_ind_C2 - length(ind_C2);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mesh_WriteSMF(of, name, Coord, Elem);
nnp = size(Coord,1); 
nel = size(Elem,1);
fid = fopen([of name], 'w');
fprintf(fid,'v %f %f %f \r\n',Coord(:,2:4)');
fprintf(fid,'t %d %d %d \r\n',Elem(:,2:4)');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ind_C1, ind_C2] = find_close_regions(C1, E1, C2, E2, ratio_lmr);

% for triangular meshes
% ind_C1 is the index of the nodes 
% ratio_lmr : ratio of maximum edge length to distance of two meshes in
% that region

ratio_lmr = 1/ratio_lmr; 

Nc1 = length(C1); % number of coordinates for the first mesh
Ne1 = length(E1); % number of elements for the first mesh

Nc2 = length(C2); % number of coordinates for the second mesh
Ne2 = length(E2); % number of elements for the second mesh

ind_C1 = [];
ind_C2 = [];
for i = 1:Nc1
    p1 = C1(i,2:4);
    % find the edges of the edges connected to p1
    [e1 j k] = find(E1(:,2:4)==i); % e1: elements of the ith node
    nn = E1(e1,2:4); nn = nn(:); 
    nn = unique(nn); nn = setdiff(nn,i); % nodes connected to ith node
    
    for j = 1:length(nn)
        el(j) = norm(p1 - C1(nn(j),2:4)); % edge length
    end

    mean_el = mean(el);
    M = C2(:,2:4) - ones(Nc2,1)*p1;
    K = sqrt(sum(M.*M,2)); % distance of all nodes of C2 to p1
    ind_small = find(K < ratio_lmr * mean_el);
    
    if ~isempty(ind_small)
        ind_C1 = [ind_C1; i];
        ind_C2 = [ind_C2; ind_small];
    end
end
ind_C2 = unique(ind_C2);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function E=ElementsOfTheNodes(Coord,Elem,A);
% A is the list of the nodes
% E is the list of the elements that have nodes in A
E=[];
for i=1:length(A)
   n1=find(Elem(:,2)==A(i));
   n2=find(Elem(:,3)==A(i));
   n3=find(Elem(:,4)==A(i));
   n4=union(n1,n2);
   n5=union(n3,n4);
   E=union(E,n5);
   clear n1 n2 n3 n4 n5
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [C,E] = Local_mesh_refine(C,E,El);

% C, E is a triangular mesh
% El is the list of elements that will be divided
% it divides the elements into 3 by placing a node in the center

Nc = length(C); % number of coordinates
Ne = length(E); % number of elements

ne = length(El); % number of elements that will be divided
for i=1:ne
    
    nn = Nc+i; % new coord number that will be added
    ec = E(El(i),2:4);
    e1 = [ec(1) ec(2) nn];
    e2 = [ec(1) nn ec(3)];
    e3 = [nn ec(2) ec(3)];
    E(El(i),2:4) = e1;
    E(Ne+1,:) = [Ne+1 e2];
    E(Ne+2,:) = [Ne+2 e3];
    new_coord = mean(C(ec,2:4));
    C(Nc+i,:) = [Nc+i new_coord];
%    Nc = length(C); % number of coordinates
    Ne = length(E); % number of elements
    
end

