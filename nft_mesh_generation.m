% nft_mesh_generation() - Mesh generation from segmentation
%
% Usage:
%         >> nft_mesh_generation(subject_name, of, nl,'Key1',Value1,...);
%
% Inputs:
%
%   subject_name : subject name as in main NFT window 
%   of : output folder  
%   nl : number of layers
%
% Optional keywords:
%
%   Segm : Segm structure obtained in NFT segmentation.
%   mesh_name : string to set a different mesh name (default: subject_name)
%   nnpl: integer value for number of nodes per layer (default = 7000)
%   lmr :  [0 or 1] for local mesh refinement (default = 1)
%   ratio_lmr :  ratio for lmr (default = 2.1)
%   lin_femmesh : [0 or 1] for linear FEM mesh generation (default = 0)
%   quad_femmesh : [0 or 1] for quadratic FEM mesh generation (default = 0)

function nft_mesh_generation(subject_name, of, nl, varargin)

% default values
mesh_name = subject_name;
nnpl = 7000;
lmr = 1;
ratio_lmr = 2.1;
lin_femmesh = 0;
quad_femmesh = 0;

current_folder = pwd;
lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end
cd(of)

file = [of subject_name '_segments'];
load(file)


for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};

      if ~isstr(Keyword)
         fprintf('keywords must be strings')
         return
      end

      if strcmp(Keyword,'Segm')
         if isstr(Value)
            fprintf('Segm must be a segmentation structure');
            return
         else
            Segm = Value;
         end
      elseif strcmp(Keyword,'mesh_name')
         if ~isstr(Value)
            fprintf('mesh_name must be a string');
            return
         else
             mesh_name = Value;
         end
      elseif strcmp(Keyword,'nnpl')
         if isstr(Value)
            fprintf('nnpl must be a positive integer');
            return
         else
             nnpl = Value;
         end
      elseif strcmp(Keyword,'lmr')
         if isstr(Value)
            fprintf('lmr must be 0 or 1');
            return
         else
             lmr = Value;
         end
      elseif strcmp(Keyword,'ratio_lmr')
         if isstr(Value)
            fprintf('ratio_lmr must be a positive real number');
            return
         else
             ratio_lmr = Value;
         end
      elseif strcmp(Keyword,'lin_femmesh')
         if isstr(Value)
            fprintf('lin_femmesh must be 0 or 1');
            return
         else
             lin_femmesh = Value;
         end
      elseif strcmp(Keyword,'quad_femmesh')
         if isstr(Value)
            fprintf('quad_femmesh must be 0 or 1');
            return
         else
             quad_femmesh = Value;
         end
      end
      
end


Sca = Segm.scalpmask;
Bra = Segm.brainmask;
Csf = Segm.innerskullmask;
Skull = Segm.outerskullmask;

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end

transform = size(Sca) / 2;
save([of 'transform'], 'transform', '-ascii'); 


disp('Saving volumes in raw...'); pause(0.5)
Mesh_writeraw(Sca, [of 'Scalp']);
Mesh_writeraw(Bra, [of 'Brain']);
Mesh_writeraw(Csf, [of 'Csf']);
Mesh_writeraw(Skull, [of 'Skull']);

[K,L,M] = size(Sca);

% load mesh configuration for path names
conf = nft_get_config;

tis_type = cell(1,4);
tis_type{1} = 'Scalp';
tis_type{2} = 'Brain';
tis_type{3} = 'Csf';
tis_type{4} = 'Skull';

disp('Triangulating volumes...'); pause(0.5)
hh = waitbar(0,'Triangulating volumes...');
for tis = 1:4
    waitbar(tis/4);
    tt = char(tis_type(tis));
    a= sprintf('"%s" -t 10 -dr1 "%s%s.raw" %d %d %d -f "%s%s.asc" -ot', conf.asc, of,tt, K, L, M,of,tt);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
end
close(hh);

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

NumberNodes = nnpl;
Quad = 0;

csi = 300000; i = 1;
while NumberNodes < csi
    i = i+1;
    csi(i) = round(csi(i-1) / 1.5);
end
csi(i) = NumberNodes;

nsteps = length(csi); % number of coarsening steps
ntis = 4;  % number of tissues (4)

%%%%%%%%%%%%%%%%%% Coarsening and correcting...
hh = waitbar(0,'Coarsening and correcting...');
for tis = 1:ntis
    tt = char(tis_type(tis));
        
    a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.asc"',conf.showmesh,of,of,tt);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
    
     for iter = 1:nsteps
         waitbar(((tis-1)*nsteps + iter) / (nsteps*ntis));
         a = sprintf('"%s" -c 0.5 -m 5 -o "%s%s.smf" -t %d "%sScS.smf"', conf.qslim, of, tt, csi(iter), of);
         [status, result] = system(a);
         if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
 
         a = sprintf('"%s" -c "%sStepSc.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
         [status, result] = system(a);
         if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
         if Quad & iter == 6
             copyfile([of 'ScS.smf'], [of tt 'f.smf'])
         end

     end
     a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of, of, tt);
     [status, result] = system(a);
     if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
     
     movefile([of 'ScS.smf'], [of tt '.smf']); % XXX yeni
end
close(hh)


param.no_layers = nl;
param.linear = 1;
disp('Final mesh correction...'); pause(0.5)

% Final correction with Matlab functions
mesh_final_correction_v2(of, nl);

%if Quad == 0
    % Local mesh refinement
    if lmr == 1
        disp('Local mesh refinement...'); pause(0.5)
        mesh_local_refinement(of, nl, ratio_lmr);
        handles.param.lmr = ratio_lmr;
    end
%end
mesh_final_correction_v2(of, nl);


% generate quadratic meshes
if Quad == 1
    handles.param.linear = 0;
    for tis = 1:ntis
        tt = char(tis_type(tis));
        a = sprintf('"%s" "%s%s.smf" "%s%sf.smf" %s%sq', conf.quad, of,tt,of,tt,of,tt);
        [status, result] = system(a);
        if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
    end
end
       

% delete unnecessary files (.raw, .asc, Scs.smf, and StepSc.txt)
%delete([of 'ScS.smf']);
delete([of 'StepSc.txt']);
delete([of 'StepSc2.txt']);
for tis = 1:ntis
    tt = char(tis_type(tis));
    delete([of tt '.raw']);
    delete([of tt '.asc']);
end


disp('Mesh generation done! Saving...'); pause(0.5)
mesh_read_write(of, mesh_name, nl, Quad); % generate nl-layer head model
param.mesh_name = mesh_name;
parameters = param;
save([of mesh_name '_mesh'], '-STRUCT', 'parameters')

if lin_femmesh | quad_femmesh
    generate_FEM_mesh(mesh_name, of, quad_femmesh);
end

cd(current_folder)

function Mesh_writeraw(file, fn);
A = single(file);
norm = max(max(max(A)));
A = A * 20 / norm;
filename = [fn '.raw'];
f=fopen(filename, 'w+');
fwrite(f, A, 'uint8');
fclose(f);
clear A norm f 


function generate_FEM_mesh(mesh_name, of, quad)


conf = nft_get_config;
fn = [of mesh_name '.bei'];
 % read BEM mesh
mesh = bem_load_mesh([of mesh_name]);
R = mesh_find_regions(mesh);
Coord(:,2:4) = mesh.coord;
Elem(:,2:4) = mesh.elem;
Coord(:,1) = [1:length(Coord)]';
Elem(:,1) = [1:length(Elem)]';
nl= mesh.num_boundaries;

% convert to Tetgen compatible input
fn = [of mesh_name '.smesh'];
WriteSMESH(fn, Coord, Elem, R);

% convert into metu-fem mesh
if nl == 4
    cnd_str = '1=C1 2=C2 3=C3 4=C4';
elseif nl == 3
    cnd_str = '1=C1 2=C2 3=C3';
end

if quad
    % call tetgen - generate linear tetrahedral mesh
    % a = sprintf('"%s" -pq1.4a12A "%s"', conf.tetgen, fn);
    a = sprintf('"%s" -pq1.4a120A "%s"', conf.tetgen, fn);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
    
    fn = [of mesh_name '.1'];
    a = sprintf('"%s" "%s" %s', conf.tetgen2msh, fn, cnd_str);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
    
    % convert into quadratic mesh
    fn_q = [of mesh_name '.q'];
    a = sprintf('"%s" -o "%s".msh "%s".msh', conf.lin2quad, fn_q, fn);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
else
    a = sprintf('"%s" -pq1.4a5A "%s"', conf.tetgen, fn);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
    
    fn = [of mesh_name '.1'];
    a = sprintf('"%s" "%s" %s', conf.tetgen2msh, fn, cnd_str);
    [status, result] = system(a);
    if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end
end




function WriteSMESH(name,Coord,Elem,Regions)
% Saves the mesh in .SMESH format for use with Tetgen
% The format is described in:
% http://tetgen.berlios.de/fformats.html
% The .smesh format is slightly simpler than more general .poly
% The Regions parameter must specify a point inside each region
% TetGen will mark output tetrahedra using these region markers

nnp=size(Coord,1);
nel=size(Elem,1);
if ~isempty(Regions)
	nreg = size(Regions,1);
	Reg = zeros(nreg, 6);
	Reg(:,2:4) = Regions;
	Reg(:,1) = 1:nreg;
	Reg(:,5) = 1:nreg;
	Reg(:,6) = -1;
else
	nreg = 0;
	Reg = [];
end

% make sure Node indices are correct
Coord(:,1) = 1:nnp;

fid=fopen(name, 'w');

fprintf(fid, '# Part 1 - node list\n');
fprintf(fid, '# node count, 3 dim, no attr, no boundary \n');
fprintf(fid, '%d 3 0 0\n', nnp);
fprintf(fid, '# Node index, node coordinates\n');
fprintf(fid, '%f %f %f %f\n',Coord');

fprintf(fid, '# Part 2 - facet list\n');
fprintf(fid, '# facet count, no boundary marker\n');
fprintf(fid, '%d\n',nel);
fprintf(fid, '# facets\n');
fprintf(fid, '3 %d %d %d\n', Elem(:,2:4)');

fprintf(fid, '# Part 3 - hole list\n');
fprintf(fid, '0            # no hole\n');

fprintf(fid, '# Part 4 - region list\n');
fprintf(fid, '%d\n', nreg);

if nreg > 0
	fprintf(fid, '%d %f %f %f %d %d\n', Reg');
end

fclose(fid);
