% nft_dsl_forward_model_generation() - Distributed source localization
% (dsl) forward model generation
%
% Usage:
%   >> nft_dsl_forward_model_generation(subject_name, session_name, of)
%
% Inputs:
%   subject_name - subject name as entered in NFT main window
%   session_name - session name as entered in NFT main window
%   of - output folder where all the outputs are saved
%
% Optional keywords:
%
%   cond : conductivity vector (default = [0.33 0.0132 1.79 0.33])
%   mesh_name : mesh name that will be loaded (default: [subject_name '.1.msh'])
%   sensor_name:  sensor name (default: [subject_name '_' session_name '.sensors'])
%   ss_name : sourcespace name (default: [subject_name '_sourcespace.dip'])
%   LFM_name :  LFM name (default: session_name_LFM)
%
% Author: Zeynep Akalin Acar, SCCN, 2021

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


function nft_dsl_forward_model_generation(subject_name, session_name, of, selection, mri, fs_dir)

cd(of)

NumberNodes = 80000; % default;

% selection of geodesic patch sizes
% selection = 2 -> 10, 6, 3 mm 
% selection = 3 -> 10 mm 
% selection = 4 -> 6 mm
% selection = 5 -> 3 mm

% Get Freesurfer meshes, run recon-all if needed
% test data
% cd /data/cta/zeynep/MRIs/140506JC/FS/surf
ofFS = fullfile(fs_dir, 'surf');
conf = nft_get_config;

if exist(ofFS, 'dir') == 0
    % fs_dir is not valid
    if ~isempty(fs_dir)
        % warn only if user passed something
        fprintf('No Freesurfer surface in fs_dir: %s\n',fs_dir);
    end
    disp('Running Freesurfer...'); pause(1);
    a = sprintf('%s -subject FS -sd "%s" -i "%s" -all', conf.freesurfer, of, mri);
    [status, result] = system(a);
    if status ~= 0; error('FreeSurfer:system','Failed to execute: %s',result); end
    disp('Freesurfer completed!'); pause(1);
    ofFS = fullfile(of, 'FS', 'surf');
end

% Co-register the NFT mesh with the Freesurfer mesh

cd(ofFS)

[vl, fl] = freesurfer_read_surf('lh.pial');
[vr, fr] = freesurfer_read_surf('rh.pial');

cd(of)

lof = length(of);
if of(lof) ~= filesep
    of(lof+1) = filesep;
end

disp('Co-registering the FS mesh with the NFT mesh...')
vertices = [vr; vl];
faces = [fr; fl+max(max(fr))];
[nf, nv] = reducepatch(faces, vertices, NumberNodes*2);
% [d, nv] = warping_distafterwarping([0 0 0 0 0 0], nv, nv);
nv = nv + 128;

[Cb, Eb] = ReadSMF('Brain.smf', 0, 0, 0, 1);
plotmesh(Eb(:,2:4), Cb(:,2:4)); view(0,90); hold
plot3(nv(:,1),nv(:,2),nv(:,3),'r.'); pause(1);

ma = mean(nv);
nv1 = nv - ones(length(nv),1) * ma;
rat = 0.95;
nv1 = rat * nv1;
nv1 = nv1 + ones(length(nv1),1) * ma;

nvt = [nv; nv1]; nft = [nf; nf+max(max(nf))];
Vfs = mesh2vol3(nvt, nft);
load([subject_name '_segments'])
Segm = new_segm(Segm, Vfs);
nsname = [subject_name 'FS'];
save([nsname '_segments.mat'],'Segm');

disp('Generating a new mesh...')
nft_mesh_generation(nsname, of, 4, 'Segm', Segm, 'lin_femmesh', 1)

ma = mean(nv);
F = nv - ones(length(nv),1) * ma;
rat = 0.95;
F = rat * F;
F = F + ones(length(F),1) * ma;

[Coord,Elem] = ReadSMF('Brain.smf',0,0,0,1); % Brain model
disp('Correcting the source space...');
[so2, k1,k2] = correct_source_space(F, Coord, Elem);

% load mesh configuration for path names
conf = nft_get_config;

%save ss_cor k2 k1 so2
disp('Running procmesh for correction...')
WriteSMF2('FSss.smf', so2, nf); % run showmesh once and save XXX

f=fopen(sprintf('%sStepSc.txt',of), 'w');
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

a = sprintf('"%s" -c "%sStepSc.txt" FSss.smf', conf.showmesh2, of);
[status, result] = system(a);
if status ~= 0; error('Mesh_Generation:system','Failed to execute: %s',result); end

[Css, Ess] = ReadSMF('ScS.smf',0,0,0,1); % FS sourcespace

save FSss_cor Css Ess   % sourcespace
disp('Calculating the area of nodes...');
[Ae, An] = area_of_nodes(Css,Ess);
save Node_area Ae An

% find node normals for Css Ess mesh
disp('Calculating node normals...')
Nn = NodeNormals(Css,Ess,1);
ss = [Css(:,2:4) Nn];

save([nsname '_ss.dip'],'ss','-ascii');


disp('Checking scalp sensor locations...')
sensor_file = [subject_name '_' session_name '.sensors'];
a = dir(sensor_file);
if size(a,1) > 0
    se = load(sensor_file,'-mat');
    handles.eloc = se.eloc;
    % check sensor locations
    [Cs,Es] = ReadSMF('Scalp.smf',0,0,0,1); % Brain model
    [F2, dmi] = warping_distmeshafterwarping([0 0 0 0 0 0], se.pnt, Cs, Es);
    se.pnt = F2;
    save(sensor_file, '-STRUCT', 'se')
else
    disp('Please co-register electrode locations with the head model')
end

disp('Cortical source space is saved!')




% generate geodesic patches


% selection = 2 -> 10, 6, 3 mm 
% selection = 3 -> 10 mm 
% selection = 4 -> 6 mm
% selection = 5 -> 3 mm

curr_dir = pwd;
cd(of)

global geodesic_library;                
geodesic_library = 'geodesic_debug';      %"release" is faster and "debug" does additional checks
rand('state', 0);                         %comment this statement if you want to produce random mesh every time

load FSss_cor
mesh = geodesic_new_mesh(Css(:,2:4),Ess(:,2:4));
algorithm = geodesic_new_algorithm(mesh, 'exact');      %initialize new geodesic algorithm

if selection == 2
    %gaussian patches
    disp('Calculating patches with 10 mm...')
    ss_g10 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 10); 
    save ss_g10 ss_g10
    disp('Calculating patches with 10 mm is done!')
    disp('Calculating patches with 6 mm..')
    ss_g6 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 6);
    save ss_g6 ss_g6
    disp('Calculating patches with 6 mm is done!')
    disp('Calculating patches with 3 mm...')
    ss_g3 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 3);
    save ss_g3 ss_g3
    disp('Calculating patches with 3 mm is done!')
elseif selection == 3
    disp('Calculating patches with 10 mm...')
    ss_g10 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 10); 
    save ss_g10 ss_g10
    disp('Calculating patches with 10 mm is done!')
elseif selection ==4
    disp('Calculating patches with 6 mm..')
    ss_g6 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 6);
    save ss_g6 ss_g6
    disp('Calculating patches with 6 mm is done!')
elseif selection ==5    
    disp('Calculating patches with 3 mm...')
    ss_g3 = ss_cortex_gaus(mesh, algorithm, Css(:,2:4), Ess(:,2:4), 3);
    save ss_g3 ss_g3
    disp('Calculating patches with 3 mm is done!')
end
    
cd(curr_dir)











function V = mesh2vol3(C, E)

vsize = floor(max(C)) + 10;
V = zeros(vsize, 'int8');
Elist = cell(vsize(3),1);

ne = size(E,1);

% put elements into Z buckets
for e = 1:ne
    Z = floor(C(E(e,:),3)) + 1;
    Z0 = min(Z);
    Z1 = max(Z) + 1;
    for z = Z0:Z1
        Elist{z} = [Elist{z} e];
    end
end

% Draw element-plane intersections for each z-plane
for z = 1:vsize(3)
    Ez = Elist{z};
    Zp = z - 0.5;
    for Zp = z-1:0.1:z
        for e = Ez
            % intersect element e with plane z = Zp
            ed = [ C(E(e,1),:) C(E(e,2),:)
                   C(E(e,2),:) C(E(e,3),:)
                   C(E(e,3),:) C(E(e,1),:)];
            pl = [];
            for i = 1:3
                z0 = ed(i,3);
                z1 = ed(i,6);
                if ((z0 < Zp && z1 > Zp) || (z0 > Zp && z1 < Zp))
                    p0 = ed(i,1:3);
                    p1 = ed(i,4:6);
                    dp = p1 - p0;
                    pi = p0 + dp * (Zp - z0)/dp(3);
                    pl = [pl; pi];
                end
                if (ed(i,3) == Zp)
                    pl = [pl; ed(i,1:3)];
                end
                if (ed(i,6) == Zp)
                    pl = [pl; ed(i,4:6)];
                end
            end
            np = size(pl,1);
            if (np > 1)
                for i = 1:np - 1
                    v1 = floor(pl(i, 1:2)) + 1;
                    v2 = floor(pl(i+1, 1:2)) + 1;
                    % draw line on z
                    dv = v2 - v1;
                    dvmax = max(abs(dv));
                    if (v1(1) > 0 && v1(2) > 0)
                        V(v1(1),v1(2),z) = 1;
                    end
                    for j = 1:dvmax
                        v = round(v1 + (dv * j /dvmax));
                        if (v(1) > 0 && v(2) > 0)
                            V(v(1),v(2),z) = 1;
                        end
                    end
                end
            end
        end
    end
end
[K,L,M]=size(V);
for i = 1:size(V,3)
    V(:,:,i) = imfill(V(:,:,i), 'holes');
end
for i = 1:size(V,2)
    V(:,i,:) = imfill(reshape(V(:,i,:),K,M), 'holes');
end
for i = 1:size(V,1)
    V(i,:,:) = imfill(reshape(V(i,:,:),L,M), 'holes');
end



    

function [Coord, Elem]=ReadSMF(name,x,y,z,sc)
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

% Vfs = mesh2vol2(v5, E1(:,2:4)); % freesurfer brain surface
function Segm = new_segm(Segm,Vfs);
Vbr0 = Segm.brainmask;
Vc0 = Segm.innerskullmask;
Vsk0 = Segm.outerskullmask;
Vsc0 = Segm.scalpmask;

s(1,:) = size(Vfs);
s(2,:) = size(Vbr0);
sz = max(s);

Vfs0 = zeros(sz, 'int8');
Vfs0(1:s(1,1), 1:s(1,2), 1:s(1,3)) = Vfs;

se1 = strel('ball',3, 3);
se = strel('ball', 3, 3, 0);
%se = strel('ball', 3, 3);

% new brain volume not intersecting with FS brain
Vfs2 = imclose3D(Vfs0, se1);
Vfs2 = imdilate3D(int8(Vfs2), se);
Vfs2 = Vfs2 - min(min(min(Vfs2)));
Vbr2 = Vfs2 | Vbr0; 
clear Vfs2

% new csf volume not intersecting with brain
Vbr3 = imdilate3D(int8(Vbr2), se);
Vbr3 = Vbr3 - min(min(min(Vbr3)));
Vc2 = Vc0 | Vbr3;
clear Vbr3

% new skull volume not intersecting with csf
Vc3 = imdilate3D(int8(Vc2), se);
Vc3 = Vc3 - min(min(min(Vc3)));
Vsk2 = Vsk0 | Vc3;
clear Vc3

% new skull volume not intersecting with csf
Vsk3 = imdilate3D(int8(Vsk2), se);
Vsk3 = Vsk3 - min(min(min(Vsk3)));
Vsc2 = Vsc0 | Vsk3;
clear Vsk3

% farklara bak
A = int8(Vbr0) - int8(Vbr2); 
clear Vbr0
fbr = sum(sum(sum(abs(A))));
A = int8(Vc0) - int8(Vc2);
clear Vc0
fc = sum(sum(sum(abs(A))));
A = int8(Vsk0) - int8(Vsk2);
clear Vsk0
fsk = sum(sum(sum(abs(A))));
A = int8(Vsc0) - int8(Vsc2);
clear Vsc0
fsc = sum(sum(sum(abs(A))));
[fbr fc fsk fsc]
% set the new volumes in Segm structure
Segm.innerskullmask = Vc2;
Segm.outerskullmask = Vsk2;
Segm.brainmask = Vbr2;
Segm.scalpmask = Vsc2;
clear Vsc2 Vbr2 Vsk2 Vc2
% Vfs = mesh2vol2(v5, E1(:,2:4)); % freesurfer brain surface

function [so2, k1,k2] = correct_source_space(so, C1, E1)

% so : sourcespace C1, E1 : mesh
% corrects source space
[N, M] = element_normals(C1, E1);
so2 = so;
[dim, inm] = CheckSourceSpace(so(:,1:3), C1, E1, 2);
k1 = find(inm == 0); %k1x = find(inm==1);
%if length(k1x) < length(k1)
%    k1 = k1x;
%end
no_intnodes = length(k1)

k2 = find(dim < 2); % find the nodes closer than 1mm
k2 = setdiff(k2, k1);
no_closenodes = length(k2)

if length(no_intnodes)>0
    
    for iter=1:2
    for i = 1:no_intnodes
        p1 = so2(k1(i),1:3);
        [dm, Pm, el, in] = DistMeshPoint2(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        % pn = Pm+2*nor;
        pn = Pm-2*N(el,:);
        [dm1, Pm1, el1, in1] = DistMeshPoint2(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k1(i)
        end
        so2(k1(i),1:3) = pn;
    end
    end
end

if length(no_closenodes)>0
    for i = 1:no_closenodes
        p1 = so2(k2(i),1:3);
        [dm, Pm, el, in] = DistMeshPoint2(p1, C1, E1);
        nor = (Pm-p1)/norm(Pm-p1);
        pn = Pm-nor*2;
        [dm1, Pm1, el1, in1] = DistMeshPoint2(pn, C1, E1);
        if in1 == 0
            disp('failed!')
            k2(i)
        end
        so2(k2(i),1:3) = pn;
    end
end

function [dim, inm] = CheckSourceSpace(so,C,E,thr)
% so is the source space
% C,E is the linear mesh
% dim gives the vector of distances
% inm gives the vector of inside sources (logical)

%hh = waitbar(0,'computing...');
M = size(so,1);
for i = 1 : M
 %   waitbar(i/M)
 if round(i/1000) == i/1000
     i
 end
    [dm, Pm, el, in] = DistMeshPoint2(so(i,:), C, E);
    dim(i) = dm;
    inm(i) = in;
end
%close(hh);
k = find(inm == 0);   % dipoles outside the mesh
l = find(dim < thr);    % dipoles closer to the mesh less than 1mm
m = setdiff(k, l);   % dipoles outside the mesh, closer to the mesh less than 1mm 
    

function [dm,Pm,el,in]=DistMeshPoint2(P,Coord,Elem);
% looks for if P is inside the mesh Coord, Elem or not
% Pm is the point on the mesh
% dm is the distance
% el is the element of Pm
% in = inside (bool)
% works for LINEAR MESH

% Coord: mesh vertex coordinates (first colum shows the indices)
% Elem: mesh connectivity (first colum shows the indices)


nnp=size(Coord,1);

r=5;
N=0;
Xo=Coord(:,2:4)-ones(nnp,1)*P;
Rad=sum(Xo.*Xo,2);
while length(N)<3
    r=r+5;
    N=find(Rad<r^2);
end

% find the neighbour elements of the closest nodes
E = ElementsOfTheNodes(Coord,Elem,N);


% find the intersection of the line PP1 with the elements of E
Pint=[]; dis=[];
for i=1:length(E)
   Pa=Coord(Elem(E(i),2),2:4);
   Pb=Coord(Elem(E(i),3),2:4);
   Pc=Coord(Elem(E(i),4),2:4);
   [D,Pp]=DistTrianglePoint2(Pa,Pb,Pc,P);
   dis(i)=D;
   Pint(i,:)=Pp;
end
[j,k]=min(abs(dis));
Pm=Pint(k,:);
dm=abs(dis(k));
el=E(k);

N2 = Elem(el,2:4);
el2 = ElementsOfTheNodes(Coord,Elem,N2);

% check if P is inside or outside

% find the normal vector of the element
v1=Coord(Elem(el2,3),2:4)-Coord(Elem(el2,2),2:4);
v2=Coord(Elem(el2,4),2:4)-Coord(Elem(el2,3),2:4);
n=cross(v1, v2); n2=mean(n);
Norm = n2/norm(n2);

% if the vector PPm.Norm > 0 inside, < 0 outside
if dot(Pm-P, Norm) > 0 
    in = 1; % P is inside the mesh
else
    in = 0;
end



function [D,Pp]=DistTrianglePoint2(Pa,Pb,Pc,Px)
% finds the minimum distance of a point Px to triangle Pa, Pb, Pc
% difference from DistTrianglePoint 
% doesn't look if the projection of the point is in the triangle or on the edge 
% of the triangle otherwise MinD is 1000

% find the minimum distance of the point with the 
% plane which is formed by the triangle
% find the normal of the plane
eps=1e-4;
v1=Pa-Pb;
v2=Pc-Pb;
n=[v1(2)*v2(3)-v2(2)*v1(3) v2(1)*v1(3)-v1(1)*v2(3) v1(1)*v2(2)-v1(2)*v2(1)];
n=n/norm(n);
d=-n(1)*Pa(1)-n(2)*Pa(2)-n(3)*Pa(3);
% the plane equation is n(1)*x+n(2)*y+n(3)*z+d=0
% distance of the point to the plane is D
D=(n(1)*Px(1)+n(2)*Px(2)+n(3)*Px(3)+d)/sqrt(n(1)^2+n(2)^2+n(3)^2);
% point on the plane
Pp=Px-D*n;
% check if the point is on the triangle
% Determine whether or not the intersection point is bounded by pa,pb,pc 
Pa1=Pa-Pp;
normPa1=norm(Pa1);
if normPa1>eps
   % normalize the unit vectors
   Pa1=Pa1/normPa1;  
end
Pa2 = Pb - Pp;
normPa2=norm(Pa2);
if normPa2>eps
   Pa2=Pa2/normPa2; 
end
Pa3 = Pc - Pp;
normPa3=norm(Pa3);
if normPa3>eps
   Pa3=Pa3/normPa3;
end
%the angles are 
a1 = acos(Pa1(1)*Pa2(1) + Pa1(2)*Pa2(2) + Pa1(3)*Pa2(3));
a2 = acos(Pa2(1)*Pa3(1) + Pa2(2)*Pa3(2) + Pa2(3)*Pa3(3));
a3 = acos(Pa3(1)*Pa1(1) + Pa3(2)*Pa1(2) + Pa3(3)*Pa1(3));

total = a1+a2+a3;
% if total is 2*pi then the point is in the triangle or on the edges
if (abs(total - 2*pi) < eps)
   MinD=abs(D);
else
    % find the distance of Px to the edges
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pb);
   Dix(1)=Di; Ppi(1,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pa,Pc);
   Dix(2)=Di; Ppi(2,:)=Ppl;
   [Di,Ppl]=DistPointLineSegment(Px,Pc,Pb);
   Dix(3)=Di; Ppi(3,:)=Ppl;
   [u,v]=min(Dix);
   MinD=u;
   Pp=Ppi(v,:);
end
D=MinD;

function [D,Pp]=DistPointLineSegment(P,P1,P2);
% Finds distance between the point P and the line segment P1P2

v1=P2-P1;
n1 = v1/norm(v1);

v2 = P-P1;

% projection (component) of v2 along v1 (n1)
vp = dot(n1, v2);

% make sure it is in the range
if (vp < 0)
    vp = 0;
elseif (vp > norm(v1))
    vp = norm(v1);
end


% closest point on edge
Pp = P1 + vp * n1;

D = norm(P-Pp);


function E=ElementsOfTheNodes(Coord,Elem,A);
% A is the list of the nodes
% E is the list of the elements that have nodes in A
E=[];
nop=size(Elem,2);
if nop==4
   for i=1:length(A)
   	n1=find(Elem(:,2)==A(i));
	   n2=find(Elem(:,3)==A(i));
   	n3=find(Elem(:,4)==A(i));
	   n4=union(n1,n2);
   	n5=union(n3,n4);
	   E=union(E,n5);
      clear n1 n2 n3 n4 n5
   end
elseif nop==7
   for i=1:length(A)
   	n1=find(Elem(:,2)==A(i));
	   n2=find(Elem(:,3)==A(i));
      n3=find(Elem(:,4)==A(i));
      n4=find(Elem(:,5)==A(i));
   	n5=find(Elem(:,6)==A(i));
   	n6=find(Elem(:,7)==A(i));
	   n7=union(n1,n2);
      n8=union(n7,n3);
      n9=union(n8,n4);
      n10=union(n9,n5);
      n11=union(n10,n6);
	   E=union(E,n11);
      clear n1 n2 n3 n4 n5 n6 n7 n8 n9 n10 n11
   end
end


function [N, M] = element_normals(Coord, Elem)

if size(Elem,2)==7

n1 = Coord(Elem(:,2),2:4);
n2 = Coord(Elem(:,4),2:4);
n3 = Coord(Elem(:,6),2:4);
elseif size(Elem,2)==4
n1 = Coord(Elem(:,2),2:4);
n2 = Coord(Elem(:,3),2:4);
n3 = Coord(Elem(:,4),2:4);
end    

M = (n1 + n2 + n3) / 3;

n1 = n1 - n3;
n2 = n2 - n3;

N = cross(n1, n2);


for i=1:size(Elem,1)
    N(i,:)=N(i,:)/norm(N(i,:));
end



function WriteSMF2(name,Coord,Elem);
nnp=size(Coord,1); nel=size(Elem,1);
fid=fopen(name, 'w');
fprintf(fid,'v %f %f %f \r\n',Coord');
fprintf(fid,'t %d %d %d \r\n',Elem');
fclose(fid);


function [Ae, An] = area_of_nodes(C,E);
% Ae : area of elements
% An : area of nodes

Ne = length(E);

for i = 1:Ne
    X = C(E(i,2:4),2:4);
  
    AB = X(1,:)-X(2,:);
    AC = X(1,:)-X(3,:);
    Ae(i) = 1/2 * norm(cross(AB,AC));
end

Nn = length(C);
for i = 1:Nn
    e1 = ElementsOfTheNodes(C,E,i);
    An(i) = mean(Ae(e1));
end
   


function Nn = NodeNormals(Coord,Elem,ccw);
% if ccw=1, then the mesh is CCW
% if ccw=0, mesh is CW

% calculates normals of nodes for a linear mesh

En = ElemNormal(Elem, Coord);
nr = length(Coord(:,1));

hh = waitbar(0,'computing node normals...');
for i = 1:nr
    waitbar(i/nr)

    %ns = findNeigNodes(Coord, Elem, i);
    ns = ElementsOfTheNodes(Coord,Elem,i);
    if length(ns)<3
        nn = findNeigNodes(Coord,Elem,i);
        ns = ElementsOfTheNodes(Coord,Elem,nn);
        if length(ns)<3
            K = Coord(:,2:4)-ones(length(Coord),1)*Coord(i,2:4);
            M = sqrt(sum(K.*K,2));
            [Y,I] = sort(M);
            j=1;
            while length(ns)<3
                ns = ElementsOfTheNodes(Coord,Elem,I(j));
                j=j+1;
            end
        end
    end
    Nn(i,:) = mean(En(ns,:),1);
    Nn(i,:) = Nn(i,:)/norm(Nn(i,:));
end
close(hh)
if ccw==0
    Nn = -Nn;
end

% correct if there are NaNs in the calculation of normals

w = isnan(Nn);
[I,J] = find(w==1);
if length(I)>0
    for i = 1:length(I)
        iw = [];
        nn = I(i);
        while length(iw) == 0
            nn = findNeigNodes(Coord,Elem,nn);
            Nni = Nn(nn,:);
            w =~ isnan(Nni);
            [iw, jw] = find(w==1);
        end
        Nn(I(i),:) = Nni(iw(1),:);
    end
end


function Norm=ElemNormal(Elem,Coord);

ne = length(Elem(:,1));
nr = length(Coord(:,1));

for k = 1:ne
   v1 = Coord(Elem(k,3),2:4)-Coord(Elem(k,2),2:4);
   v2 = Coord(Elem(k,4),2:4)-Coord(Elem(k,3),2:4);
   n = cross(v1, v2);
   Norm(k,1:3) = n/norm(n);
end

function N=findNeigNodes(Coord,Elem,n);
% find the neighbour nodes of the nth node of the mesh

nop=size(Elem,2);
% find the neighbouring elements
if nop==4
    E=[];
    for i=1:length(n)
        E2=find(Elem(:,2)==n(i) | Elem(:,3)==n(i) | Elem(:,4)==n(i));
        E = union(E,E2);
    end

    f=Elem(E,2:4);
    N=unique(f(:));
    
elseif nop==7
    E=[];
    for i=1:length(n)
        E2=find(Elem(:,2)==n(i) | Elem(:,3)==n(i) | Elem(:,4)==n(i) | Elem(:,5)==n(i) | Elem(:,6)==n(i) | Elem(:,7)==n(i));
        E = union(E,E2);
    end

    f=Elem(E,2:7);
    N=unique(f(:));
end
    

