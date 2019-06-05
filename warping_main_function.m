% warping_main_function() - Main warping function
%
% Usage:
%   >> [Ptm, ind, Cscalp_w,Cskull_w,CCSF_w,W,A,e,LMm2] =
%   warping_main_function(Cscalp, Escalp, Cskull, CCSF,LMm, Fm, Fd, pos);
%
% Inputs:
%   Cscalp, Escalp - scalp mesh
%   Cskull         - coordinates of skull mesh
%   CCSF           - coordinates of CSF mesh
%   LMm            - landmarks on the template mesh
%   Fm             - fiducials of the template mesh
%   Fd             - fiducials from the digitizer data
%   pos            - electrode locations
%
% Outputs:
%   Ptm      - electrode locations on the mesh
%   ind      - index of the electrodes on the mesh
%   Cscalp_w - warped scalp mesh coordinates
%   Cskull_w - warped skull mesh coordinates
%   CCSF_w   - warped CSF mesh coordinates
%   W, A, e  - warping transform parameters
%   LMm2     - warped landmarks 
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

function [Ptm, ind, Cscalp_w,Cskull_w,CCSF_w,Cbrain_w, W,A,e,LMm2,back,Escalp,Eskull,ECSF,Ebrain] = warping_main_function(of,Cscalp, Escalp, Eskull, ECSF, Ebrain,Cskull, CCSF, Cbrain, LMm, Fm, Fd, pos, index_kdm);

% index_kdm is the index of the electrodes that are close to the scalp, to
% do the warping

% find warping parameters
[W,A,e,LMm2, Pt, back] = find_warping(Cscalp, Escalp, LMm, Fm, Fd, pos,index_kdm);

% warp the mesh
Cscalp_w = warped_mesh(Cscalp,A,W,LMm2);
Cskull_w = warped_mesh(Cskull,A,W,LMm2);
CCSF_w = warped_mesh(CCSF,A,W,LMm2);
Cbrain_w = warped_mesh(Cbrain,A,W,LMm2);

%%%%%%%

% generate a file for StepSc2.txt for final improvement
f=fopen(sprintf('%sStepSc2.txt',of), 'w');
fprintf(f, 'correct 2\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 5\n');
fprintf(f, 'improve 2 0.1 0.05\n'); % 2=# of iter, 0.1=elem aspect ratio, 3edge<0.05*mean edge length=>delete
fprintf(f, 'improve 2 0.3 0.2\n');
fprintf(f, 'correct 5\n');
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'split intersect\n'); % XXX yeni
fprintf(f, 'correct 2\n');       % XXX yeni
fprintf(f, 'prune all\n');         % XXX yeni
fprintf(f, 'save %sScS.smf\n',of);
fprintf(f, 'quit\n');
fclose(f);

conf = nft_get_config;
Mesh_WriteSMF(of, 'temp.smf', Cscalp_w, Escalp);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cscalp_w,Escalp] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', Cskull_w, Eskull);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cskull_w,Eskull] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', CCSF_w, ECSF);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[CCSF_w,ECSF] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

Mesh_WriteSMF(of, 'temp.smf', Cbrain_w, Ebrain);
a = sprintf('"%s" -c "%sStepSc2.txt" "%s%s.smf"', conf.showmesh, of,of,'temp');
[status, result] = system(a);
if status ~=0; error('Mesh_generation:system', 'Failed to execute: %s', result); end
movefile([of 'ScS.smf'], [of  'temp.smf'])
[Cbrain_w,Ebrain] = mesh_readsmf([of 'temp.smf'],0,0,0,1); 

[so2, k1,k2] = mesh_check_intersection(Cbrain_w(:,2:4), CCSF_w, ECSF);
Cbrain_w(:,2:4) = so2;
[so2, k1,k2] = mesh_check_intersection(CCSF_w(:,2:4), Cskull_w, Eskull);
CCSF_w(:,2:4) = so2;
[so2, k1,k2] = mesh_check_intersection(Cskull_w(:,2:4), Cscalp_w, Escalp);
Cskull_w(:,2:4) = so2;


%%%%%%%


ind = warping_scalp_eloc_index(Pt, Cscalp_w,index_kdm);

% new electrode loc. and the chosen indices
Pt1 = Pt(ind,:);
[Ptm, dmi] = funrstp2([0 0 0 0 0 0], Pt1, Cscalp_w, Escalp);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Mesh_WriteSMF(of, name, Coord, Elem);
nnp = size(Coord,1); 
nel = size(Elem,1);
fid = fopen([of name], 'w');
fprintf(fid,'v %f %f %f \r\n',Coord(:,2:4)');
fprintf(fid,'t %d %d %d \r\n',Elem(:,2:4)');
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Coordw = warped_mesh(Coord,A,W,p);

r = Coord(:,2:4);
rw = warp_lm(r,A,W,p) + r;
Coordw = Coord;
Coordw(:,2:4) = rw;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [rw] = warp_lm(r,A,W,p)
rw = r * A(1:3,1:3) + repmat(A(4,:), size(r,1), 1);
for i = 1 : size(p,1)
    U = sqrt(sum((r - repmat(p(i,:), size(r,1),1)).^2, 2));  
    rw = rw + U * W(i,:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, A, e, LMm2, Ptbu, back] = find_warping(Coord, Elem, LMm, Fm, Fd, pos,index_kdm);
% Coord, Elem is the mesh that will be warped 
% LMm : Landmarks on the mesh (n2 x 3)
% Fm : fiducials on the mesh (n1 x 3)
% Fd : fiducials on the digitizer data (n1 x 3)
% pos : digitizer data (ne x 3)

ne = size(pos, 1); % number of electrodes
n1 = size(Fm, 1);  % number of fiducials
n2 = size(LMm, 1); % number of landmarks

options = optimset('MaxFunEvals', 100000, 'MaxIter', 100000, 'TolFun',1e-6);
Xo = [0 0 0 0 0 0];

% find translation and rotation for fiducials
X = fminsearch(@(X) funrstPP(X, Fd, Fm), Xo, options);
X = fminsearch(@(X) funrstPP(X, Fd, Fm), X, options);

% find translated digitizer fiducials
[d, Fdt] = warping_distafterwarping(X, Fd, Fm); 

% find translated and rotated digitizer locations
[d, Pt] = warping_distafterwarping(X, pos, ones(ne,3));

% Ptm are the rotated and translated digitizer locations 
% and moved to the closest point on the mesh 

Ptbu = Pt;
Pt = Pt(index_kdm,:);

[Ptm, dmi] = funrstp2(X, pos(index_kdm,:), Coord, Elem);
%[Ptm, dmi] = funrstp2(X, pos, Coord, Elem);

% find the index and distance between minimum distance Ptm and LMm  
for i = 1 : n2
    K = Ptm - ones(size(Ptm,1),1) * LMm(i,:);
    L = sum(K.*K,2);
    [k,l] = min(L);
    indm(i) = l;
    minm(i) = k;
end

% find the landmarks on the translated and rotated digitizer locations
LMd = Pt(indm,:); % on digitizer

LMm2 = Ptm(indm,:); % on model

% calculate the warping transformation
[W,A,e] = warp_transform(LMm2, LMd);

% warp_back parameters
[W2,A2,e2] = warp_transform(LMd, LMm2);
back.W = W2;
back.A = A2;
back.e = e2;
back.LMd = LMd;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [W, A, e] = warp_transform(p, q)
% [W,A,e]=warp_transform(p,q)
% calculates nonlinear transformatin coefficents see Ermer's Thesis
% p... = Landmarks in system 1
% q... = landmarks in system 2
% e =warp energy
K = zeros(size(p,1),size(p,1));
for i = 1:size(K,1)
    for j = 1:size(K,2)
        K(i,j) = norm(p(i,:)-p(j,:));    
    end
end
P = [p ones(size(p,1),1)];
L = [K P;P' zeros(4,4)];
D = [q-p;zeros(4,3)];
H = pinv(L)*D; % instead of L\D - June 13, 2013 Zeynep
W = H(1:size(p,1),:);
A = H(size(p,1)+1:end,:);
e = sum(diag(W'*K*W));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function d = funrstPP(X, F, Fe);
% F is the point set of the digitizer
% d is the distance between translated and rotated F and Fe
% X is the vector of translation and rotation parameters

tx = X(1);
ty = X(2);
tz = X(3);
alpx = X(4)*pi/180;
alpy = X(5)*pi/180;
alpz = X(6)*pi/180;

x = F(:,1);
y = F(:,2);
z = F(:,3);

% rotation around x-axis
x1 = x;
y1 = y * cos(alpx) - z * sin(alpx);
z1 = y * sin(alpx) + z * cos(alpx);

% rotation around y-axis
x2 = z1 * sin(alpy) + x1 * cos(alpy);
y2 = y1;
z2 = z1 * cos(alpy) - x1 * sin(alpy);

% rotation around z-axis
x3 = x2 * cos(alpz) - y2 * sin(alpz);
y3 = x2 * sin(alpz) + y2 * cos(alpz);
z3 = z2;

% translation
x4 = x3 + tx;
y4 = y3 + ty;
z4 = z3 + tz;

N = size(Fe,1);

Ma = Fe - [x4 y4 z4];
d = sum(sqrt(sum(Ma.*Ma,2)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [F2, dmi] = funrstp2(X, F, Coord, Elem);
% F : digitizer points
% F2 : are the rotated and translated version of F according to X
% they are moved to the closest points on the mesh.

tx = X(1);
ty = X(2);
tz = X(3);
alpx = X(4) * pi / 180;
alpy = X(5) * pi / 180;
alpz = X(6) * pi / 180;

x2 = F(:,1);
y2 = F(:,2);
z2 = F(:,3);

% rotation around x-axis
x3 = x2;
y3 = y2 * cos(alpx) - z2 * sin(alpx);
z3 = y2 * sin(alpx) + z2 * cos(alpx);

% rotation around y-axis
x4 = z3 * sin(alpy) + x3 * cos(alpy);
y4 = y3;
z4 = z3 * cos(alpy) - x3 * sin(alpy);

% rotation around z-axis
x5 = x4 * cos(alpz) - y4 * sin(alpz);
y5 = x4 * sin(alpz) + y4 * cos(alpz);
z5 = z4;

% translation
x1 = x5 + tx;
y1 = y5 + ty;
z1 = z5 + tz;

F2 = [x1 y1 z1];

hh = waitbar(0,'calculating the distance between digitizer and mesh');
for i = 1 : length(F2);
    waitbar(i/length(F2));
    [dm, Pm] = warping_distmeshpoint(F2(i,:), Coord, Elem);
    dmi(i) = dm;    Pmi(i,:) = Pm;
end; 
close(hh);
F2 = Pmi;
