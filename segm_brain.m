% Segm_brain() - Performs brain segmentation
%
% Usage:
%   >> [Bra] = Segm_brain(b,Sca, sli, WMp,sl, st);
%
% Inputs:
%   b - input image (filtered MR image)
%   sli - lowest point for crebellum
%   WMp - White matter point
%   sl,st - fill level and threshold for watershed segmentation
%
% Outputs:
%   Bra - brain mask% Segm_brain() - Performs brain segmentation
%
% Usage:
%   >> [Bra] = Segm_brain(b,Sca, sli, WMp,sl, st);
%
% Inputs:
%   b - input image (filtered MR image)
%   sli - lowest point for crebellum
%   WMp - White matter point
%   sl,st - fill level and threshold for watershed segmentation
%
% Outputs:
%   Bra - brain mask
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


function [Bra] = segm_brain(b,Sca, sli, WMp,sl, st);
% sl: set level for watershed segm.
% st: set threshold for watershed segm.
%save brain_params b Sca sli WMp sl st
N = 3;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;
                  
% mask the image with head
b2 = double(Sca).*b;
%b2=double(Sca).*imc;
[K,L,M] = size(b);
% mask the lower part of the head non-brain regions
b3 = zeros(K,L,M);
b3(:,sli:L,:)=b2(:,sli:L,:);

mb = max(max(max(b3)));
c = matitk('SWS', [sl st], mb-double(b3),[],WMp);

% initial brain
bint = c(WMp(1),WMp(2),WMp(3));
bi = c==bint;

clear c b2 

[f2 n2 v2 z2] = matitk('FOMT',[4 100],single(b3));
%f2 = b3 > 60; f2  = not(f2);
clear n2 v2 z2
bb = bi & not(f2);
clear f2


se = strel('disk',5);
Bra = imclose3D(bb, se);
clear bb


A = Bra.*b;
[k1max, k2max, h] = utilsegm_thresh(A, 2);
at = max(max(max(b)));
thre = k1max + at*0.01;

se3 = strel('disk',10);
Brd = imerode3D(Bra,se3);

Bra2 = A > thre;
Bra2 = imopen3D(Bra2, ON);
R1 = utilsegm_regiongrow(Bra2, WMp(3), WMp(1)-2:WMp(1)+10, WMp(2)-2:WMp(2)+10, 4);
Bra3 = R1 | Brd;

% to get rid of the bone marrow 09/16/2008
d1 = (Bra2.*b);  
clear Bra2;
d1 = (d1 > (k1max + at*0.05));
d2 = utilsegm_regiongrow(d1, WMp(3), WMp(1):WMp(1)+10, WMp(2):WMp(2)+10, 8);
d3 = d1 & not(d2);

% check WMp(2)
dx = zeros(K,L,M);
dx(:,WMp(2):L,:) = d3(:,WMp(2):L,:);
clear d1 d2 d3
d4 = Bra3 & not(dx);
d4 = imopen3D(d4, ON);
d4 = utilsegm_regiongrow(d4, WMp(3), WMp(1):WMp(1)+10, WMp(2):WMp(2)+10, 4);
Bra = d4|Brd;

N = 10;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;

Bra = imopen3D(Bra, ON); % 2/9/2010
Bra = imclose3D(Bra, ON); % 6/11/2011

