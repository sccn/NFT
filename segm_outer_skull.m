% Segm_Outer_skull() - Performs outer skull segmentation
%
% Usage:
%   >> [Sk_out, X_dark] = Segm_Outer_skull(b, Sca, Bra, sli_eyes);
%
% Inputs:
%   b - imput image (filtered MR image)
%   Sca - Scalp mask
%   Bra - Brain mask
%   sli_eyes - slice of the eyes
%
% Outputs:
%   Sk_out - outer skull mask
%   X_dark - dark regions of b
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

function [Sk_out, X_dark,thr] = Segm_Outer_skull(b, Sca, Bra, sli_eyes);

% outer skull extraction 
 %save segmsk b Sca Bra sli_eyes
% structuring elements
C1 = ones(3,3,3);
R1 = zeros(3,3,3); R1(2,2,:) = 1; R1(1,2,2) = 1; R1(2,1,2) = 1; R1(3,2,2) = 1; R1(2,3,2) = 1;
O2 = ones(5,5,5); O2(1,5,:) = 0; O2(1,1,:) = 0; O2(5,5,:) = 0; O2(5,1,:) = 0;
                  O2(1,:,5) = 0; O2(1,:,1) = 0; O2(5,:,5) = 0; O2(5,:,1) = 0;
                  O2(:,1,5) = 0; O2(:,1,1) = 0; O2(:,5,5) = 0; O2(:,5,1) = 0;
N = 8;
ON = ones(N,N,N); ON(1,N,:) = 0; ON(1,1,:) = 0; ON(N,N,:) = 0; ON(N,1,:) = 0;
                  ON(1,:,N) = 0; ON(1,:,1) = 0; ON(N,:,N) = 0; ON(N,:,1) = 0;
                  ON(:,1,N) = 0; ON(:,1,1) = 0; ON(:,N,N) = 0; ON(:,N,1) = 0;

N = 10;
ON10 = ones(N,N,N); ON(1,N,:) = 0; ON(1,1,:) = 0; ON(N,N,:) = 0; ON(N,1,:) = 0;
                  ON(1,:,N) = 0; ON(1,:,1) = 0; ON(N,:,N) = 0; ON(N,:,1) = 0;
                  ON(:,1,N) = 0; ON(:,1,1) = 0; ON(:,N,N) = 0; ON(:,N,1) = 0;

% brain and scalp masks are Bra and Sca

[K, L, M] = size(b);

[k1max, k2max, h] = utilsegm_thresh(b, 2);
thr = (k1max-max(max(max(b)))*0.01)
%thr = 70 % child
%thr = 40
X_dark = b < thr;

X_dark = logical(X_dark);

% select eyes
h = figure; imagesc(reshape(X_dark(:,sli_eyes,:),K,M)); colormap gray;
[xp,yp] = ginput(2); xp = round(xp); yp=round(yp);
close(h); pause(1);               

Se1 = imdilate3D(imerode3D(Sca,ones(25,25,25)),ones(25,25,25));
Se2 = imerode3D(Se1, ones(7,7,7));
%clear Se1
  
%B_dm = imdilate3D(imdilate3D(Bra, C1),C1);           
B_dm = imdilate3D(Bra, ON);           
X_u_out = X_dark | B_dm;                         
X_i_out = X_u_out & Se2;

% mask the lower part of the skull
sli = 50;
X_i_out(:,1:sli,:) = 0;

%load X_i_out
%imagesc(squeeze(X_i_out(:,80,:))); [x,y] = ginput(1); x=53, y=148
%R1 = utilsegm_regiongrow(X_i_out,56,80,133,4);
%R2 = utilsegm_regiongrow(X_i_out,53,80,148,4);
%X_i_out = X_i_out & (not(R1)); 
%X_i_out = X_i_out & (not(R2)); 
%clear Bra b

A = X_dark & Sca; % to delete the connection between scalp and eyes
G1 = utilsegm_regiongrow(A,xp(1),sli_eyes, yp(1),8);
G2 = utilsegm_regiongrow(A,xp(2),sli_eyes, yp(2),8);

% discard eyes
X_lc = X_i_out;
X_lc = X_lc & not(G1);
X_lc = X_lc & not(G2);

Xc = imdilate3D(imdilate3D(X_lc, O2), O2);
Xc = imerode3D(imerode3D(Xc, O2), O2);
Xc = imfill(Xc,'holes');  % Feb 23 2015
Sk_out = Xc & Se2;

Sk_out = logical(Sk_out);
Sk_out = imopen3D(Sk_out, ON10); 
Sk_out = imclose3D(Sk_out, ON10); % 6/11/2011
