% segm_inner_skull() - Performs inner skull segmentation
%
% Usage:
%   >> Sk_in = segm_inner_skull(b, Sk_out, X_dark, Bra);
%
% Inputs:
%   b - imput image (filtered MR image)
%   Sk_out - outer skull mask
%   X_dark - dark regions of b
%   Bra - Brain mask
%
% Outputs:
%   Sk_in - inner skull mask
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

function Sk_in = segm_inner_skull(b, Sk_out, X_dark, Bra, WMp);
% inner skull extraction 

% structuring elements
C1 = ones(3,3,3);
R1 = zeros(3,3,3); R1(2,2,:) =1; R1(1,2,2)=1; R1(2,1,2)=1; R1(3,2,2)=1; R1(2,3,2)=1;
O2 = ones(5,5,5); O2(1,5,:)=0; O2(1,1,:) = 0; O2(5,5,:)=0; O2(5,1,:)=0;
                  O2(1,:,5)=0; O2(1,:,1) = 0; O2(5,:,5)=0; O2(5,:,1)=0;
                  O2(:,1,5)=0; O2(:,1,1) = 0; O2(:,5,5)=0; O2(:,5,1)=0;
O4 = ones(10,10,10); O4(1,10,:)=0; O4(1,1,:) = 0; O4(10,10,:)=0; O4(10,1,:)=0;
                  O4(1,:,10)=0; O4(1,:,1) = 0; O4(10,:,10)=0; O4(10,:,1)=0;
                  O4(:,1,10)=0; O4(:,1,1) = 0; O4(:,10,10)=0; O4(:,10,1)=0;
           
                  
% dilate brain, discard the regions of outer skull 

Sk_oute = imerode3D(Sk_out, C1);
X_o = double(Sk_oute).*b;
clear Sk_oute
%X_bright = X_o >= t_skull;
X_bright = (X_o & not(X_dark))>0;
clear X_o

%B_dm = imdilate(Bra, C1);
B_dm = imdilate3D(Bra, O2); % O2 olarak degistirildi 06.06.08
X_u_in = X_bright | B_dm; clear B_dm;
X_open = imerode3D(imerode3D(X_u_in, O2), O2);
X_open = imdilate3D(imdilate3D(X_open, O2), O2);
Sk_out_e = imerode3D(Sk_out, O2);
Sk_out_e = imerode3D(Sk_out_e, O2);

Sk_out_e = imopen3D(Sk_out_e, O2); % new!
Sk_out_e = utilsegm_regiongrow(Sk_out_e, WMp(3), WMp(1), WMp(2),8); % new!

Sk_in = X_open | Sk_out_e;

% to delete the lower regions
% 06.05.08
se = strel('disk',5);
X1 = imopen3D(X_u_in,se);
R1 = utilsegm_regiongrow(X1, WMp(3), WMp(1), WMp(2), 4);
X2 = imdilate3D(imdilate3D(R1, O2), O2); 
Sk_in = X2 & Sk_in;

N = 10;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;

Sk_in = imopen3D(Sk_in, ON); % 2/9/2010
Sk_in = imclose3D(Sk_in, ON); % 6/11/2011


