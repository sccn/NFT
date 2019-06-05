% segm_final_skull() - Corrects scalp and outer skull masks with respect to
%                    inner skull mask
%
% Usage:
%   >> [Sca, Sk_out] = segm_final_skull(Sca, Sk_out, Sk_in, WMp);
%
% Inputs:
%   Sca    - scalp mask
%   Sk_out - outer skull mask
%   Sk_in  - inner skull mask
%   WMp    - white matter point
%
% Outputs:
%   Sca    - scalp mask
%   Sk_out - outer skull mask
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

function [Sca, Sk_out] = segm_final_skull(Sca, Sk_out, Sk_in, WMp);

O2 = ones(5,5,5); O2(1,5,:)=0; O2(1,1,:) = 0; O2(5,5,:)=0; O2(5,1,:)=0;
                  O2(1,:,5)=0; O2(1,:,1) = 0; O2(5,:,5)=0; O2(5,:,1)=0;
                  O2(:,1,5)=0; O2(:,1,1) = 0; O2(:,5,5)=0; O2(:,5,1)=0;
N = 8;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;
                  
[K,L,M] = size(Sca);
B_dm = imclose3D(Sk_in, O2);
B_dm = imdilate3D(B_dm, O2); 
Sk_out = Sk_out | B_dm;

se = strel('disk',5);
Sk_out = imopen3D(Sk_out,se);
%Sk_out = utilsegm_regiongrow(Sk_out, WMp(3), WMp(1), WMp(2), 4);
Sk_out = imclose3D(Sk_out,ON);

se1 = strel('disk',3);
for i = 1:L
    im1 = squeeze(Sk_out(:,i,:));
    im2 = imopen(im1, se1);
    Sk_out(:,i,:) = im2;
end
%Sk_out = utilsegm_regiongrow(Sk_out, WMp(3), WMp(1), WMp(2), 4);
Sk_out = imclose3D(Sk_out, ON);

B_dm = imdilate3D(B_dm, O2);
Sca = Sca | B_dm;


N = 10;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;

Sk_out = imopen3D(Sk_out, ON); 
Sk_out = imclose3D(Sk_out, ON); % 6/11/2011
