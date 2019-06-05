% segm_scalp() - Perform scalp segmentation
%
% Usage:
%   >> [Sca] = segm_scalp(b);
%
% Inputs:
%   b - input image (filtered MR image)
%
% Outputs:
%   Sca - scalp mask
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

function [Sca] = segm_scalp(b);

thrsc = 50; % consider using k1max in thresh2.m  rl=5
[K,L,M] = size(b);
[d2 g2] = matitk('FOMT',[2 thrsc],b); 
%d2 = b < 95;

d2(:,:,1) = ones(K,L); d2(:,:,M) = ones(K,L); 
d2(:,1,:) = ones(K,M); d2(:,L,:) = ones(K,M);
d2(1,:,:) = ones(L,M); d2(K,:,:) = ones(L,M);

N = 4;
ON = ones(N,N,N); ON(1,N,:)=0; ON(1,1,:) = 0; ON(N,N,:)=0; ON(N,1,:)=0;
                  ON(1,:,N)=0; ON(1,:,1) = 0; ON(N,:,N)=0; ON(N,:,1)=0;
                  ON(:,1,N)=0; ON(:,1,1) = 0; ON(:,N,N)=0; ON(:,N,1)=0;
c2 = imclose3D(not(d2),ON);

q1 = round(K/2); q2 = round(L/2); q3 = round(M/2);
[w1,w2] = find(c2(q1-10:q1+10,q2-10:q2+10,q3)==1);
R1 = c2;%utilsegm_regiongrow(c2, q3, w1+q1, w2+q2, 8);
R2 = R1;
for i = 1:L
    ek = (reshape(R2(:,i,:),K,M));
    ek = imfill(ek,'holes');
    R2(:,i,:) = ek;
end
R3 = R2;
clear R2;
for i = 1:M
    ek = R3(:,:,i);
    ek = imfill(ek,'holes');
    R3(:,:,i) = ek;
end
R4 = R3;
clear R3;
for i = 1:K
    ek = (reshape(R4(i,:,:),L,M));
    ek = imfill(ek,'holes');
    R4(i,:,:) = ek;
end

% convert sagittal to axial
% clean the small pieces using bwselect around the scalp
for i=1:L
    a = reshape(R4(:,i,:),K,M);
    Le = bweuler(a,4);
    ax(:,:,i) = a;
    if Le>1
        a = bwselect(a, round(M/2), round(K/2), 4);
        ax(:,:,i) = a;
    end
end

% convert axial to sagittal
for i=1:M
    Sca(:,:,i) = reshape(ax(:,i,:),K,L);
end
Sca = R4;
Sca = logical(Sca);
Sca(:,:,1) = zeros(K,L); Sca(:,:,M) = zeros(K,L); 
Sca(:,1,:) = zeros(K,M); Sca(:,L,:) = zeros(K,M);
Sca(1,:,:) = zeros(L,M); Sca(K,:,:) = zeros(L,M);
 
N = 5;
ON = ones(N,N,N); ON(1,N,:) = 0; ON(1,1,:) = 0; ON(N,N,:) = 0; ON(N,1,:) = 0;
                  ON(1,:,N) = 0; ON(1,:,1) = 0; ON(N,:,N) = 0; ON(N,:,1) = 0;
                  ON(:,1,N) = 0; ON(:,1,1) = 0; ON(:,N,N) = 0; ON(:,N,1) = 0;

Sca = imclose3D(Sca, ON); % 6/11/2011
Sca = imopen3D(Sca, ON); 

Sca(:,:,1) = zeros(K,L); Sca(:,:,M) = zeros(K,L); 
Sca(:,1,:) = zeros(K,M); Sca(:,L,:) = zeros(K,M);
Sca(1,:,:) = zeros(L,M); Sca(K,:,:) = zeros(L,M);


%%%%%%%%%%%%%%
