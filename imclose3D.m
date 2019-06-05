% imclose3D() - Applies the closing operation to the image in 3D
%
% Usage:
%   >> B = imclose3D(A, se);
%
% Inputs:
%   A - input image
%   se - structuring element
%
% Outputs:
%   B - output image
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

function B = imclose3D(A, se);

sz = size(A);

% axial 
A1 = imclose(A, se);
A1 = A1 - min(min(min(A1)));

% sagittal
for i = 1:sz(2); R1(:,:,i) = (reshape(A(:,i,:),sz(1),sz(3))); end
R1 = imclose(R1, se);
for i=1:sz(3);   R2(:,:,i) = reshape(R1(:,i,:),sz(1),sz(2)); end
R2 = R2 - min(min(min(R2)));

% coronal
for i = 1:sz(1); R3(:,:,i) = (reshape(A(i,:,:),sz(2),sz(3))); end
R3 = imclose(R3, se);
for i=1:sz(3);   R4(:,:,i) = squeeze(R3(:,i,:))'; end
R4 = R4 - min(min(min(R4)));

B = A1 | R2 | R4;
