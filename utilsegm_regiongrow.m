% utilsegm_regiongrow() - Performs region growing
%
% Usage:
%   >> R = utilsegm_regiongrow(A,msli,xi,yi,n);
%
% Inputs:
%   A - input image (binary)
%   msli - slice number for seed
%   xi, yi - x,y locations for the seed
%   n - connectivity (4 or 8)
%
% Outputs:
%   R - output image
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


function R = utilsegm_regiongrow(A,msli,xi,yi,n);
% 3D region growing by bwselect (non zero terms in A)
% returns the object in R
% n =4 or 8 determines connectivity

% enlarge the seed region if it is a single point % 10-26-2010
if (length(xi) < 1) & (xi > 10) & (xi < sz(1)-10)
    xi = xi - 5 : xi + 5;
end
if (length(yi) < 1) & (yi > 10) & (yi < sz(2)-10)
    yi = yi - 5 : yi + 5;
end

sz=size(A);
if length(sz) ~= 3
   error('input must be a 3D volume');
end
if xi<1 | xi>sz(2)
   error('x value out of range');
end
if yi<1 | yi>sz(1)
   error('y value out of range');
end
if A(yi,xi,msli)==0
    error('seed is empty');
end


csli = bwselect(A(:,:,msli), xi, yi, n);

R = uint8(A);
R(:) = 0;

R(:,:,msli) = csli;
csli = bwmorph(csli,'open');
for sli = msli-1:-1:1
   map = A(:,:,sli)&csli;
   ind = find(map==1);
   y = mod(ind,sz(1));
   x = (ind-y)/sz(1)+1;
   if ~isempty(x)&~isempty(y)
       csli=bwselect(A(:,:,sli),x,y,n);      
   else
       csli=zeros(sz(1),sz(2));
   end
   R(:,:,sli)=csli;
   csli=bwmorph(csli,'open');
end

csli=bwselect(A(:,:,msli),xi,yi,n);
csli=bwmorph(csli,'open');
for sli=msli+1:sz(3)
   map=A(:,:,sli)&csli;
   ind=find(map==1);
   y=mod(ind,sz(1));
   x=(ind-y)/sz(1)+1;
   if ~isempty(x)&~isempty(y)
       csli=bwselect(A(:,:,sli),x,y,n);      
   else
       csli=zeros(sz(1),sz(2));
   end
   R(:,:,sli)=csli;
   csli=bwmorph(csli,'open');
end


