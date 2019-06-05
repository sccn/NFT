% warping_scalp_eloc_index() - Finds the indices of the
%           electrodes outside the head region (electrodes
%           corresponding the jaw area etc.).
%
% Usage:
%   >> ind = warping_scalp_eloc_index(pos, Coord, Elem);
%
% Inputs:
%   Coord - coordinate matrix of input mesh
%   pos   - positions of the electrodes
%
% Outputs:
%   ind - indices
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

function ind = warping_scalp_eloc_index(pos, Coord, index_kdm);

% returns the index of the electrodes in "pos" that are outside of the
% mesh (Coord,Elem)
ind = [];
d3 = mean(Coord(:,2));
sd = std(Coord(:,2));
sd = 5;
xd1 = d3-sd;
xd2 = d3+sd;

% find the line that defines the bottom of the head
d = find(Coord(:,2) < xd1);
[k,l] = min(Coord(d,4));
x1 = Coord(d(l),2);
z1 = Coord(d(l),4);

d = find(Coord(:,2) < xd2);
[k,l] = min(Coord(d,4));
x2 = Coord(d(l),2);
z2 = Coord(d(l),4);

% aX+b=Z find a and b
if (abs(x1-x2)) < eps
    a = 0;
    b = z1;
else
    a = (z1-z2)/(x1-x2);
    b = z1-a*x1;
end


z_b = max(Coord(:,4))-min(Coord(:,4)); 
% aX+b-z < 0 
n=0;
for i=1:length(pos)
    w = a * pos(i,1) + b - pos(i,3);
    if w < 1
    %if w < -z_b * 0.1 
        n = n + 1;
        ind(n) = i;
    end
end
ind = intersect(ind,index_kdm);

