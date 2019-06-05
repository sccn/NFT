% warping_eloc_to_MNI() - Transforms a set of points
%                 from digitized electrode space to MNI space 
%                 using the warping parameters
%
%
% Usage:
%   >> Pnew = warping_eloc_to_MNI(warping_param, P);
%


% Inputs:
%   P  - point set in MNI space
%   warping_param - Parameters generated when warping MNI to eloc
%
% Outputs:
%   Pnew - point set in digitized electrode space
%
%
% Author: Zeynep Akalin Acar, SCCN, 2013

% Copyright (C) 2013 Zeynep Akalin Acar, SCCN, zeynep@sccn.ucsd.edu
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

function Pnew = warping_eloc_to_MNI(warping_param, P)

sc = warping_param.initial.sc;
T  = warping_param.initial.T;
R  = warping_param.initial.R;

A = warping_param.back.A;
W = warping_param.back.W;
L = warping_param.back.LMd;

Pnew = warping_apply_initial_registration(sc, T, R, P);
Pnew = Pnew + warp_lm(Pnew, A, W, L);
%Pnew = P;

