% utilbem_headplot() - plot a spherically-splined EEG field map on a semi-realistic 
%              3-D head model. Can 3-D rotate the head image using the left 
%              mouse button.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (C) Arnaud Delorme, Colin Humphries and Scott Makeig, 
%               CNL / Salk Institute, Feb. 1998
% Copyright (C) Zeynep Akalin Acar, SCCN, March 2008 
%
% Spherical spline method: Perrin et al. (1989) Electroenceph clin Neurophys
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

% $Log: headplot.m,v $
% Based on Revision 1.76  of headplot.m from EEGLAB

function P = utilbem_headplot(pot, mesh, eloc)

ma = mean(mesh.coord);
%ma = mean(eloc);
mesh.coord = mesh.coord-ones(size(mesh.coord,1),1) * ma;
eloc = eloc-ones(size(eloc,1),1) * ma;
[HeadCenter,radius,residuals] = spherefit(eloc(:,1),eloc(:,2),eloc(:,3));
HeadCenter = HeadCenter';

%Xeori = eloc(:,1);
%Yeori = eloc(:,2);
%Zeori = eloc(:,3);
   
newcoordsnorm      = eloc - ones(size(eloc,1),1)*HeadCenter;
tmpnorm            = sqrt(sum(newcoordsnorm.^2,2));
%Xe = Xeori./tmpnorm;
%Ye = Yeori./tmpnorm;
%Ze = Zeori./tmpnorm;
Xe = newcoordsnorm(:,1)./tmpnorm;
Ye = newcoordsnorm(:,2)./tmpnorm;
Ze = newcoordsnorm(:,3)./tmpnorm;

eloc2 = eloc; % to generate the spline matrix
mine = min(min(eloc2)); % make it positive
eloc2 = eloc - mine;
Xeori2 = eloc2(:,1);
Yeori2 = eloc2(:,2);
Zeori2 = eloc2(:,3);
newcoordsnorm2      = eloc2 - ones(size(eloc2,1),1)*HeadCenter;
tmpnorm2            = sqrt(sum(newcoordsnorm2.^2,2));
Xe2 = Xeori2./tmpnorm2;
Ye2 = Yeori2./tmpnorm2;
Ze2 = Zeori2./tmpnorm2;

fprintf('Setting up splining matrix.\n');    
enum = length(Xe2);
onemat = ones(enum,1);
G = zeros(enum,enum);
for i = 1:enum
    ei = onemat-sqrt((Xe2(i)*onemat-Xe2).^2 + (Ye2(i)*onemat-Ye2).^2 + ...
                     (Ze2(i)*onemat-Ze2).^2); % default was /2 and no sqrt
    gx = zeros(1,enum);
    for j = 1:enum
        gx(j) = calcgx(ei(j));
    end
    G(i,:) = gx;
end
  
fprintf('Calculating splining matrix...\n')

TRI1 = mesh.elem(1 : mesh.bnd(1,1), :);
POS  = mesh.coord(1 : max(max(TRI1)), :);
center = mean(POS);

if exist('index1') ~= 1, index1 = sort(unique(TRI1(:))); end;
if exist('TRI2')   ~= 1, TRI2 = []; end;
if exist('NORM')   ~= 1, NORM = []; end;
    
newPOS = POS(index1,:);
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Project head vertices onto unit sphere
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
spherePOS      = newPOS-ones(size(newPOS,1),1)*HeadCenter; % recenter
nPOSnorm       = sqrt(sum(spherePOS.^2,2));
spherePOS(:,1) = spherePOS(:,1)./nPOSnorm;
spherePOS(:,2) = spherePOS(:,2)./nPOSnorm;
spherePOS(:,3) = spherePOS(:,3)./nPOSnorm;
x = spherePOS(:,1);
y = spherePOS(:,2);
z = spherePOS(:,3);
    
 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate g(x) for sphere mesh vertices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

gx = fastcalcgx(x,y,z,Xe,Ye,Ze);

%%%%%%%%%%%%%%%%%%%%%%%%%% Make the plot %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   



% Perform interpolation
values=pot;
meanval = mean(values); values = values - meanval; % make mean zero
onemat = ones(enum,1);
lamd = 0.1;
C = pinv([(G + lamd);ones(1,enum)]) * [values(:);0]; % fixing division error
P = zeros(1,size(gx,1));
for j = 1:size(gx,1)
   P(j) = dot(C,gx(j,:));
end
P = P + meanval;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot surfaces
%%%%%%%%%%%%%%%%%%%%%%%%%%
cla % clear axis
HeadAxes = gca;
Coord = zeros(length(POS),4);
Elem = zeros(length(TRI1),4);
Coord(:,2:4) = POS + ones(size(POS,1),1) * ma;
Elem(:,2:4) = TRI1;
Coord(:,1) = [1:length(POS)]';
Elem(:,1) = [1:length(TRI1)]';
utilbem_PlotPot(Coord, Elem, P)
 
if ~isempty(TRI2)
    FCmap = [g.colormap; g.colormap(end,:); FaceColor; FaceColor; FaceColor];
    colormap(FCmap)
    W = ones(1,size(POS,1))*(m+2);
    p2 = patch('Vertices',POS,'Faces',TRI2,'FaceColor','interp',...
               'FaceVertexCdata',W(:)); %%%%%%%% Plot face and lower head %%%%%%
else 
    p2 = [];
end;

eloc = eloc + ones(size(eloc,1),1) * ma;
plotelec(eloc);

rotate3d on;   % Allow 3-D rotation of the plot by dragging the

function [out] = calcgx(in)

out = 0;
m = 4;       % 4th degree Legendre polynomial
for n = 1:7  % compute 7 terms
  L = legendre(n,in);
    out = out + ((2*n+1)/(n^m*(n+1)^m))*L(1);
end
out = out/(4*pi);

    

function gx = fastcalcgx(x,y,z,Xe,Ye,Ze)
onemat = ones(length(x),length(Xe));
EI = onemat - sqrt((repmat(x,1,length(Xe)) - repmat(Xe',length(x),1)).^2 +... 
                    (repmat(y,1,length(Xe)) - repmat(Ye',length(x),1)).^2 +...
                    (repmat(z,1,length(Xe)) - repmat(Ze',length(x),1)).^2);
gx = zeros(length(x),length(Xe));
m = 4;
hwbend = 7;
for n = 1:7
    L = legendre(n,EI);
    gx = gx + ((2*n+1)/(n^m*(n+1)^m))*squeeze(L(1,:,:));
end
gx = gx/(4*pi);    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  distance() - function used in 'setup'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [out] = distance(w,p)
% w is a matrix of row vectors
% p is a matrix of column vectors

l1 = size(w,1);
l2 = size(p,2);
out = zeros(l1,l2);

for i = 1:l1
  x = w(i,:)'*ones(1,l2);
  out(i,:) = sum((x-p).^2).^.5;
end

% %%%%%%%%%%%%%%%
% plot electrodes
% %%%%%%%%%%%%%%%
function plotelec(newElect);
    
for i = 1:size(newElect,1)
    if newElect(i,:) ~= [0 0 0]  % plot radial lines to electrode sites
    % plot electrode markers
        line(newElect(:,1),newElect(:,2),newElect(:,3),'marker',...
            '.','markersize',20,'color','r','linestyle','none');
    end
end;

function [center,radius,residuals] = spherefit(x,y,z)
%SPHEREFIT find least squares sphere
%
% Fit a sphere to a set of xyz data points
% [center,radius,residuals] = shperefit(X)
% [center,radius,residuals] = spherefit(x,y,z);
% Input
% x,y,z Cartesian data, n x 3 matrix or three vectors (n x 1 or 1 x n)
% Output
% center least squares sphere center coordinates, == [xc yc zc]
% radius radius of curvature
% residuals residuals in the radial direction
%
% Fit the equation of a sphere in Cartesian coordinates to a set of xyz
% data points by solving the overdetermined system of normal equations,
% ie, x^2 + y^2 + z^2 + a*x + b*y + c*z + d = 0
% The least squares sphere has radius R = sqrt((a^2+b^2+c^2)/4-d) and
% center coordinates (x,y,z) = (-a/2,-b/2,-c/2)
%
% http://www.mathworks.de/matlabcentral/newsreader/view_thread/157211

error(nargchk(1,3,nargin)); % check input arguments
if nargin == 1 % n x 3 matrix
   if size(x,2) ~= 3
      error ('input data must have three columns')
   else
      z = x(:,3); % save columns as x,y,z vectors
      y = x(:,2);
      x = x(:,1);
   end
elseif nargin == 3 % three x,y,z vectors
   x = x(:); % force into columns
   y = y(:);
   z = z(:);
   if ~isequal(length(x),length(y),length(z)) % same length ?
      error('input vectors must be same length');
   end
else % must have one or three inputs
   error('invalid input, n x 3 matrix or 3 n x 1 vectors expected');
end

% need four or more data points
if length(x) < 4
   error('must have at least four points to fit a unique sphere');
end

% solve linear system of normal equations
A = [x y z ones(size(x))];
b = -(x.^2 + y.^2 + z.^2);
a = A \ b;

% return center coordinates and sphere radius
center = -a(1:3)./2;
radius = sqrt(sum(center.^2)-a(4));

% calculate residuals
if nargout > 2
   residuals = radius - sqrt(sum(bsxfun(@minus,[x y z],center.').^2,2));
end


