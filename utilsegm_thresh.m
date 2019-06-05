% utilsegm_thresh() - Finds thresholds
%
% Usage:
%   >> [k1max, k2max, h] = utilsegm_thresh(A,n);
%
% Inputs:
%   A - imput image
%   n - number of thresholds (1 or 2)
%
% Outputs:
%   k1max, k2max - thresholds
%   h - histogram
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

function [k1max, k2max, h] = utilsegm_thresh(A,n);
% if n==1 there is one threshold
% if n==2 there are two thresholds
% A is the image to be segmented
% thrk1 and thrk2 are the threshold values

L=round(max(max(max(A)))-min(min(min(A)))+1);
h=hist(A(:),L);
%h=histog(A,L);
h=h(2:L); L=L-1; % omit the pixels that are 0


N=sum(h); p=h/N; ip=[1:L].*p; muT=sum(ip);

Bmax=-1;
k1max=0;
k2max=0;

if n==1
   sigB=zeros(1,L-1);
   for k=1:L-1
   	w0=sum(p(1:k));
	 	w1=sum(p(k+1:L));
  		mu0=sum(ip(1:k))/w0;
	   mu1=sum(ip(k+1:L))/w1;
   	sigB(k)=w0*w1*(mu1-mu0)^2;
	end
	[Bmax k1max]=max(sigB);
end

if n==2
   m0=0;
   w0=0;
   w2base=sum(p(1:L));
   m2base=sum(ip(1:L));
   
   for k1=1:L-2
      w0=w0+p(k1);
      m0=m0+ip(k1);
      mu0=m0/w0;
      
      w2base=w2base-p(k1);
      m2base=m2base-ip(k1);
      
      w1=0;
      m1=0;
      w2=w2base;
      m2=m2base;
      for k2=k1+1:L-1
         w1=w1+p(k2);
         w2=w2-p(k2);
         m1=m1+ip(k2);
         m2=m2-ip(k2);
         
         if w1 ~= 0 & w2 ~= 0
            mu1=m1/w1;
            mu2=m2/w2;
            sigB=w0*(mu0-muT)^2+w1*(mu1-muT)^2+w2*(mu2-muT)^2;
            if sigB > Bmax
               k1max=k1;
               k2max=k2;
               Bmax=sigB;
            end
         end
      end
   end
end

