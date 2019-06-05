% utilmesh_modify() - Removes or modifies bad elements in a given mesh
%
% Usage:
%   >> [Coord, Elem, Check]=utilmesh_modify(Coord,Elem,Check);
%
% Inputs:
%   Coord - coordinate matrix of input mesh
%   Elem  - connectivity matrix of input mesh
%   Check - (optional)
%
% Outputs:
%   Coord - coordinate matrix of output mesh
%   Elem  - connectivity matrix of output mesh
%   Check - elements that are deleted
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

function [Coord, Elem, Check]=utilmesh_modify(Coord,Elem,Check);
%function [Elem, Coord]=CorrectMesh(Elem,Coord);

if exist('Check') ~= 1
    Check = [];
end

%threshold for bad elements
atresh=0.3;

AR=ElemAspect(Elem, Coord);

done = 0;

while done == 0
    [ARy,ARi]=sort(AR);
    ne=length(AR);
    done = 1;
    for t=1:ne
        if ARy(t) > atresh
            break;
        end
        k=ARi(t);
      
        % find shortest edge
        d1=EdgeLength(Elem(k,2), Elem(k,3), Coord);
        d2=EdgeLength(Elem(k,3), Elem(k,4), Coord);
        d3=EdgeLength(Elem(k,4), Elem(k,2), Coord);
         
        if d1 < d2 && d1 < d3
            n1=Elem(k,2); n2=Elem(k,3);
        elseif d2 < d1 && d2 < d3
            n1=Elem(k,3); n2=Elem(k,4);
        else
            n1=Elem(k,4); n2=Elem(k,2);
        end
        N1=unique([find(Elem(:,2)==n1)' find(Elem(:,3)==n1)' find(Elem(:,4)==n1)']);
        N2=unique([find(Elem(:,2)==n2)' find(Elem(:,3)==n2)' find(Elem(:,4)==n2)']);
     
     	% elements to be deleted
        D=intersect(N1,N2);
        check=intersect(Check,D);
      
     	% elements to be modified
        M=setdiff(union(N1,N2),D);
      
        % aspect ratios already calculated
        Ma0=AR(M);
        % calculate normals
        Mn0=ElemNormal(Elem(M,:), Coord);
        c1=Coord(n1,2:4);
        c2=Coord(n2,2:4);
        % set both points to the mid point
        c0=(c1+c2)/2;
        Coord(n1,2:4)=c0;
        Coord(n2,2:4)=c0;
        % calculate with new coordinate
        Ma1=ElemAspect(Elem(M,:),Coord);
        Mn1=ElemNormal(Elem(M,:),Coord);
     	     
        Mn=(Mn0.*Mn1) * [1 1 1]';
     
        overlap = 0;
        ovElem=[];
        ovNode=[];
        ovSel=[];
        lenM=length(M);
        for mi=1:lenM-1
            for mj=mi+1:lenM
                ei=Elem(M(mi),2:4);
                ej=Elem(M(mj),2:4);
                fi=find(ei == n2);
                ei(fi)=n1;
                fi=find(ej == n2);
                ej(fi)=n1;
                ei=sort(ei);
                ej=sort(ej);
                if ei == ej
                    on1=unique([find(Elem(:,2)==ei(1))' find(Elem(:,3)==ei(1))' find(Elem(:,4)==ei(1))']);
                    on2=unique([find(Elem(:,2)==ei(2))' find(Elem(:,3)==ei(2))' find(Elem(:,4)==ei(2))']);
                    on3=unique([find(Elem(:,2)==ei(3))' find(Elem(:,3)==ei(3))' find(Elem(:,4)==ei(3))']);
                    on1=setdiff(on1,D);
                    on2=setdiff(on2,D);
                    on3=setdiff(on3,D);
                    
                    on=[];
                    if (ei(1) ~= n1 && length(on1) == 2)
                        on=union(on,ei(1));
                    end
                    if (ei(2) ~= n1 && length(on2) == 2)
                        on=union(on,ei(2));
                    end
                    if (ei(3) ~= n1 && length(on3) == 2)
                        on=union(on,ei(3));
                    end
                    if isempty(on)
                        overlap=1;
                        M;
                        D;
                    else
                        ovNode=union(ovNode, on);
                        ovSel=union(ovSel, [mi mj]);
                    end
                end
            end
        end
        
        Mn(ovSel)=1;
        Ma1(ovSel)=1;
        ovElem=M(ovSel);
        if overlap == 0 && min(Mn) > 0 && norm(Ma1) > norm(Ma0)
            done=0;
            AR(:,M)=Ma1;    
            Elem=ReplaceNode(Elem,n2,n1);
            D = union(D, ovElem);
            Elem=DeleteRow(Elem, D);
            AR=DeleteRow(AR', D)';
            ovNode=union(ovNode, n2);
            ovNode=sort(ovNode, 'descend');
            for on=ovNode
                [Elem, Coord]=RemoveNode(Elem, Coord, on);
            end
            
            [ARy,ARi]=sort(AR);
            ne=length(AR);
        else
            % backoff
            Coord(n1,2:4)=c1;
            Coord(n2,2:4)=c2;
        end
    end
end

Coord(:,1)=[1:length(Coord(:,1))]';
Elem(:,1)=[1:length(Elem(:,1))]';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Aspect=ElemAspect(Elem,Coord);
% Elem'deki elemanlarin en kucuk edgelength'inin en buyugune orani.

ne=length(Elem(:,1));

nr=length(Coord(:,1));

for k=1:ne
   d(1)=EdgeLength(Elem(k,2), Elem(k,3), Coord);
   d(2)=EdgeLength(Elem(k,3), Elem(k,4), Coord);
   d(3)=EdgeLength(Elem(k,4), Elem(k,2), Coord);
   a=sort(d);
   Aspect(k)=a(1)/a(3);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Len=EdgeLength(n1,n2,Coord);

p1=Coord(n1,2:4);
p2=Coord(n2,2:4)-p1;
Len=sqrt(p2*p2');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Norm=ElemNormal(Elem,Coord);

ne=length(Elem(:,1));
nr=length(Coord(:,1));

for k=1:ne
   v1=Coord(Elem(k,3),2:4)-Coord(Elem(k,2),2:4);
   v2=Coord(Elem(k,4),2:4)-Coord(Elem(k,3),2:4);
   n=cross(v1, v2);
   Norm(k,1:3)=n/norm(n);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Elem=ReplaceNode(Elem,n2,n1);
i=find(Elem(:,2)==n2);
Elem(i,2)=n1*ones(length(i),1);

i=find(Elem(:,3)==n2);
Elem(i,3)=n1*ones(length(i),1);

i=find(Elem(:,4)==n2);
Elem(i,4)=n1*ones(length(i),1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Mat=DeleteRow(Mat, R);
nc=length(Mat(:,1));

R = sort(unique(R), 'descend');
%R = sort(R, 1, 'ascend');
for k=R
    if k < 1 | k > nc
   	error('node out of range');
	end

	if k < nc
	   Mat(k,:)=Mat(nc,:);
	end

	nc = nc-1;

	if nc > 0
		Mat=Mat(1:nc,:);
    else
   	    Mat=[];
    end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Elem, Coord] = RemoveNode(Elem, Coord, n2);
ne = length(Elem(:,1));
nc = length(Coord(:,1));

if n2 < 1 | n2 > nc
   error('node out of range');
end

elems = union(find(Elem(:,2)==n2), find(Elem(:,3)==n2));
elems = union(elems, find(Elem(:,4)==n2));
if ~ isempty(elems)
   error('node is in use');
end

if n2 < nc
   Elem = ReplaceNode(Elem,nc, n2);
   Coord(n2,:) = Coord(nc,:);
end

nc = nc-1;

if nc > 0
	Coord=Coord(1:nc,:);
else
   Coord=[];
end
