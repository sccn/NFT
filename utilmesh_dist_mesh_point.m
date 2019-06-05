function [dm,Pm,el,in] = utilmesh_dist_mesh_point(P,Coord,Elem);
% looks for if P is inside the mesh Coord, Elem or not
% Pm is the point on the mesh
% dm is the distance
% el is the element of Pm
% in = inside (bool)
% works for LINEAR MESH


nnp=size(Coord,1);

r=1;
N=0;
Xo=Coord(:,2:4)-ones(nnp,1)*P;
Rad=sum(Xo.*Xo,2);
while length(N)<3
    r=r+5;
    N=find(Rad<r^2);
end

% find the neighbour elements of the closest nodes
E=mesh_elementsofthenodes(Coord,Elem,N);


% find the intersection of the line PP1 with the elements of E
Pint=[]; dis=[];
for i=1:length(E)
   Pa=Coord(Elem(E(i),2),2:4);
   Pb=Coord(Elem(E(i),3),2:4);
   Pc=Coord(Elem(E(i),4),2:4);
   [D,Pp]=warping_disttrianglepoint(Pa,Pb,Pc,P);
   dis(i)=D;
   Pint(i,:)=Pp;
end
[j,k]=min(abs(dis));
Pm=Pint(k,:);
dm=abs(dis(k));
el=E(k);

N2 = Elem(el,2:4);
el2 = mesh_elementsofthenodes(Coord,Elem,N2);

% check if P is inside or outside

% find the normal vector of the element
v1=Coord(Elem(el2,3),2:4)-Coord(Elem(el2,2),2:4);
v2=Coord(Elem(el2,4),2:4)-Coord(Elem(el2,3),2:4);
n=cross(v1, v2); n2=mean(n);
Norm = n2/norm(n2);

% if the vector PPm.Norm > 0 inside, < 0 outside
if dot(Pm-P, Norm) > 0 
    in = 1; % P is inside the mesh
else
    in = 0;
end