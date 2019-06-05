% nft_source_space_generation() - Source space generation
%
% Usage:
%         >> nft_source_space_generation(subject_name, of, 'Key1',Value1,...);
%
% Inputs:
%
%   subject_name : subject name as in main NFT window 
%   of : output folder  
%
% Optional keywords:
%
%   sp : grid spacing - spacing between two dipoles (default = 8 mm)
%   th : minimum distance from the mesh (default = 2 mm)
%   regular : [0 or 1] Regular or symmetric source space (default = 1 (regular))

function nft_source_space_generation(subject_name, of, varargin)

% default values
sp = 8;
th = 2;
regular = 1;

for i = 1:2:length(varargin) % for each Keyword
      Keyword = varargin{i};
      Value = varargin{i+1};

      if ~isstr(Keyword)
         fprintf('keywords must be strings')
         return
      end

      if strcmp(Keyword,'sp')
         if isstr(Value)
            fprintf('sp must be a positive real number');
            return
         else
            sp = Value;
         end
      elseif strcmp(Keyword,'th')
         if isstr(Value)
            fprintf('th must be a positive real number');
            return
         else
             th = Value;
         end
      elseif strcmp(Keyword,'regular')
         if isstr(Value)
            fprintf('regular must be 0 or 1');
            return
         else
             regular = Value;
         end
      end
end

lof = length(of);
if of(lof) ~= '/'
    of(lof+1) = '/';
end


[C1,E1] = mesh_readsmf([of 'Brain.smf'],0,0,0,1); % subject's scalp mesh
% C1 and E1 is the brain mesh

if regular
    [so, ss] = Create_regular_source_space(C1, E1, sp, th);
    save([of subject_name '_sourcespace.dip'], 'so', '-ascii'); 
else
    [so, ss] = Create_symmetric_source_space(C1, E1, sp, th);
    save([of subject_name '_sourcespace.sdip'], 'so', '-ascii'); 
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [so, ss] = Create_regular_source_space(Coord, Elem, spacing, thr);

ma = max(Coord(:,2:4));
mi = min(Coord(:,2:4));
md = round((ma - mi)/spacing);

n=0;
for i=1:md(1)+1
    a(i) = mi(1)+spacing*(i-1)-spacing/2;
    for j=1:md(2)+1
        b(j) = mi(2)+spacing*(j-1)-spacing/2;
        for k=1:md(3)+1
            c(k) = mi(3)+spacing*(k-1)-spacing/2;
            n=n+1;
            rw(n,1:3)=[a(i) b(j) c(k)];
        end
    end
end

[dim, inm] = utilmesh_check_source_space(rw, Coord, Elem);

k = find(inm == 1);   % dipoles inside the mesh
l = find(dim < thr);  % dipoles closer to the mesh less than thr
m = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr
rw = rw(m,:);

ne=size(rw,1);    
ss=zeros(ne*3,7);
ss(1:ne,2:4)=rw;
ss(ne+1:2*ne,2:4)=rw;
ss(2*ne+1:3*ne,2:4)=rw;
ss(1:ne,5)=1;
ss(ne+1:2*ne,6)=1;
ss(2*ne+1:3*ne,7)=1;
ss(:,1)=[1:3*ne]';

Ns = size(rw, 1);
so = zeros(3*Ns, 6);
so(1:Ns, 1:3) = rw;
so(1:Ns, 4) = 1;
so(1+Ns:2*Ns, 1:3) = rw;
so(1+Ns:2*Ns, 5) = 1;
so(1+Ns*2:3*Ns, 1:3) = rw;
so(1+Ns*2:3*Ns, 6) = 1;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [so, ss] = Create_symmetric_source_space(Coord, Elem, spacing, thr);

ma = max(Coord(:,2:4));
mi = min(Coord(:,2:4));
md = round((ma - mi)/spacing);

n=0;
for i=1:md(1)+1
    a(i) = mi(1)+spacing*(i-1)-spacing/2;
    for j=1:md(2)+1
        b(j) = mi(2)+spacing*(j-1)-spacing/2;
        for k=1:md(3)+1
            c(k) = mi(3)+spacing*(k-1)-spacing/2;
            n=n+1;
            ss(n,1:3)=[a(i) b(j) c(k)];
        end
    end
end

% symmetry according to y axis
smy = (max(Coord(:,3))+min(Coord(:,3)))/2;

ksmy = find(ss(:,2) > smy);
ss1 = ss(ksmy,:);
ss2 = ss1;
ss2(:,2) = -ss2(:,2);
ss2(:,2) = ss2(:,2)+smy*2;

[dim1, inm1] = utilmesh_check_source_space(ss1, Coord, Elem);
[dim2, inm2] = utilmesh_check_source_space(ss2, Coord, Elem);


k = find(inm1 == 1);   % dipoles inside the mesh
l = find(dim1 < thr);  % dipoles closer to the mesh less than thr
m1 = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr

k = find(inm2 == 1);   % dipoles inside the mesh
l = find(dim2 < thr);  % dipoles closer to the mesh less than thr
m2 = setdiff(k, l);    % dipoles inside the mesh, closer to the mesh less than thr

mi = intersect(m1,m2);

ss1 = ss1(mi,:);
ss2 = ss2(mi,:);

% first Ns y directed dipoles
Ns = size(ss1, 1);
so1 = zeros(3*Ns, 6);
so1(1:Ns, 1:3) = ss1;
so1(1:Ns, 4) = 1;
so1(1+Ns:2*Ns, 1:3) = ss1;
so1(1+Ns:2*Ns, 5) = 1;
so1(1+Ns*2:3*Ns, 1:3) = ss1;
so1(1+Ns*2:3*Ns, 6) = 1;

% second Ns dipoes are the symmetric ones
so2 = zeros(3*Ns, 6);
so2(1:Ns, 1:3) = ss2;
so2(1:Ns, 4) = 1;
so2(1+Ns:2*Ns, 1:3) = ss2;
so2(1+Ns:2*Ns, 5) = -1;
so2(1+Ns*2:3*Ns, 1:3) = ss2;
so2(1+Ns*2:3*Ns, 6) = 1;

so = [so1;so2]; 
