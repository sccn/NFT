function [dim, inm] = utilmesh_check_source_space(so,C,E);
% so is the source space
% C,E is the linear mesh
% dim gives the vector of distances
% inm gives the vector of inside sources (logical)

hh = waitbar(0,'computing...');
M = size(so,1);
for i = 1 : M
    waitbar(i/M)
    [dm, Pm, el, in] = utilmesh_dist_mesh_point(so(i,:), C, E);
    dim(i) = dm;
    inm(i) = in;
end
close(hh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

