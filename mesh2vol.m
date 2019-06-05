function vol = mesh2vol(mesh)

vol = [];
ncoordp = 0;
nelemp = 0;

for ii = 1:mesh.num_boundaries
    nelem = mesh.bnd(ii,1) + nelemp;
    vol.bnd(ii).tri = mesh.elem(nelemp+1:nelem,:) - ncoordp;
    ncoord = max(max(vol.bnd(ii).tri)) + ncoordp;
    vol.bnd(ii).pnt = mesh.coord(ncoordp+1:ncoord,:);
    ncoordp = ncoord;
    nelemp = nelem;
end
