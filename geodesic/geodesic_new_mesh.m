function [mesh, edge_to_vertex, edge_to_face] = geodesic_new_mesh(points, tri)

dim = find(size(points) == 3);
if dim == 2
    points = points';
end;

dim = find(size(tri) == 3);
if dim == 2
    tri = tri';
end;

[mesh_id, edges] = geodesic('new_mesh', points, tri - 1);

mesh.id = mesh_id;
mesh.object_type = 'mesh';
edge_to_vertex = (edges(1:2,:) + 1)';
edge_to_face = (edges(3:4,:) + 1)';
