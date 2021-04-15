function geo_distance = mesh_geodesic_distancep(mesh, algorithm, vertices, faces);

Nn = size(vertices,1);
geo_distance = zeros(Nn,Nn);
hh = waitbar(0, 'computing geodesic distance ...');
parfor i = 1:Nn
    waitbar(i/Nn)
    vertex_from = i;
    source_points = {geodesic_create_surface_point('vertex',vertex_from,vertices(vertex_from,:))};
    %stop_points = {geodesic_create_surface_point('vertex',vertex_to,vertices(vertex_id,:))};
        
    stop_points = [];
    stop_distance = 10000000;
    geodesic_propagate(algorithm, source_points, stop_points, stop_distance); 
    
    %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1...    
    [source_id, distances] = geodesic_distance_and_source(algorithm);     
    geo_distance(:,i) = distances;
end
close(hh)


