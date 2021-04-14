function ss_sparse = ss_cortex_gaus(mesh,algorithm, vertices, faces, stop_distance);
Nn = size(vertices,1);
%ss_sparse = sparse(Nn,Nn);
ss_sparse = sparse([],[],[],Nn,Nn,Nn*1000); % allocate 1000 nonzero spaces

sigma = stop_distance/3; % standard deviation

hh = waitbar(0,'computing patches ...');
for vertex_id = 1:Nn
    waitbar(vertex_id/Nn)

    source_points = {geodesic_create_surface_point('vertex',vertex_id,vertices(vertex_id,:))};

    stop_points = [];

    geodesic_propagate(algorithm, source_points, stop_points, stop_distance); 

    [source_id, distances] = geodesic_distance_and_source(algorithm);     %find distances to all vertices of the mesh; in this example we have a single source, so source_id is always equal to 1
    k1 = find(distances < 1.0000e+100);
    k2 = find(distances == 1.0000e+100);

    max_d = max(distances(k1));
    deg = distances;
    deg(k2) = 0;
    deg(k1) = 1/sigma/sqrt(2*pi) * exp(-0.5/sigma^2*(max_d - distances(k1)).^2);
    deg(k1) = max(deg(k1)) - deg(k1);
    ss_sparse(vertex_id,:) = sparse(deg);
    
end
close(hh)


