function path = geodesic_trace_back(algorithm, destination)

global geodesic_library;

tmp{1} = destination;
d = geodesic_convert_surface_points(tmp);

path = geodesic_convert_surface_points(geodesic('trace_back', algorithm.id, d));