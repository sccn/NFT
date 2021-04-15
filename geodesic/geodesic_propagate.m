function geodesic_propagate(algorithm, source_points, stop_points, max_distance)

global geodesic_library;

if nargin < 4
    max_distance = 1e100;
end

if nargin < 3
    stop_points = [];
end

sources = geodesic_convert_surface_points(source_points);
stops = geodesic_convert_surface_points(stop_points);

geodesic('propagate', algorithm.id, sources, stops, max_distance);
