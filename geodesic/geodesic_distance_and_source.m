%finds best source and distance to the best source
% if distance is negative, the best source cannot be found (for example, because the propagation was stopped before it reached this point)
% Danil Kirsanov, 09/2007 

function [source_id, distance] = geodesic_distance_and_source(algorithm, destination)

if nargin == 2
  d = geodesic_convert_surface_points({destination});
  [source_id, distance] = geodesic('distance_and_source', algorithm.id, d);
else                                    %return distances and sources for all vertices
  [source_id, distance] = geodesic('distance_and_source', algorithm.id);
end

source_id = source_id + 1;
