% delete mesh, algorithm, or everything at once
% Danil Kirsanov, 09/2007 

function object = geodesic_delete(object)

if nargin == 0          % the simplest way to delete everything is
  clear geodesic;                        % to unload the library
elseif strcmp(object.object_type, 'mesh')       
  geodesic('delete_mesh', object.id);      % delete mesh and corresponding algorithms
else                                        
  geodesic('delete_algorithm', object.id); % delete a single algorithm
end

object = [];
