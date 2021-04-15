#ifndef GEODESIC_DLL_HEADER_HPP_234232
#define GEODESIC_DLL_HEADER_HPP_234232

#ifndef GEODESIC_DLL_IMPORT 
#define GEODESIC_DLL_IMPORT
#endif

#ifdef __cplusplus
extern "C" {
#endif

GEODESIC_DLL_IMPORT int new_mesh(int num_points,		//creates new mesh
								  double* points,	
								  int num_triangles,
								  int* triangles, 
								  int* num_edges, 
								  double** edges);

GEODESIC_DLL_IMPORT int new_algorithm(int mesh_id,	//creates a geodesic algorithm for a given mesh
						               int type,
									   int subdivision);

GEODESIC_DLL_IMPORT void delete_algorithm(int id);

GEODESIC_DLL_IMPORT void delete_mesh(int id);			//delete mesh and all associated algorithms

GEODESIC_DLL_IMPORT void propagate(int algorithm_id,		//compute distance field for given source points
									double* source_points,	
									int num_sources,
									double* stop_points,	//limitations on distance field propagation
									int num_stop_points,
									double max_propagation_distance);

GEODESIC_DLL_IMPORT int trace_back(int algorithm_id,		//using procomputed distance field, compute a shortest path from destination to the closest source
									double* destination,
									double** path);

GEODESIC_DLL_IMPORT int distance_and_source(int algorithm_id,		//quickly find what source this point beints to and what is the distance to this source
											 double* destination,			
											 double* best_source_distance);

GEODESIC_DLL_IMPORT int distance_and_source_for_all_vertices(int algorithm_id,	//same idea as in the previous function
															  double** distances,	//list distance/source info for all vertices of the mesh
															  int** sources);

#ifdef __cplusplus
}
#endif

#endif
