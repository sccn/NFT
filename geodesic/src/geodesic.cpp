#include <iostream>
#include <fstream>
#include <string>
#include <memory>

#include "mex.h"

#define Xassert(expr) ((expr) ? (void) (0) : Xassert_fail(#expr, __FILE__, __LINE__))

void Xassert_fail(const char *, const char *, unsigned int);

#include "geodesic_mesh.h"
#include "geodesic_algorithm_dijkstra_alternative.h"
#include "geodesic_algorithm_dijkstra.h"
#include "geodesic_algorithm_subdivision.h"
#include "geodesic_algorithm_exact.h"

#ifdef _MSC_VER 
#define strcasecmp _stricmp
#endif

extern "C" {
size_t strlcpy(char *dst, const char *src, size_t siz);
}

typedef std::shared_ptr<geodesic::Mesh> mesh_shared_pointer;
std::vector<mesh_shared_pointer> meshes;

typedef std::shared_ptr<geodesic::GeodesicAlgorithmBase> algorithm_shared_pointer;
std::vector<algorithm_shared_pointer> algorithms;

void gd_help( int nlhs, mxArray *plhs[], int nrhs, const mxArray*prhs[]);

void
Xassert_fail(const char *expr, const char *file, unsigned int line)
{
	char buf[256];
	snprintf(buf, sizeof(buf), "%s:%u assert(%s) failed!", file, line, expr);
	mexErrMsgTxt(buf);
}

static
std::size_t find_mesh_id(geodesic::Mesh* mesh)
{
	for(std::size_t i=0; i<meshes.size(); ++i) {
		if(meshes[i].get() == mesh)
			return i;
	}
	return -1;
}

// [mesh_id, edges] = geodesic('new_mesh', points, triangles);
void
gd_new_mesh(int nlhs, mxArray *plhs[], 
	    int nrhs, const mxArray*prhs[])
{
	const mwSize *size;
	double *points, *triangles, *edges;
	long num_points, num_triangles, num_edges;

	if (nrhs != 2)
		mexErrMsgTxt("expecting two input arguments");

	if (nlhs > 2)
		mexErrMsgTxt("expecting at most two output arguments");

	if ( ! mxIsDouble(prhs[0]) || mxGetNumberOfDimensions(prhs[0]) != 2 )
		mexErrMsgTxt("expecting points (arg 1) to be a matrix");

	if ( ! mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 )
		mexErrMsgTxt("expecting tri (arg 1) to be a matrix");

	size = mxGetDimensions(prhs[0]);
	if (size[0] != 3)
		mexErrMsgTxt("points: invalid size (expected 3 columns)");

	num_points = size[1];
	points = mxGetPr(prhs[0]);

	size = mxGetDimensions(prhs[1]);
	if (size[0] != 3)
		mexErrMsgTxt("tri: invalid size (expected 3 columns)");
	num_triangles = size[1];
	triangles = mxGetPr(prhs[1]);
	
	mesh_shared_pointer new_mesh = mesh_shared_pointer(new geodesic::Mesh);
	meshes.push_back(new_mesh);

	new_mesh->initialize_mesh_data(num_points,
				       points,	
				       num_triangles,
				       triangles);

	// mesh_id
	plhs[0] = mxCreateDoubleScalar(meshes.size() - 1);

	if (nlhs < 2 )
		return;

	num_edges = new_mesh->edges().size();

	plhs[1] = mxCreateDoubleMatrix(4, num_edges, mxREAL);

	edges = mxGetPr(plhs[1]);

	if (edges == NULL)
		mexErrMsgTxt("Failed to allocate edges!");

	for(std::size_t i=0; i<num_edges; ++i) {
		geodesic::Edge& edge = new_mesh->edges()[i];
		double* buffer = edges + 4*i;

		buffer[0] = edge.adjacent_vertices()[0]->id();
		buffer[1] = edge.adjacent_vertices()[1]->id();
		buffer[2] = edge.adjacent_faces()[0]->id();
		buffer[3] = edge.adjacent_faces().size() == 2 ? 
			edge.adjacent_faces()[1]->id() : -1; 
	}
}

// algorithm_id = geodesic('new_algorithm', mesh_id, type_id, subdivision);
void
gd_new_algorithm(int nlhs, mxArray *plhs[], 
		 int nrhs, const mxArray*prhs[])
{
	long mesh_id, type_id, subdivision;

	if (nrhs != 3)
		mexErrMsgTxt("expecting three input arguments");

	if (nlhs != 1)
		mexErrMsgTxt("expecting one output argument");

	if (! (mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting mesh_id (arg 1) to be a scalar");

	if (! (mxIsNumeric(prhs[1]) && mxIsScalar(prhs[1])))
		mexErrMsgTxt("expecting type_id (arg 2) to be a scalar");

	if (! (mxIsNumeric(prhs[2]) && mxIsScalar(prhs[2])))
		mexErrMsgTxt("expecting subdividion (arg 3) to be a scalar");

	mesh_id = mxGetScalar(prhs[0]);
	type_id = mxGetScalar(prhs[1]);
	subdivision = mxGetScalar(prhs[2]);

	if (mesh_id < 0 || mesh_id >= meshes.size())
		mexErrMsgTxt("invalid mesh_id");

	if (subdivision == 0 && type_id == 1)
		type_id = 2;

	geodesic::Mesh* mesh = meshes[mesh_id].get();
	geodesic::GeodesicAlgorithmBase* algorithm;

	switch(type_id) {
		case 2: //DIJKSTRA
		{
			algorithm = new geodesic::GeodesicAlgorithmDijkstra(mesh);
			break;
		}
		case 1: //SUBDIVISION:
		{
			algorithm = new geodesic::GeodesicAlgorithmSubdivision(mesh, subdivision);
			break;
		}
		default: // EXACT
		{
			algorithm = new geodesic::GeodesicAlgorithmExact(mesh);
			break;
		}
	}

	algorithms.push_back(algorithm_shared_pointer(algorithm));

	// mesh_id
	plhs[0] = mxCreateDoubleScalar(algorithms.size() - 1);
}

// geodesic('delete_mesh', mesh_id);
void
gd_delete_mesh(int nlhs, mxArray *plhs[], 
	       int nrhs, const mxArray*prhs[])
{
	long mesh_id;

	if (nrhs != 1)
		mexErrMsgTxt("expecting one input argument");

	if (! (mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting mesh_id (arg 1) to be a scalar");

	mesh_id = mxGetScalar(prhs[0]);

	if (mesh_id < 0 || mesh_id >= meshes.size())
		mexErrMsgTxt("invalid mesh_id");

	geodesic::Mesh* mesh = meshes[mesh_id].get();
	for(std::size_t i=0; i<algorithms.size(); ++i)
	{
		geodesic::GeodesicAlgorithmBase* algorithm = algorithms[i].get();
		if(algorithm && algorithm->mesh() == mesh)
		{
			algorithms[i] = algorithm_shared_pointer();
		}
	}

	meshes[mesh_id] = mesh_shared_pointer();
}

// geodesic('delete_algorithm', algorithm_id);
void
gd_delete_algorithm(int nlhs, mxArray *plhs[], 
		    int nrhs, const mxArray*prhs[])
{
	long id;

	if (nrhs != 1)
		mexErrMsgTxt("expecting one input argument");

	if (! (mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting algorithm_id (arg 1) to be a scalar");

	id = mxGetScalar(prhs[0]);

	if (id < 0 || id >= algorithms.size())
		mexErrMsgTxt("invalid algorithm_id");

	algorithms[id] = algorithm_shared_pointer();
}

// geodesic('propagate', algorithm, source_points, stop_points, max_distance);
void
gd_propagate(int nlhs, mxArray *plhs[], 
	    int nrhs, const mxArray*prhs[])
{
	const mwSize *size;
	mwSize num_sources, num_stop_points;
	double *source_points, *stop_points = NULL, max_distance;
	long id;

	if (nrhs != 4)
		mexErrMsgTxt("expecting four input arguments");

	if (! (mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting algorithm_id (arg 1) to be a scalar");

	if (! mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 )
		mexErrMsgTxt("expecting source_points (arg 2) to be a matrix/vector");
	else {
		size = mxGetDimensions(prhs[1]);
		num_sources = size[0] * size[1];
		if ((num_sources % 5) != 0)
			mexErrMsgTxt("source_points: invalid size (expected multiple of 5 )");
		num_sources = num_sources / 5;
	}

	if (mxIsEmpty(prhs[2]))
		num_stop_points = 0;
	else if (! mxIsDouble(prhs[2]) || mxGetNumberOfDimensions(prhs[2]) != 2 )
		mexErrMsgTxt("expecting stop_points (arg 3) to be a matrix/vector");
	else {
		size = mxGetDimensions(prhs[2]);
		num_stop_points = size[0] * size[1];
		if ((num_stop_points % 5) != 0)
			mexErrMsgTxt("stop_points: invalid size (expected multiple of 5 )");
		num_stop_points = num_stop_points / 5;
	}

	if (! (mxIsNumeric(prhs[3]) && mxIsScalar(prhs[3])))
		mexErrMsgTxt("expecting max_distance (arg 4) to be a scalar");

	id = mxGetScalar(prhs[0]);
	source_points = mxGetPr(prhs[1]);
	if (num_stop_points > 0)
		stop_points = mxGetPr(prhs[2]);

	max_distance = mxGetScalar(prhs[3]);

	if (id < 0 || id >= algorithms.size())
		mexErrMsgTxt("invalid algorithm_id");

	std::vector<geodesic::SurfacePoint> sources(num_sources);

	geodesic::Mesh* mesh = algorithms[id]->mesh();
	for(std::size_t i=0; i<num_sources; ++i)
	{
		geodesic::fill_surface_point_structure(&sources[i], 
						       source_points + 5*i, 
						       mesh);
	}

	std::vector<geodesic::SurfacePoint> stop(num_stop_points);
	for(std::size_t i=0; i<num_stop_points; ++i)
	{
		geodesic::fill_surface_point_structure(&stop[i], 
						       stop_points + 5*i, 
						       mesh);
	}

	algorithms[id]->propagate(sources, 
					    max_distance,
					    &stop);
}

// path = geodesic('trace_back', algorithm_id, destination);
void
gd_trace_back(int nlhs, mxArray *plhs[], 
	      int nrhs, const mxArray*prhs[])
{
	const mwSize *size;
	double *destination, *path;
	long id;

	if (nrhs != 2)
		mexErrMsgTxt("expecting two input arguments");

	if (nlhs != 1)
		mexErrMsgTxt("expecting one output argument");

	if (!(mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting algorithm_id (arg 1) to be a scalar");

	id = mxGetScalar(prhs[0]);
	if (id < 0 || id >= algorithms.size())
		mexErrMsgTxt("invalid algorithm_id");

	if (! mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 )
		mexErrMsgTxt("expecting destination (arg 2) to be a matrix/vector");

	size = mxGetDimensions(prhs[1]);
	long num_dest = size[0] * size[1];
	if (num_dest != 5)
		mexErrMsgTxt("destination: must be a single point");
	destination = mxGetPr(prhs[1]);

	geodesic::SurfacePoint point;
	geodesic::GeodesicAlgorithmBase* algorithm = algorithms[id].get();
	std::vector<geodesic::SurfacePoint> output_path;

	geodesic::fill_surface_point_structure(&point, 
					       destination, 
					       algorithm->mesh());

	algorithm->trace_back(point, output_path);

	std::size_t mesh_id = find_mesh_id(algorithm->mesh());

	plhs[0] = mxCreateDoubleMatrix(1, output_path.size()*5, mxREAL);

	path = mxGetPr(plhs[0]);

	for(std::size_t i=0; i<output_path.size(); ++i)
		geodesic::fill_surface_point_double(&output_path[i], path + 5 * i, mesh_id);
}

// [source_id, distances] = geodesic('distance_and_source', algorithm_id, destination);
void
gd_distance_and_source(int nlhs, mxArray *plhs[], 
		       int nrhs, const mxArray*prhs[])
{
	const mwSize *size;
	double *destination, *path;
	long id;

	if (nrhs < 1)
		mexErrMsgTxt("expecting at least one input argument");

	if (nrhs > 2)
		mexErrMsgTxt("expecting at most two input arguments");

	if (nlhs > 2)
		mexErrMsgTxt("expecting at most two output arguments");

	if (!(mxIsNumeric(prhs[0]) && mxIsScalar(prhs[0])))
		mexErrMsgTxt("expecting algorithm_id (arg 1) to be a scalar");

	id = mxGetScalar(prhs[0]);
	if (id < 0 || id >= algorithms.size())
		mexErrMsgTxt("invalid algorithm_id");

	geodesic::GeodesicAlgorithmBase* algorithm = algorithms[id].get();
	geodesic::Mesh* mesh = algorithm->mesh();

	if (nrhs < 2) { // for all vertices
		plhs[0] = mxCreateDoubleMatrix(1, mesh->vertices().size(), mxREAL);
		plhs[1] = mxCreateDoubleMatrix(1, mesh->vertices().size(), mxREAL);
		double *sources = mxGetPr(plhs[0]);
		double *distances = mxGetPr(plhs[1]);

		for(std::size_t i = 0; i < mesh->vertices().size(); ++i) {
			geodesic::SurfacePoint point(&mesh->vertices()[i]);
			sources[i] = algorithm->best_source(point, distances[i]);
		}
	} else {
		// for destination vertex
		if (! mxIsDouble(prhs[1]) || mxGetNumberOfDimensions(prhs[1]) != 2 )
			mexErrMsgTxt("expecting destination (arg 2) to be a matrix/vector");

		size = mxGetDimensions(prhs[1]);
		long num_dest = size[0] * size[1];
		if (num_dest != 5)
			mexErrMsgTxt("destination: must be a single point");
		destination = mxGetPr(prhs[1]);

		geodesic::SurfacePoint point;
		geodesic::fill_surface_point_structure(&point, destination, mesh);
		double best_source_distance;
		std::size_t best_source = algorithm->best_source(point, best_source_distance);

		plhs[0] = mxCreateDoubleScalar(best_source);
		plhs[1] = mxCreateDoubleScalar(best_source_distance);
	}
}

struct command {
	const char *cm_name;
	const char *cm_proto;
	const char *cm_desc;
	void  (*cm_func)(int, mxArray *[], int, const mxArray*[]);
};


struct command gd_commands[] = {
	{"new_mesh", "[mesh_id, edges] = geodesic('new_mesh', points, triangles)",
	 "creates new mesh", gd_new_mesh},
	{"new_algorithm", "algorithm_id = geodesic('new_algorithm', mesh_id, type_id, subdivision)",
	 "creates a geodesic algorithm for a given mesh", gd_new_algorithm},
	{"delete_mesh", "geodesic('delete_mesh', mesh_id)",
	 "delete mesh and all associated algorithms", gd_delete_mesh},
	{"delete_algorithm","geodesic('delete_algorithm', algorithm_id)",
	 "deletes a given algorithm", gd_delete_algorithm},
	{"propagate", "geodesic('propagate', algorithm, source_points, stop_points, max_distance)",
	 "compute distance field for given source points", gd_propagate}, 
	{"trace_back", "path = geodesic('trace_back', algorithm_id, destination)",
	 "using procomputed distance field, compute a shortest path from destination to the closest source",
	 gd_trace_back},
	{"distance_and_source",
	 "[source_id, distances] = geodesic('distance_and_source', algorithm_id, destination)",
	 "quickly find what source this point beints to and what is the distance to this source",
	 gd_distance_and_source},
	{"help", "geodesic('help')", "show help text", gd_help},
	{"", "", "", NULL}};

void
gd_help( int nlhs, mxArray *plhs[], 
	int nrhs, const mxArray*prhs[])
{
	mexPrintf("usage: geodesic('command', ...);\n\n");
	mexPrintf("where 'command' is:\n");

	for (struct command *c = gd_commands; c->cm_func; c++) {
		mexPrintf("  '%s': %s\n", c->cm_name, c->cm_desc);
		mexPrintf("    %s\n", c->cm_proto);
	}
}


void
mexFunction( int nlhs, mxArray *plhs[], 
	     int nrhs, const mxArray*prhs[] )
{
	static int need_setup = 1;
	const char *errmsg = NULL;
	struct command *c;
	char cmd[32], *cp;

	cp = NULL;

	if (nrhs < 1)
		mexErrMsgTxt("Usage: geodesic('command', ...)");

	cp = mxArrayToString(prhs[0]);

	if (cp == NULL) {
		errmsg = "expecting command (arg 1) to be a string";
		goto exit;
	}
	if (strlcpy(cmd, cp, sizeof(cmd)) >= sizeof(cmd)) {
		errmsg = "invalid command, use 'help' for a list.";
		mxFree(cp);
		goto exit;
	}
	
	mxFree(cp);

	prhs++;
	nrhs--;

	for (c = gd_commands; c->cm_func; c++) {
		if (strcasecmp(cmd, c->cm_name) == 0) {
			c->cm_func(nlhs, plhs, nrhs, prhs);
			goto exit;
		}
	}

	errmsg = "unknown command, use 'help' for a list.";

exit:
	if (errmsg)
		mexErrMsgTxt(errmsg);
}
