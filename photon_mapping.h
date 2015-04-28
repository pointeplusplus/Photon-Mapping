#ifndef _PHOTON_MAPPING_H_
#define _PHOTON_MAPPING_H_

#include <mutex>
#include <fstream>
#include <thread>
#include <vector>


#include "vectors.h"
#include "photon.h"
#include "vbo_structs.h"
#include "face.h"

class Material;
class Mesh;
class ArgParser;
class KDTree;
class Ray;
class Hit;
class RayTracer;
class Radiosity;

// =========================================================================
// The basic class to shoot photons within the scene and collect and
// process the nearest photons for use in the raytracer

class LineSegment {

public:
	// CONSTRUCTOR
	LineSegment(Vec3f start_loc, Vec3f end_loc, Vec4f col) {
		// first clamp the segment to "reasonable" values 
		// to make sure it is drawn correctly in OpenGL
		a = start_loc;
		b = end_loc; 
		color = col;
	}
	Vec4f getColor() const { return color; }
	const Vec3f& getStart() const { return a; }
	const Vec3f& getEnd() const { return b; }
private:
	// REPRESENTATION
	Vec3f a;
	Vec3f b;
	Vec4f color;
};

class PhotonMapping {

 public:

	// CONSTRUCTOR & DESTRUCTOR
	PhotonMapping(Mesh *_mesh, ArgParser *_args) {
		mesh = _mesh;
		args = _args;
		raytracer = NULL;
		kdtree = NULL;
		last_radius = args->default_radius;
	}
	~PhotonMapping();
	void setRayTracer(RayTracer *r) { raytracer = r; }
	void setRadiosity(Radiosity *r) { radiosity = r; }

	void initializeVBOs(); 
	void setupVBOs(); 
	void drawVBOs();
	void cleanupVBOs();

	void printEscapingFacePhoton();
	void printOutputFile();

	// step 1: send the photons throughout the scene
	void TracePhotons();
	void TracePhotonsWorker(
		const std::vector<Face*>& lights, 
		double total_lights_area,
		int num_threads
	);
	
	// step 2: collect the photons and return the contribution from indirect illumination
	Vec3f GatherIndirect(const Vec3f &point, const Vec3f &normal, const Vec3f &direction_from);

	bool CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const;

	// trace a single photon
	void TracePhoton(const Vec3f &position, const Vec3f &direction, const float wavelength, int iter, Vec4f viz_color, Material* current_material, float current_n_val, bool single_photon);

 private:

	

	// REPRESENTATION
	KDTree *kdtree;
	Mesh *mesh;
	ArgParser *args;
	RayTracer *raytracer;
	Radiosity *radiosity;

	std::vector<LineSegment> visualization_line_segments;
	
	// For gathering
	double last_radius;
	std::mutex kdtree_lock;
	std::mutex ray_tree_lock;
	
	// VBO
	GLuint photon_verts_VBO;
	GLuint photon_direction_indices_VBO;
	GLuint kdtree_verts_VBO;
	GLuint kdtree_edge_indices_VBO;
	std::vector<VBOPosColor> photon_verts; 
	std::vector<VBOIndexedEdge> photon_direction_indices;
	std::vector<VBOPos> kdtree_verts; 
	std::vector<VBOIndexedEdge> kdtree_edge_indices;

	GLuint raytree_verts_VBO;
	GLuint raytree_edge_indices_VBO;
	std::vector<VBOPosColor4> raytree_verts; 
	std::vector<VBOIndexedEdge> raytree_edge_indices;
};

// =========================================================================

#endif
