#ifndef _FACE_H_
#define _FACE_H_

#include "edge.h"
#include "ray.h"
#include "vertex.h"
#include "hit.h"
#include <vector> 
class Material;
class Hit;

// ==============================================================
// Simple class to store quads for use in radiosity & raytracing.

class Face {

public:

	// ========================
	// CONSTRUCTOR & DESTRUCTOR
	Face(Material *m) {
		edge = NULL;
		material = m;
		num_rays_leaving_face = 0;	
		num_rays_entering_face = 0;	 
		num_interior_bounces = 0;
		num_rays_reflected = 0;
		cached_normal = NULL;
	}

	// =========
	// ACCESSORS
	Vertex* operator[](int i) const { 
		assert (edge != NULL);
		if (i==0) return edge->getStartVertex();
		if (i==1) return edge->getNext()->getStartVertex();
		if (i==2) return edge->getNext()->getNext()->getStartVertex();
		if (i==3) return edge->getNext()->getNext()->getNext()->getStartVertex();
		assert(0);
		exit(0);
	}
	Edge* getEdge() const { 
		assert (edge != NULL);
		return edge; 
	}
	Vec3f computeCentroid() const {
		return 0.25 * ((*this)[0]->get() +
									 (*this)[1]->get() +
									 (*this)[2]->get() +
									 (*this)[3]->get());
	}
	Material* getMaterial() const { return material; }
	double getArea() const;

	Vec3f RandomPoint() const;
	//Rebecca defined
	Vec3f RandomStratifiedPoint(int x_iter, int y_iter, int num_samples) const;
	Vec3f computeNormal();
	int getNumRaysLeavingFace() const { return num_rays_leaving_face; }
	int getNumRaysEnteringFace() const {return num_rays_entering_face; }
	int getNumInteriorBounces() const {return num_interior_bounces; }
	int getNumRaysReflected() const {return num_rays_reflected; }
	const std::vector<Vec3f>& getLeavingDirections() const {return light_leaving_directions;}

	// =========
	// MODIFIERS
	void setEdge(Edge *e) {
		assert (edge == NULL);
		assert (e != NULL);
		edge = e;
	}
	void incrementNumRaysLeaving(){ 
		num_rays_leaving_face++;
	}

	void incrementNumRaysEntering(){ 
		num_rays_entering_face++;
	}

	void incrementNumInteriorBounces(){ 
		num_interior_bounces++;
	}

	void incrementNumRaysReflected(){
		num_rays_reflected++;
	}

	void addLightLeavingDirection(Vec3f leaving_direction){
		light_leaving_directions.push_back(leaving_direction);
	}

	//Debug Functions
	bool isConvex() const;
	void printVertices() const;
	double shortestEdge() const;


 
	// ==========
	// RAYTRACING
	bool intersect(const Ray &r, Hit &h, bool intersect_backfacing, bool* backfacing_hit);

	// =========
	// RADIOSITY
	int getRadiosityPatchIndex() const { return radiosity_patch_index; }
	void setRadiosityPatchIndex(int i) { radiosity_patch_index = i; }

protected:

	// helper functions
	bool triangle_intersect(const Ray &r, Hit &h, Vertex *a, Vertex *b, Vertex *c, bool intersect_backfacing, bool* backfacing_hit);
	bool plane_intersect(const Ray &r, Hit &h, bool intersect_backfacing, bool* backfacing_hit);

	void computeCachedNormal();

	// don't use this constructor
	Face& operator= (const Face&) { assert(0); exit(0); }
	
	// ==============
	// REPRESENTATION
	Edge *edge;
	// NOTE: If you want to modify a face, remove it from the mesh,
	// delete it, create a new copy with the changes, and re-add it.
	// This will ensure the edges get updated appropriately.
	
	int radiosity_patch_index;	// an awkward pointer to this patch in the Radiosity patch array
	Material *material;

	//Rebecca Added
	int num_rays_leaving_face;
	int num_rays_entering_face;
	int num_interior_bounces;
	int num_rays_reflected;
	std::vector<Vec3f> light_leaving_directions;
	Vec3f* cached_normal; // cache normal
};

// ===========================================================

#endif
