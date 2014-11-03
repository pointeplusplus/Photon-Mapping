#ifndef _RAY_TREE_H
#define _RAY_TREE_H

#include "glCanvas.h"
#include <vector>
#include "ray.h"
#include "vbo_structs.h"

// ====================================================================
// ====================================================================
// data structure to store a segment

class Segment {

public:
	// CONSTRUCTOR
	Segment(const Ray &ray, double tstart, double tstop) {
		// first clamp the segment to "reasonable" values 
		// to make sure it is drawn correctly in OpenGL
		if (tstart < -1000) tstart = -1000;
		if (tstop	>	1000) tstop	=	1000;
		a = ray.pointAtParameter(tstart);
		b = ray.pointAtParameter(tstop); }
	const Vec3f& getStart() const { return a; }
	const Vec3f& getEnd() const { return b; }
private:
	// REPRESENTATION
	Vec3f a;
	Vec3f b;
};

// ====================================================================
// ====================================================================
//
// This class only contains static variables and static member
// functions so there is no need to call the constructor, destructor
// etc.	It's just a wrapper for the ray tree visualization data.
//

class RayTree {

public:

	// most of the time the RayTree is NOT activated, so the segments are not updated
	static void Activate() { Clear(); activated = 1; }
	static void Deactivate() { activated = 0; }

	// when activated, these function calls store the segments of the tree
	static void AddMainSegment(const Ray &ray, double tstart, double tstop) {
		if (!activated) return;
		main_segments.push_back(Segment(ray,tstart,tstop));
	}
	static void AddShadowSegment(const Ray &ray, double tstart, double tstop) {
		if (!activated) return;
		shadow_segments.push_back(Segment(ray,tstart,tstop));
	}
	static void AddReflectedSegment(const Ray &ray, double tstart, double tstop) {
		if (!activated) return;
		reflected_segments.push_back(Segment(ray,tstart,tstop));
	}
	static void AddTransmittedSegment(const Ray &ray, double tstart, double tstop) {
		if (!activated) return;
		transmitted_segments.push_back(Segment(ray,tstart,tstop));
	}

	static void AddGeneralSegment(const Ray &ray, double tstart, double tstop, Vec4f color) {
		if (!activated) return;
		general_segments.push_back(Segment(ray,tstart,tstop));
		general_segment_colors.push_back(color);
		std::cout << "number of segments so far: " << general_segments.size() << std::endl;
	}

	static void initializeVBOs(); 
	static void setupVBOs(); 
	static void drawVBOs();
	static void cleanupVBOs();
	
private:
	
	// HELPER FUNCTIONS
	static void Clear() {
//		std::cout << "Clearing segments!" << std::endl;
//		std::cout << "    there were " << general_segments.size() << " general segments." << std::endl;
		main_segments.clear();
		shadow_segments.clear();
		reflected_segments.clear();
		transmitted_segments.clear();
		general_segments.clear();
		general_segment_colors.clear();
	}
	
	// REPRESENTATION
	static int activated;
	static std::vector<Segment> main_segments;
	static std::vector<Segment> shadow_segments;
	static std::vector<Segment> reflected_segments;
	static std::vector<Segment> transmitted_segments;
	static std::vector<Segment> general_segments;
	static std::vector<Vec4f> general_segment_colors;

	// VBO
	static GLuint raytree_verts_VBO;
	static GLuint raytree_edge_indices_VBO;
	static std::vector<VBOPosColor4> raytree_verts; 
	static std::vector<VBOIndexedEdge> raytree_edge_indices;
};

// ====================================================================
// ====================================================================

#endif
