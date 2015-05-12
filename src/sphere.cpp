#include "utils.h"
#include "material.h"
#include "argparser.h"
#include "sphere.h"
#include "vertex.h"
#include "mesh.h"
#include "ray.h"
#include "hit.h"
#include <cmath>

#define eps 0.1
// ====================================================================
// ====================================================================

bool Sphere::intersect(const Ray &r, Hit &h) const {
//	std::cout << " current t value " << h.getT() << std::endl;; 

	// ==========================================
	// IMPLEMENT SPHERE INTERSECTION
	// ==========================================

	// plug the explicit ray equation into the implict sphere equation and solve



	// return true if the sphere was intersected, and update
	// the hit data structure to contain the value of t for the ray at
	// the intersection point, the material, and the normal
	Vec3f ray_origin = r.getOrigin();
	Vec3f ray_direction = r.getDirection();

	double b = 2.0 * (ray_direction.Dot3(ray_origin - center));
	double c = (ray_origin - center).Dot3(ray_origin - center) - pow(radius, 2.0);

	double d = (pow(b,2.0) - 4.0 * c);

	//case where roots are imaginary
	if(d < 0){
		return false;
	}

//	std::cout << "b: " << b << "	c: " << c << std::endl; 
	double root1 = ((-1)* b + sqrt(d))/2.0;
	double root2 = ((-1)* b - sqrt(d))/2.0;
//	std::cout << "root1: " << root1 << "	root2: " << root2 << std::endl; 

	Vec3f root1_normal = ray_origin + root1*ray_direction;
	root1_normal -=center;
	root1_normal.Normalize();
	Vec3f root2_normal = ray_origin + root2*ray_direction;
	root2_normal -=center;
	root2_normal.Normalize();


	//sphere is behind us
	if(root1 < 0 && root2 < 0){
//		std::cout << "both are negative" << std::endl;
		return false;
	}
	//only root1 is in front
	if(root1 >= eps && root2 < 0 && (h.getT() < eps || root1 < h.getT())) {
//		if(root1 < h.get_t()){
			h.set(root1, material, root1_normal,NULL);
//		}
		
		return true;
	}
	//only root2 is in front
	if(root2 >= eps && root1 < 0 && (h.getT() < eps || root1 < h.getT())) {
//		if(root2 < h.get_t()){
			h.set(root2, material, root2_normal,NULL);
//		}
		
		return true;
	}
	//both are in front
	if(root1 >= eps && root2 >= eps && (root1 < h.getT() || root2 < h.getT() || h.getT() < eps) ) {
		if(root1 < root2 ){

			h.set(root1, material, root1_normal,NULL);
		}
		else {
			h.set(root2, material, root2_normal,NULL);
		}

		return true;
	}


	//If it gets to this case it means that we did hit the sphere, but there was already a closer hit.
	//std::cout << "All cases exhausted. This shouldn't happen" << std::endl;
	return false;
} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
Vec3f ComputeSpherePoint(double s, double t, const Vec3f center, double radius) {
	double angle = 2*M_PI*s;
	double y = -cos(M_PI*t);
	double factor = sqrt(1-y*y);
	double x = factor*cos(angle);
	double z = factor*-sin(angle);
	Vec3f answer = Vec3f(x,y,z);
	answer *= radius;
	answer += center;
	return answer;
}

void Sphere::addRasterizedFaces(Mesh *m, ArgParser *args) {
	
	// and convert it into quad patches for radiosity
	int h = args->sphere_horiz;
	int v = args->sphere_vert;
	assert (h % 2 == 0);
	int i,j;
	int va,vb,vc,vd;
	Vertex *a,*b,*c,*d;
	int offset = m->numVertices(); //vertices.size();

	// place vertices
	m->addVertex(center+radius*Vec3f(0,-1,0));	// bottom
	for (j = 1; j < v; j++) { // middle
		for (i = 0; i < h; i++) {
			double s = i / double(h);
			double t = j / double(v);
			m->addVertex(ComputeSpherePoint(s,t,center,radius));
		}
	}
	m->addVertex(center+radius*Vec3f(0,1,0)); // top

	// the middle patches
	for (j = 1; j < v-1; j++) {
		for (i = 0; i < h; i++) {
			va = 1 +	i	 + h*(j-1);
			vb = 1 + (i+1)%h + h*(j-1);
			vc = 1 +	i	 + h*(j);
			vd = 1 + (i+1)%h + h*(j);
			a = m->getVertex(offset + va);
			b = m->getVertex(offset + vb);
			c = m->getVertex(offset + vc);
			d = m->getVertex(offset + vd);
			m->addRasterizedPrimitiveFace(a,b,d,c,material);
		}
	}

	for (i = 0; i < h; i+=2) {
		// the bottom patches
		va = 0;
		vb = 1 +	i;
		vc = 1 + (i+1)%h;
		vd = 1 + (i+2)%h;
		a = m->getVertex(offset + va);
		b = m->getVertex(offset + vb);
		c = m->getVertex(offset + vc);
		d = m->getVertex(offset + vd);
		m->addRasterizedPrimitiveFace(d,c,b,a,material);
		// the top patches
		va = 1 + h*(v-1);
		vb = 1 +	i	 + h*(v-2);
		vc = 1 + (i+1)%h + h*(v-2);
		vd = 1 + (i+2)%h + h*(v-2);
		a = m->getVertex(offset + va);
		b = m->getVertex(offset + vb);
		c = m->getVertex(offset + vc);
		d = m->getVertex(offset + vd);
		m->addRasterizedPrimitiveFace(b,c,d,a,material);
	}
}
