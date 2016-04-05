#include "face.h"
#include "matrix.h"
#include "utils.h"

#include <cmath>
#include <cassert>

// =========================================================================
// =========================================================================

double Face::getArea() const {
	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();
	return 
		AreaOfTriangle(DistanceBetweenTwoPoints(a,b),
									 DistanceBetweenTwoPoints(a,c),
									 DistanceBetweenTwoPoints(b,c)) +
		AreaOfTriangle(DistanceBetweenTwoPoints(c,d),
									 DistanceBetweenTwoPoints(a,d),
									 DistanceBetweenTwoPoints(a,c));
}

// =========================================================================

Vec3f Face::RandomPoint() const {
	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();

	double s = GLOBAL_mtrand.rand(); // random real in [0,1]
	double t = GLOBAL_mtrand.rand(); // random real in [0,1]

	Vec3f answer = s*t*a + s*(1-t)*b + (1-s)*t*d + (1-s)*(1-t)*c;
	return answer;
}

//Rebecca defined this
//takes in which iteration you want in the x and y direction, and how many samples there are in each direction
//Note: there will be (num_samples)^2 total stratified samples
Vec3f Face::RandomStratifiedPoint(int x_iter, int y_iter, int num_samples) const {
	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();

	double s = GLOBAL_mtrand.rand(); // random real in [0,1]
	double t = GLOBAL_mtrand.rand(); // random real in [0,1]

	double s_pos = (s + x_iter)/(double)num_samples;
	double t_pos = (t + y_iter)/(double)num_samples;

	Vec3f answer = s_pos*t_pos*a + s_pos*(1-t_pos)*b + (1-s_pos)*t_pos*d + (1-s_pos)*(1-t_pos)*c;
	return answer;
}

// =========================================================================
// the intersection routines

bool Face::intersect(const Ray &r, Hit &h, bool intersect_backfacing, bool* backfacing_hit) {
	// intersect with each of the subtriangles
	Vertex *a = (*this)[0];
	Vertex *b = (*this)[1];
	Vertex *c = (*this)[2];
	Vertex *d = (*this)[3];
	return triangle_intersect(r,h,a,b,c,intersect_backfacing, backfacing_hit) || triangle_intersect(r,h,a,c,d,intersect_backfacing, backfacing_hit);
}

bool Face::triangle_intersect(const Ray &r, Hit &h, Vertex *a, Vertex *b, Vertex *c, bool intersect_backfacing, bool* backfacing_hit) {

	*backfacing_hit = false;
	// compute the intersection with the plane of the triangle
	Hit h2 = Hit(h);
	if (!plane_intersect(r,h2,intersect_backfacing, backfacing_hit)) return 0;	

	// figure out the barycentric coordinates:
	Vec3f Ro = r.getOrigin();
	Vec3f Rd = r.getDirection();
	// [ ax-bx	 ax-cx	Rdx ][ beta	]		 [ ax-Rox ] 
	// [ ay-by	 ay-cy	Rdy ][ gamma ]	=	[ ay-Roy ] 
	// [ az-bz	 az-cz	Rdz ][ t		 ]		 [ az-Roz ] 
	// solve for beta, gamma, & t using Cramer's rule
	
	double detA = Matrix::det3x3(a->get().x()-b->get().x(),a->get().x()-c->get().x(),Rd.x(),
						 a->get().y()-b->get().y(),a->get().y()-c->get().y(),Rd.y(),
						 a->get().z()-b->get().z(),a->get().z()-c->get().z(),Rd.z());
	
	if (fabs(detA) <= 0.000001) return 0;
	assert (fabs(detA) >= 0.000001);

	double beta	= Matrix::det3x3(a->get().x()-Ro.x(),a->get().x()-c->get().x(),Rd.x(),
				a->get().y()-Ro.y(),a->get().y()-c->get().y(),Rd.y(),
				a->get().z()-Ro.z(),a->get().z()-c->get().z(),Rd.z()) / detA;
	
	double gamma = Matrix::det3x3(a->get().x()-b->get().x(),a->get().x()-Ro.x(),Rd.x(),
				a->get().y()-b->get().y(),a->get().y()-Ro.y(),Rd.y(),
				a->get().z()-b->get().z(),a->get().z()-Ro.z(),Rd.z()) / detA;

	//Case of an intersection
	if (beta >= -0.00001 && beta <= 1.00001 &&
			gamma >= -0.00001 && gamma <= 1.00001 &&
			beta + gamma <= 1.00001) {
		h = h2;
		// interpolate the texture coordinates
		double alpha = 1 - beta - gamma;
		double t_s = alpha * a->get_s() + beta * b->get_s() + gamma * c->get_s();
		double t_t = alpha * a->get_t() + beta * b->get_t() + gamma * c->get_t();
		h.setTextureCoords(t_s,t_t);
		assert (h.getT() >= EPSILON);
		return 1;
	}

	return 0;
}


bool Face::plane_intersect(const Ray &r, Hit &h, bool intersect_backfacing, bool* backfacing_hit) {

	// insert the explicit equation for the ray into the implicit equation of the plane

	// equation for a plane
	// ax + by + cz = d;
	// normal . p + direction = 0
	// plug in ray
	// origin + direction * t = p(t)
	// origin . normal + t * direction . normal = d;
	// t = d - origin.normal / direction.normal;

	Vec3f normal = computeNormal();
	double d = normal.Dot3((*this)[0]->get());

	double numer = d - r.getOrigin().Dot3(normal);
	double denom = r.getDirection().Dot3(normal);

	if (denom == 0) return 0;	// parallel to plane

	if (!intersect_backfacing && normal.Dot3(r.getDirection()) >= 0) 
		return 0; // hit the backside

	double t = numer / denom;
	if (t > EPSILON && t < h.getT()) {
		h.set(t,this->getMaterial(),normal,this);
		assert (h.getT() >= EPSILON);
		//hit the backside but that's okay in this case
		if (normal.Dot3(r.getDirection()) >= 0){
			*backfacing_hit = true;
		}
		return 1;
	}
	return 0;
}


inline Vec3f ComputeNormal(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3) {
	Vec3f v12 = p2;
	v12 -= p1;
	Vec3f v23 = p3;
	v23 -= p2;
	Vec3f normal;
	Vec3f::Cross3(normal,v12,v23);
	normal.Normalize();
	return normal;
}

bool isLinear(const Vec3f &p1, const Vec3f &p2, const Vec3f &p3){

	double triangle_area = AreaOfTriangle(DistanceBetweenTwoPoints(p1,p2),
											DistanceBetweenTwoPoints(p1, p3),
											DistanceBetweenTwoPoints(p2,p3));

	//Non-zero triangle area means that the points are not colinear
	if(triangle_area > EPSILON){
		return false;
	}
	return true;
}

void Face::computeCachedNormal()
{
	// note: this face might be non-planar, so average the two triangle normals
	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();

	//Checks to make sure we don't have three points that we are using here to check the normal be linear -- that could cause errors!! 
	if(isLinear(a,b,c)){
		cached_normal = new Vec3f(ComputeNormal(a,c,d));
		return;
	}
	if(isLinear(a,c,d)){
		cached_normal = new Vec3f(ComputeNormal(a,b,c));
		return;
	}
	cached_normal = new Vec3f(0.5 * (ComputeNormal(a,b,c) + ComputeNormal(a,c,d)));
	return;
}

Vec3f Face::computeNormal() {
	if(cached_normal == NULL)
	{
		computeCachedNormal();	
	}
	return *cached_normal;
}

bool Face::isConvex() const {

	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();

	// no D
	double t1 = AreaOfTriangle(DistanceBetweenTwoPoints(a,b),
									 DistanceBetweenTwoPoints(a,c),
									 DistanceBetweenTwoPoints(b,c));

	// no B
	double t2 = AreaOfTriangle(DistanceBetweenTwoPoints(a,c),
									 DistanceBetweenTwoPoints(a,d),
									 DistanceBetweenTwoPoints(c,d));

	// no C
	double t3 = AreaOfTriangle(DistanceBetweenTwoPoints(a,d),
									 DistanceBetweenTwoPoints(a,b),
									 DistanceBetweenTwoPoints(b,d));

	// no A
	double t4 = AreaOfTriangle(DistanceBetweenTwoPoints(c,d),
									 DistanceBetweenTwoPoints(b,d),
									 DistanceBetweenTwoPoints(b,c));
	

	

	if (!fabs((t1+t2)-(t3+t4) < .001)){
		std::cout << "        t1: " << t1 << std::endl;
		std::cout << "        t2: " << t2 << std::endl;
		std::cout << "        t3: " << t3 << std::endl;
		std::cout << "        t4: " << t4 << std::endl;
	}

	//assert(fabs((t1+t2)-(t3+t4) < .001));

	if (fabs((t1+t2)-(t3+t4) < .001)){
		return true;
	}
	else{
		return false;	
	} 
}

void Face::printVertices() const {

	Edge* edge_ptr = edge;
	Vertex* vert = edge_ptr->getEndVertex();
	std::cout << "        v1: " << vert->x() << ", " << vert->y() << ", " << vert->z() << std::endl;
	edge_ptr = edge_ptr->getNext();
	vert = edge_ptr->getEndVertex();
	std::cout << "        v2: " << vert->x() << ", " << vert->y() << ", " << vert->z() << std::endl;
	edge_ptr = edge_ptr->getNext();
	vert = edge_ptr->getEndVertex();
	std::cout << "        v3: " << vert->x() << ", " << vert->y() << ", " << vert->z() << std::endl;
	edge_ptr = edge_ptr->getNext();
	vert = edge_ptr->getEndVertex();
	std::cout << "        v4: " << vert->x() << ", " << vert->y() << ", " << vert->z() << std::endl;
}

double Face::shortestEdge() const {
	Vec3f a = (*this)[0]->get();
	Vec3f b = (*this)[1]->get();
	Vec3f c = (*this)[2]->get();
	Vec3f d = (*this)[3]->get();

	double e1 = DistanceBetweenTwoPoints(a,b);
	double e2 = DistanceBetweenTwoPoints(b,c);
	double e3 = DistanceBetweenTwoPoints(c,d);
	double e4 = DistanceBetweenTwoPoints(d,a);

	double smallest = e1;
	if(e2 < smallest) smallest = e2;
	if(e3 < smallest) smallest = e3;
	if(e4 < smallest) smallest = e4;

	return smallest;

}
