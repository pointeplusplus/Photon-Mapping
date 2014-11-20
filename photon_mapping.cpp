#include "glCanvas.h"

#include <iostream>
#include <algorithm>
#include <cmath>

#include "argparser.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "matrix.h"
#include "face.h"
#include "primitive.h"
#include "raytree.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"

#define PHOTON_VISUALIZATION_ALPHA 0.7
//#define DRAW_VISUALIZATION false
#define DRAW_PHOTON_PATHS false	
#define DRAW_COLORED_NORMALS true
#define DRAW_ESCAPING_PHOTONS false
//#define ORIGINAL_N_VAL 1.0  // TODO, this should be passed in because it won't always be coming from air
#define REFRACTIVE_INDEX_OF_AIR 1.000293
#define NORMAL_VISUALIZATION_LENGTH .3

Vec3f wavelengthToRGB(double wavelength);

// ==========
// DESTRUCTOR
PhotonMapping::~PhotonMapping() {
	// cleanup all the photons
	delete kdtree;
}

// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const Vec3f &position, const Vec3f &direction, 
				const float wavelength, int iter, Vec4f viz_color, float current_n_val) {

	//std::cout << "In trace photon" << std::endl;

	//FOR DEBUG ONLY
	if(viz_color[0] == 1.0 && viz_color[1] == 0.5 && viz_color[2] == 0.0){
		//std::cout << "In trace photon with refraction color" << std::endl;	
	}


	// TODO: something else with this vector other than making it again every iteration
	std::vector<Vec4f> colors;
	
	colors.push_back(Vec4f(1.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // magenta
	//colors.push_back(Vec4f(1.0, 0.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // red
	colors.push_back(Vec4f(0.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // blue
	colors.push_back(Vec4f(0.0, 1.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // green
	colors.push_back(Vec4f(1.0, 1.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // yellow
	colors.push_back(Vec4f(0.0, 1.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // cyan
	
	//colors.push_back(Vec4f(1.0, 1.0, 1.0, PHOTON_VISUALIZATION_ALPHA));

	//Add this photon to thee kd tree
	Photon* this_iteration_photon = new Photon(position, direction, wavelength, iter);

	// kdtree is false during the t-press situation
	if(kdtree){
		kdtree->AddPhoton(*this_iteration_photon);
	}
	
	//If we can bounce no more, return
	if (iter >= args->num_bounces){
		return;
	}

	Ray ray = Ray(position, direction);
	Hit hit;
	//find the next thing or the photon to bounce off of
	if(CastRay(ray, hit, false)){
		Vec3f hit_normal = hit.getNormal();
//		std::cout << "    cast ray == true" << std::endl;
		if(hit.getIsBackfacing()){
			hit_normal = -1 * hit_normal;
//			std:: cout << "    backfacing hit" << std::endl;
		}
		/*
		//FOR DEBUG ONLY
		if(viz_color[0] == 1.0 && viz_color[1] == 0.5 && viz_color[2] == 0.0){
			std::cout << "Cast ray worked with refraction color.  Iter: " << iter << std::endl;	
		}
		*/
		Material* material = hit.getMaterial();
		Vec3f bounce_location = ray.getOrigin() + (hit.getT() * ray.getDirection());

		//If it is not refracted, set the viz color.  Otherwise, keep orange.
		if(!(viz_color[0] == 1.0 && viz_color[1] == 0.5 && viz_color[2] == 0.0)){
			viz_color = colors[iter%colors.size()];
		}
		
		//Only add the photon to the visualization if it isn't the first bounce (iter is already incremented)
		// checking kdtree determines whether or not this was called from pressing t.  If so, draw main segment also.
		if(iter != 0 || !kdtree){
			//visualization_line_segments.push_back(LineSegment(position, bounce_location, viz_color));	
			Vec3f photon_color = wavelengthToRGB(wavelength);
			Vec4f photon_color_with_alpha = Vec4f(photon_color[0], photon_color[1], photon_color[2], PHOTON_VISUALIZATION_ALPHA);

			if(DRAW_PHOTON_PATHS){
				RayTree::AddGeneralSegment(ray,0,hit.getT(), photon_color_with_alpha);
			}
			//Normal of the hit (white and fully opaque) -- make the time constant just to draw a line
			Ray hit_normal_ray = Ray(bounce_location, hit_normal);
			if(DRAW_COLORED_NORMALS){
				RayTree::AddGeneralSegment(hit_normal_ray, 0, NORMAL_VISUALIZATION_LENGTH, photon_color_with_alpha);
			}
			else{
				RayTree::AddGeneralSegment(hit_normal_ray, 0, NORMAL_VISUALIZATION_LENGTH, Vec4f(1.0, 1.0, 1.0, 1.0));
			}
			
		}

		//Change the color again to catch the refractive case (so it doesn't continue to be orange)
		viz_color = colors[iter%colors.size()];
		iter++;


		//If we are reflective
		//if((material->getReflectiveColor())[0] != 0.0 || 
		//	(material->getReflectiveColor())[1] != 0.0 ||
		//	(material->getReflectiveColor())[2] != 0.0){
		if(true){ //TODO: when to reflect
			Vec3f bounce_direction = ray.getDirection() - 2.0* (ray.getDirection().Dot3(hit_normal) * hit_normal);
			//This was for energy and RGB
			/*
			Vec3f bounce_reflective_color = Vec3f(energy[0]* (material->getReflectiveColor()[0]),
													energy[1]* (material->getReflectiveColor()[1]),
													energy[2]* (material->getReflectiveColor()[2]));
			*/
			//TODO: double check that reflection always leads to the photon going into air
			TracePhoton(bounce_location,bounce_direction,wavelength,iter, viz_color, current_n_val);
		}
		/*
		//Diffuse
		//TODO: maybe remove the raytree segment
		if((material->getDiffuseColor())[0] != 0.0 || 
			(material->getDiffuseColor())[1] != 0.0 ||
			(material->getReflectiveColor())[2] != 0.0){

			Vec3f bounce_direction = Vec3f(2*GLOBAL_mtrand.rand()-1, 2*GLOBAL_mtrand.rand()-1, 2*GLOBAL_mtrand.rand()-1);
			Vec3f bounce_diffuse_color = Vec3f(energy[0]* (material->getDiffuseColor()[0]),
												energy[1]* (material->getDiffuseColor()[1]),
												energy[2]* (material->getDiffuseColor()[2]));

			bounce_direction.Normalize();
			TracePhoton(bounce_location,bounce_direction,bounce_diffuse_color,iter, viz_color);
		}	
		*/
		// REFRACTION
		//TODO: when to refract
		Vec3f incoming_direction = ray.getDirection();
		float next_n_val = material->getRefractiveIndex();
		if(hit.getIsBackfacing()){
			next_n_val = REFRACTIVE_INDEX_OF_AIR;
		}
		double n = current_n_val / next_n_val;
		double cosI = -1 * hit_normal.Dot3(incoming_direction);
		double sinT2 = n * n * (1.0 - cosI * cosI);
		if(sinT2 <= 1.0){
			//std::cout << "    Should also refract" << std::endl;
			double cosT = sqrt(1.0 - sinT2);
			Vec3f outgoing_direction = n * incoming_direction + (n * cosI - cosT) * hit_normal;
			//TODO: if we are leaving the material, send the refractive index of air instead of the material's refractive index
			//	to do this, check angle of the normal and the hit										
			TracePhoton(bounce_location,outgoing_direction,wavelength,iter, Vec4f(1.0, 0.5, 0.0, PHOTON_VISUALIZATION_ALPHA), next_n_val);
		}
		else{
			//TODO: total internal refraction
			std::cout << "    Hitting the TIR case" << std::endl;
		}

		//Old way
		/*
		incoming_direction.Negate(); //Negate such that the angle between them is proper
		float incoming_angle = incoming_direction.AngleBetween(hit.getNormal());
		float outgoing_angle = asin(sin(incoming_angle) * current_n_val / material->getRefractiveIndex());

		float r = current_n_val / material->getRefractiveIndex();

		Vec3f find_c =  hit.getNormal();
		find_c.Negate();
		float c = find_c.Dot3(ray.getDirection());`
		
		//Check for total internal refraction

		float radical = 1 - ( pow(r,2) * (1 - pow(c,2) ) );
		if(radical >= 0){
			Vec3f outgoing_direction =  ( r * ray.getDirection() ) +
										( sqrt(radical) * hit.getNormal() ) ;

			//TODO: if we are leaving the material, send the refractive index of air instead of the material's refractive index
			//	to do this, check angle of the normal and the hit										
			TracePhoton(bounce_location,outgoing_direction,wavelength,iter, Vec4f(1.0, 0.5, 0.0, PHOTON_VISUALIZATION_ALPHA), material->getRefractiveIndex());
			//std::cout << "This is supposed to be a refractive thing" << std::endl;

			//Debug segment just to see what's up
			//visualization_line_segments.push_back(LineSegment(bounce_location, bounce_location + (outgoing_direction * 5), Vec4f(1.0, 1.0, 0.0, PHOTON_VISUALIZATION_ALPHA)));

		}
		else{
			// TODO: do something if total internal refraction
		}
		*/
		

		//std::cout << "    number of segments: "  << visualization_line_segments.size() << std::endl;
	}
	
	//Visualize photons even when they don't hit anything.
	else if(DRAW_ESCAPING_PHOTONS){
		Vec3f bounce_direction = ray.getDirection() - 2.0* (ray.getDirection().Dot3(hit.getNormal()) * hit.getNormal());
		//Ray bounce = Ray(bounce_location, bounce);
		RayTree::AddGeneralSegment(ray, 0, NORMAL_VISUALIZATION_LENGTH * 5, Vec4f(0.0,0.0,0.0,1.0));
	}
}

// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {
	std::cout << "trace photons" << std::endl;

	RayTree::Activate();

	// first, throw away any existing photons
	delete kdtree;

	// consruct a kdtree to store the photons
	BoundingBox *bb = mesh->getBoundingBox();
	Vec3f min = bb->getMin();
	Vec3f max = bb->getMax();
	Vec3f diff = max-min;
	min -= 0.001*diff;
	max += 0.001*diff;
	kdtree = new KDTree(BoundingBox(min,max));

	// photons emanate from the light sources
	const std::vector<Face*>& lights = mesh->getLights();

	// compute the total area of the lights
	double total_lights_area = 0;
	for (unsigned int i = 0; i < lights.size(); i++) {
		total_lights_area += lights[i]->getArea();
	}

	// shoot a constant number of photons per unit area of light source
	// (alternatively, this could be based on the total energy of each light)
	for (unsigned int i = 0; i < lights.size(); i++) {	
		double my_area = lights[i]->getArea();
		int num = args->num_photons_to_shoot * my_area / total_lights_area;
		// the initial energy for this photon 
		//Vec3f energy = my_area/double(num) * lights[i]->getMaterial()->getEmittedColor();
		//replace energy with photon color 
		Vec3f normal = lights[i]->computeNormal();
		for (int j = 0; j < num; j++) {
			Vec3f start = lights[i]->RandomPoint();
			// the initial direction for this photon (for diffuse light sources)
			Vec3f direction = RandomDiffuseDirection(normal);
			//Vec4f photon_color = Vec4f(1.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA);
			Vec4f photon_color = Vec4f(1.0, 1.0, 1.0, PHOTON_VISUALIZATION_ALPHA);
			float wavelength = (GLOBAL_mtrand.rand() * 400) + 380; //random number between 380 and 780 (visible light)
			TracePhoton(start,direction,wavelength,0, photon_color, REFRACTIVE_INDEX_OF_AIR);
		}
	}

	RayTree::Deactivate();
}

// ======================================================================
// During ray tracing, when a diffuse (or partially diffuse) object is
// hit, gather the nearby photons to approximate indirect illumination

Vec3f PhotonMapping::GatherIndirect(const Vec3f &point, const Vec3f &normal, const Vec3f &direction_from) const {


	if (kdtree == NULL) { 
		std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
		return Vec3f(0,0,0); 
	}


	// ================================================================
	// ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
	// ================================================================


	// collect the closest args->num_photons_to_collect photons
	// determine the radius that was necessary to collect that many photons
	// average the energy of those photons over that radius


	// return the color
	return Vec3f(0,0,0);
}



// ======================================================================
// PHOTON VISUALIZATION FOR DEBUGGING
// ======================================================================


void PhotonMapping::initializeVBOs() {
	glGenBuffers(1, &photon_verts_VBO);
	glGenBuffers(1, &photon_direction_indices_VBO);
	glGenBuffers(1, &kdtree_verts_VBO);
	glGenBuffers(1, &kdtree_edge_indices_VBO);

	//Visualization
	glGenBuffers(1, &raytree_verts_VBO);
 	glGenBuffers(1, &raytree_edge_indices_VBO);
}

void PhotonMapping::setupVBOs() {
	photon_verts.clear();
	photon_direction_indices.clear();
	kdtree_verts.clear();
	kdtree_edge_indices.clear();

	raytree_verts.clear();
 	raytree_edge_indices.clear();

	// initialize the data
	int dir_count = 0;
	int edge_count = 0;
	BoundingBox *bb = mesh->getBoundingBox();
	double max_dim = bb->maxDim();

	if (kdtree == NULL) return;
	std::vector<const KDTree*> todo;	
	todo.push_back(kdtree);
	while (!todo.empty()) {
		const KDTree *node = todo.back();
		todo.pop_back(); 
		if (node->isLeaf()) {
			const std::vector<Photon> &photons = node->getPhotons();
			int num_photons = photons.size();
			for (int i = 0; i < num_photons; i++) {
	const Photon &p = photons[i];
	//Vec3f energy = p.getEnergy()*args->num_photons_to_shoot;
	//Changed for wavelength TODO: ask Barb about this
	Vec3f energy = wavelengthToRGB(p.getWavelength()) * args->num_photons_to_shoot;
	const Vec3f &position = p.getPosition();
	//Vec3f other = position + p.getDirectionFrom()*0.02*max_dim;
	Vec3f other = position + p.getDirectionFrom()*0.005*max_dim;
	photon_verts.push_back(VBOPosColor(position,energy));
	photon_verts.push_back(VBOPosColor(other,energy));
	photon_direction_indices.push_back(VBOIndexedEdge(dir_count,dir_count+1)); dir_count+=2;
			}

			// initialize kdtree vbo
			const Vec3f& min = node->getMin();
			const Vec3f& max = node->getMax();
			kdtree_verts.push_back(VBOPos(Vec3f(min.x(),min.y(),min.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(min.x(),min.y(),max.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(min.x(),max.y(),min.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(min.x(),max.y(),max.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(max.x(),min.y(),min.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(max.x(),min.y(),max.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(max.x(),max.y(),min.z())));
			kdtree_verts.push_back(VBOPos(Vec3f(max.x(),max.y(),max.z())));

			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count ,edge_count+1)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+1,edge_count+3)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+3,edge_count+2)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+2,edge_count	)); 

			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+4,edge_count+5)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+5,edge_count+7)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+7,edge_count+6)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+6,edge_count+4)); 

			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count ,edge_count+4)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+1,edge_count+5)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+2,edge_count+6)); 
			kdtree_edge_indices.push_back(VBOIndexedEdge(edge_count+3,edge_count+7)); 


			edge_count += 8;

		} else {
			todo.push_back(node->getChild1());
			todo.push_back(node->getChild2());
		} 
	}
	assert (2*photon_direction_indices.size() == photon_verts.size());
	int num_directions = photon_direction_indices.size();
	int num_edges = kdtree_edge_indices.size();
/*
	//visualization
	unsigned int i;
	int count = 0;
	for (i = 0; i < visualization_line_segments.size(); i++) {
		raytree_verts.push_back(VBOPosColor4(visualization_line_segments[i].getStart(),visualization_line_segments[i].getColor()));
		raytree_verts.push_back(VBOPosColor4(visualization_line_segments[i].getEnd(),visualization_line_segments[i].getColor()));
		raytree_edge_indices.push_back(VBOIndexedEdge(count,count+1)); count+=2;
	}
	*/
	assert (2*raytree_edge_indices.size() == raytree_verts.size());
	int num_viz_edges = raytree_edge_indices.size();

	// cleanup old buffer data (if any)
	cleanupVBOs();

	// copy the data to each VBO
	if (num_directions > 0) {
		glBindBuffer(GL_ARRAY_BUFFER,photon_verts_VBO); 
		glBufferData(GL_ARRAY_BUFFER,
		 sizeof(VBOPosColor) * num_directions * 2,
		 &photon_verts[0],
		 GL_STATIC_DRAW); 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,photon_direction_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * num_directions,
		 &photon_direction_indices[0], GL_STATIC_DRAW);
	} 

	if (num_edges > 0) {
		glBindBuffer(GL_ARRAY_BUFFER,kdtree_verts_VBO); 
		glBufferData(GL_ARRAY_BUFFER,
		 sizeof(VBOPos) * num_edges * 2,
		 &kdtree_verts[0],
		 GL_STATIC_DRAW); 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,kdtree_edge_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * num_edges,
		 &kdtree_edge_indices[0], GL_STATIC_DRAW);
	} 
	// copy the data to each VBO
	if (num_viz_edges > 0) {
		glBindBuffer(GL_ARRAY_BUFFER,raytree_verts_VBO); 
		glBufferData(GL_ARRAY_BUFFER,
			sizeof(VBOPosColor4) * num_viz_edges * 2,
			&raytree_verts[0],
			GL_STATIC_DRAW); 
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,raytree_edge_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
			sizeof(VBOIndexedEdge) * num_viz_edges,
			&raytree_edge_indices[0], GL_STATIC_DRAW);
	} 

}

void PhotonMapping::drawVBOs() {

	glDisable(GL_LIGHTING);
	if (args->render_photons) {
		int num_directions = photon_direction_indices.size();
		if (num_directions > 0) {
			// render directions
			glLineWidth(1);
			glBindBuffer(GL_ARRAY_BUFFER, photon_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPosColor), BUFFER_OFFSET(0));
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor), BUFFER_OFFSET(12));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, photon_direction_indices_VBO);
			glDrawElements(GL_LINES, num_directions*2, GL_UNSIGNED_INT, BUFFER_OFFSET(0));
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_VERTEX_ARRAY);
			// render hit points
			glPointSize(3);
			glBindBuffer(GL_ARRAY_BUFFER, photon_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPosColor)*2, BUFFER_OFFSET(0));
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(3, GL_FLOAT, sizeof(VBOPosColor)*2, BUFFER_OFFSET(12));
			glDrawArrays(GL_POINTS, 0, num_directions);
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
	}
	if (args->render_kdtree) {
		int num_edges = kdtree_edge_indices.size();
		if (num_edges > 0) {
			glDisable(GL_LIGHTING);
			// render edges
			glLineWidth(1);
			glColor3f(0,1,1);
			glBindBuffer(GL_ARRAY_BUFFER, kdtree_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPos), BUFFER_OFFSET(0));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, kdtree_edge_indices_VBO);
			glDrawElements(GL_LINES, num_edges*2, GL_UNSIGNED_INT, BUFFER_OFFSET(0));
			glDisableClientState(GL_VERTEX_ARRAY);
		}
	}

	/*
	if(DRAW_VISUALIZATION){
		int num_viz_edges = raytree_edge_indices.size();
		if (num_viz_edges > 0) {

		// this allows you to see rays passing through objects
		// turn off the depth test and blend with the current pixel color
			glDisable(GL_DEPTH_TEST);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

			glDisable(GL_LIGHTING);
			glLineWidth(2);
			glBindBuffer(GL_ARRAY_BUFFER, raytree_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPosColor4), BUFFER_OFFSET(0));
			glEnableClientState(GL_COLOR_ARRAY);
			glColorPointer(4, GL_FLOAT, sizeof(VBOPosColor4), BUFFER_OFFSET(12));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, raytree_edge_indices_VBO);
			glDrawElements(GL_LINES, num_viz_edges*2, GL_UNSIGNED_INT, BUFFER_OFFSET(0));
			glDisableClientState(GL_COLOR_ARRAY);
			glDisableClientState(GL_VERTEX_ARRAY);
			glEnable(GL_DEPTH_TEST);
		}
	}
	*/
}

void PhotonMapping::cleanupVBOs() {
	glDeleteBuffers(1, &photon_verts_VBO);
	glDeleteBuffers(1, &photon_direction_indices_VBO);
	glDeleteBuffers(1, &kdtree_verts_VBO);
	glDeleteBuffers(1, &kdtree_edge_indices_VBO);

	//Visualization
	glDeleteBuffers(1, &raytree_verts_VBO);
	glDeleteBuffers(1, &raytree_edge_indices_VBO);
}


// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool PhotonMapping::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const {
	bool answer = false;

	// intersect each of the quads
	for (int i = 0; i < mesh->numOriginalQuads(); i++) {
		Face *f = mesh->getOriginalQuad(i);
		bool backfacing_hit = false;
		if (f->intersect(ray,h,args->intersect_backfacing, &backfacing_hit)){
			answer = true;
			h.setIsBackfacing(backfacing_hit);
		}
	}

	// intersect each of the primitives (either the patches, or the original primitives)
	if (use_rasterized_patches) {
		for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
			Face *f = mesh->getRasterizedPrimitiveFace(i);
			bool backfacing_hit = false;
			if (f->intersect(ray,h,args->intersect_backfacing, &backfacing_hit)){
				answer = true;
				h.setIsBackfacing(backfacing_hit);
			}
		}
	} else {
		int num_primitives = mesh->numPrimitives();
		for (int i = 0; i < num_primitives; i++) {
			if (mesh->getPrimitive(i)->intersect(ray,h)){
				answer = true;	
				//TODO:fix for primitives
				h.setIsBackfacing(false);
			}
		}
	}
	return answer;
}