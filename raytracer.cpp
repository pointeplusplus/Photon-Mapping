#include "raytracer.h"
#include "material.h"
#include "vectors.h"
#include "argparser.h"
#include "raytree.h"
#include "utils.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "photon_mapping.h"

#define eps 0.1

// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool RayTracer::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches) const {
	bool answer = false;

	// intersect each of the quads
	for (int i = 0; i < mesh->numOriginalQuads(); i++) {
		Face *f = mesh->getOriginalQuad(i);
		bool backfacing_hit = false;
		if (f->intersect(ray,h,args->intersect_backfacing, &backfacing_hit)) answer = true;
	}

	// intersect each of the primitives (either the patches, or the original primitives)
	if (use_rasterized_patches) {
		for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
			Face *f = mesh->getRasterizedPrimitiveFace(i);
			bool backfacing_hit = false;
			if (f->intersect(ray,h,args->intersect_backfacing, &backfacing_hit)) answer = true;
		}
	} else {
		int num_primitives = mesh->numPrimitives();
		for (int i = 0; i < num_primitives; i++) {
			if (mesh->getPrimitive(i)->intersect(ray,h)) answer = true;
		}
	}
	return answer;
}

// ===========================================================================
// does the recursive (shadow rays & recursive rays) work
//note: bounce_count is counting down, not up
Vec3f RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count) const {
//	std::cout << "In the trace ray" << std::endl
//	 << "    Bounce count " << bounce_count << std::endl;



	// First cast a ray and see if we hit anything.
	hit = Hit();
	bool intersect_with_objects = CastRay(ray,hit,false);
		
	// if there is no intersection, simply return the background color
	if (intersect_with_objects == false) {
		return Vec3f(srgb_to_linear(mesh->background_color.r()),
		 srgb_to_linear(mesh->background_color.g()),
		 srgb_to_linear(mesh->background_color.b()));
		
	}

	// otherwise decide what to do based on the material
	Material *m = hit.getMaterial();
	assert (m != NULL);

	// rays coming from the light source are set to white, don't bother to ray trace further.
	if (m->getEmittedColor().Length() > 0.001) {
		return Vec3f(1,1,1);
	} 
 
	
	Vec3f normal = hit.getNormal();
	Vec3f point = ray.pointAtParameter(hit.getT());
	Vec3f answer;

	// ----------------------------------------------
	//	start with the indirect light (ambient light)
	Vec3f diffuse_color = m->getDiffuseColor(hit.get_s(),hit.get_t());
	if (args->gather_indirect) {
		// photon mapping for more accurate indirect light
		answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + args->ambient_light);
	} else {
		// the usual ray tracing hack for indirect light
		answer = diffuse_color * args->ambient_light;
	}			

	// ----------------------------------------------
	// add contributions from each light that is not in shadow
	int num_lights = mesh->getLights().size();
	for (int i = 0; i < num_lights; i++) {

		Face *f = mesh->getLights()[i];
		Vec3f lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
		Vec3f myLightColor;
		Vec3f lightCentroid = f->computeCentroid();
		Vec3f dirToLightCentroid = lightCentroid-point;
		dirToLightCentroid.Normalize();

		double distToLightCentroid = (lightCentroid-point).Length();
		myLightColor = lightColor / (M_PI*distToLightCentroid*distToLightCentroid);

		// ===========================================
		// ADD SHADOW & SOFT SHADOW LOGIC
		// ===========================================

		std::vector<std::vector<Vec3f> > light_positions;
		for(unsigned int l = 0; l < mesh->getLights().size(); l++){
			std::vector<Vec3f> light;
			
			//if we only have 1 shadow sample, take the middle of the light
			if(args->num_shadow_samples == 1){
				light.push_back(mesh->getLights()[l]->computeCentroid());
			}
			//soft shadow logic -> take random points in the light
			else{
				
				//This is for regular soft shadow random sampling
				if(!args->stratified_shadows){
					for(int n = 0; n < args->num_shadow_samples; n++){
						light.push_back(mesh->getLights()[l]->RandomPoint());
					}
				}
				//This is for statified sampling
				else{
					double sqrt_num_shadow_samples = sqrt((double)args->num_shadow_samples);
					int stratified_num_shadow_samples = ceil(sqrt_num_shadow_samples);
	//				std::cout << "This is the ceiling of the sqrt of the num_shadow_samples: " << stratified_num_shadow_samples << std::endl;
					for(int x = 0; x < stratified_num_shadow_samples; x++){
						for(int y = 0; y < stratified_num_shadow_samples; y++){
							Vec3f point_on_light = mesh->getLights()[l]->RandomStratifiedPoint(x, y, stratified_num_shadow_samples);
							light.push_back(point_on_light);
						}
					}
				}
			}

			light_positions.push_back(light);
		}

		for(unsigned int l = 0; l < light_positions.size(); l++){
			for(unsigned int p = 0; p < light_positions[l].size(); p++){

				//find intersect with the point on the light
				//get direction to light point
				Vec3f direction_to_light = Vec3f(light_positions[l][p].x()-point.x(),light_positions[l][p].y()-point.y(), light_positions[l][p].z()-point.z());
				direction_to_light.Normalize();
				Ray ray_to_light = Ray(point,direction_to_light);
				Hit hit_to_light = Hit();
				//use direction to get intersect
				bool backfacing_hit = false;
				bool intersect_with_light = mesh->getLights()[l]->intersect(ray_to_light, hit_to_light, false, &backfacing_hit);
				assert(intersect_with_light);
				RayTree::AddShadowSegment(ray_to_light, 0, hit_to_light.getT());

				//check intersections with everything else
				Hit other_hits = Hit();
				CastRay(ray_to_light, other_hits, false);

				//if the closest thing was the lighting
				if(hit_to_light.getT() - other_hits.getT() < eps && other_hits.getT() - hit_to_light.getT() < eps){
					Vec3f light_from_ray = m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
					light_from_ray /= args->num_shadow_samples;
					answer += light_from_ray;
				}
			}
		}
		
		// add the lighting contribution from this particular light at this point
		// (fix this to check for blockers between the light & this surface)
		//answer += m->Shade(ray,hit,dirToLightCentroid,myLightColor,args);
	}
			
	// ----------------------------------------------
	// add contribution from reflection, if the surface is shiny
	Vec3f reflectiveColor = m->getReflectiveColor();
	
	// =================================
	// 	ADD REFLECTIVE LOGIC
	// =================================
	//check to make sure we aren't overriding the main segment
	if(bounce_count != 0 ){
		RayTree::AddReflectedSegment(ray, 0.0, hit.getT());
	}
	if(bounce_count > 0){
//		std::cout << "At least I made it into the loop" << std::endl;
		Vec3f current_location = ray.getOrigin() + (hit.getT() * ray.getDirection());
		Vec3f bounce_direction = ray.getDirection() - 2.0* (ray.getDirection().Dot3(hit.getNormal()) * hit.getNormal());

		Hit bounce_hit = Hit();
		Ray bounce_ray = Ray(current_location, bounce_direction);
		Vec3f bounce_color = TraceRay(bounce_ray ,bounce_hit, bounce_count-1);
		answer+= (bounce_color * reflectiveColor); 
	}

	return answer; 
}

