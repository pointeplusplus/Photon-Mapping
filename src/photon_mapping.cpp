#include "glCanvas.h"

#include <iostream>
#include <algorithm>
#include <cmath>
#include <ctime>
#include <vector>
#include <fstream>
#include <chrono>

#include "argparser.h"
#include "photon_mapping.h"
#include "photon.h"
#include "boundingbox.h"
#include "mesh.h"
#include "matrix.h"
#include "face.h"
#include "primitive.h"
#include "raytree.h"
#include "ray.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"

#define PHOTON_VISUALIZATION_ALPHA 0.7
//#define DRAW_VISUALIZATION false
#define DIRECTED_LIGHT false
#define DRAW_PHOTON_PATHS false
#define DRAW_FIRST_BOUNCE false
#define DRAW_COLORED_PATHS false
#define DRAW_NORMALS false
#define DRAW_COLORED_NORMALS false
#define DRAW_ESCAPING_PHOTONS false
#define NUM_BOUNCE_VIZ false
//#define ORIGINAL_N_VAL 1.0  // TODO, this should be passed in because it won't always be coming from air
#define REFRACTIVE_INDEX_OF_AIR 1.000293
#define NORMAL_VISUALIZATION_LENGTH .1
#define PI 3.14159265
#define eps 0.001	
int num_shot = 0;	

Vec3f wavelengthToRGB(double wavelength);
Vec3f mixColors(const std::vector<std::pair<Photon, float> > & wavelengths);

//TODO: change wavelength?
float nValueSellmeier(float wavelength, std::vector<float> B, std::vector<float> C){

	//std::cout << "In Sellmeier" << std::endl;
	//std::cout << "    B and C length: " << B.size() << " " << C.size() << std::endl;
	// If we don't have any Sellmeier values, return the refractive index of air
	if(B.size() < 1 || C.size() < 1){
		return REFRACTIVE_INDEX_OF_AIR;
	}
	//value from -0.5 to .5 based on wavelength, centered around 580, which is the middle wavelength
	//float wavelength_factor = (wavelength - 580) / 400.0;
	//return n_val + wavelength_factor;
	float frac = 0; 
	for(unsigned int val = 0; val < B.size() && val < C.size(); val++){
		frac += B[val] * pow(wavelength, 2) / (pow(wavelength, 2) - C[val]);	
	}
	
	frac += 1;
	assert (frac > 0);
	//Test prints:
	
	//std::cout << "    Wavelength: " << wavelength << std::endl;
	//std::cout << "    N value: " << sqrt(frac) << std::endl;
	return sqrt(frac);
}

// ==========
// DESTRUCTOR
PhotonMapping::~PhotonMapping() {
	// cleanup all the photons
	delete kdtree;
}

void PhotonMapping::printEscapingFacePhoton(){
	/*
	//After all of the photons have been traced, print the light escaping each face
	std::cout << "Number of rays leaving faces:" << std::endl;
	Face* face;
	for (int f = 0; f < mesh->numFaces(); f++){
		face = mesh->getFace(f);
		std::cout << "    " << f << ": " << face->getNumRaysLeavingFace() << "   ";
	}
	std::cout << "filename: " << args->output_file << std::endl; 
	

	

	int total_interior_bounces = 0;
	int total_rays_reflected = 0; 
	int total_rays_entering = 0;
	int total_rays_leaving = 0;	
	Face* face;

	for (int f = 0; f < mesh->numFaces(); f++){

		face = mesh->getFace(f);
		//Summing 
		if(face->getMaterial()->getName() == "diamond"){
			total_interior_bounces += face->getNumInteriorBounces();
			total_rays_reflected += face->getNumRaysReflected(); 
			total_rays_entering += face->getNumRaysEnteringFace();
			total_rays_leaving += face->getNumRaysLeavingFace();
		}
	}
	std::cout << total_interior_bounces << " " << total_rays_reflected << " " << total_rays_entering << " " << total_rays_leaving << std::endl;	
	*/
}

void PhotonMapping::printOutputFile(){
	
	
	//summming variables
	int total_interior_bounces = 0;
	int total_rays_reflected = 0; 
	int total_rays_entering = 0;
	int total_rays_leaving = 0;	

	std::ofstream output(args->output_file.c_str());

	std::cout << "outputfile " << args->output_file.c_str() << std::endl;
	output << "Photons shot: " << args->num_photons_to_shoot << std::endl;
	output << "Input file: " << args->input_file << std::endl;
	output << "Number of bounces: " << args->num_bounces << std::endl;
	output << "Number of shadow samples: " << args->num_shadow_samples << std::endl;
	Face* face;
	output << "Face_Number Face_Material Normal Area Internal_Bounce Reflected_Rays Rays_Entering_Face Rays_Leaving_Face" << std::endl;
	for (int f = 0; f < mesh->numFaces(); f++){

		face = mesh->getFace(f);

		//Summing 
		if(face->getMaterial()->getName() == "diamond"){
			total_interior_bounces += face->getNumInteriorBounces();
			total_rays_reflected += face->getNumRaysReflected(); 
			total_rays_entering += face->getNumRaysEnteringFace();
			total_rays_leaving += face->getNumRaysLeavingFace();
		}
			
		//Printing 
		Vec3f normal = face->computeNormal();
		output << f << ": " 
			<< face->getMaterial()->getName() << " "
			<< "{" << normal.x() << "," << normal.y() << "," << normal.z() << "} "
			<< face->getArea() << " "
			<< face->getNumInteriorBounces() << " "
			<< face->getNumRaysReflected() << " "
			<< face->getNumRaysEnteringFace() << " "
			<< face->getNumRaysLeavingFace() << std::endl;
	}
	//print out totals
	output << total_interior_bounces << " " << total_rays_reflected << " " << total_rays_entering << " " << total_rays_leaving << std::endl;

	//print out all directions the light is leaving all faces together
	output << "leaving_directions: ";
	for(int f = 0; f < mesh->numFaces(); f++){
		face = mesh->getFace(f);
		std::vector<Vec3f> directions = face->getLeavingDirections();
		//std::cout << "size of directions vector: " << directions.size() << std::endl;
		for(int d = 0; d < directions.size(); d++){
			output << "{" << directions[d].x() << "," << directions[d].y() << "," << directions[d].z() << "} ";
		}
	}
	output << std::endl;
	std::cout << "Output file printed successfully" << std::endl;
	
}

// ========================================================================
// Recursively trace a single photon

bool PhotonMapping::TracePhoton(const Vec3f &position, const Vec3f &direction, 
				const float wavelength, int iter, Vec4f viz_color, Material* current_material, float current_n_val, bool single_photon) {

	//std::cout << "Iter for this photon: " << iter << std::endl;

	// NOTE: current_material can == NULL 
	// TODO: something else with this vector other than making it again every iteration
	std::vector<Vec4f> colors;

	colors.push_back(Vec4f(1.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // magenta
	//colors.push_back(Vec4f(1.0, 0.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // red
	colors.push_back(Vec4f(0.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // blue
	colors.push_back(Vec4f(0.0, 1.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // green
	colors.push_back(Vec4f(1.0, 1.0, 0.0, PHOTON_VISUALIZATION_ALPHA)); // yellow
	colors.push_back(Vec4f(0.0, 1.0, 1.0, PHOTON_VISUALIZATION_ALPHA)); // cyan
	colors.push_back(Vec4f(1.0, 0.5, 0.0, PHOTON_VISUALIZATION_ALPHA)); // orange

	//Add this photon to thee kd tree
	//Photon p(position, direction, wavelength, iter);
	//Photon* this_iteration_photon = &p;
	
	// Push back onto photon vector for KDTree creating later
	// -- Commented out for runs where we don't need visualization --	
	//photon_lock.lock();
	//photons.push_back(*this_iteration_photon);
	//photon_lock.unlock();



	//If we can bounce no more, return
	// Check to see if the material absorbs the light 
	// TODO:  allow there to be a more complicated absorbing pattern that just a straight percentage
	float absorb_random = GLOBAL_mtrand.rand();
	if (iter >= args->num_bounces || (current_material != NULL && absorb_random < current_material->getAbsorbed() ) ){
		//If this is the only photon, print about face data
		if(single_photon){
			printEscapingFacePhoton();
			printOutputFile();
		}
		return true;
	}

	Ray ray = Ray(position, direction);
	
	Hit hit;
	

	//find the next thing or the photon to bounce off of
	if(CastRay(ray, hit, false)){
		Vec3f hit_normal = hit.getNormal();

		Material* material = hit.getMaterial();
		Vec3f bounce_location = ray.getOrigin() + (hit.getT() * ray.getDirection());

		//If it is not refracted, set the viz color.  Otherwise, keep orange.
		//if(!(viz_color[0] == 1.0 && viz_color[1] == 0.5 && viz_color[2] == 0.0)){
			viz_color = colors[iter%colors.size()];
		//}
		
		//Only add the photon to the visualization if it isn't the first bounce (iter is already incremented)
		// checking kdtree determines whether or not this was called from pressing t.  If so, draw main segment also.
		if((iter != 0 || DRAW_FIRST_BOUNCE) || single_photon){
			//visualization_line_segments.push_back(LineSegment(position, bounce_location, viz_color));	
			Vec3f photon_color = wavelengthToRGB(wavelength);
 			Vec4f photon_color_with_alpha = Vec4f(photon_color[0], photon_color[1], photon_color[2], PHOTON_VISUALIZATION_ALPHA);

			if(DRAW_PHOTON_PATHS){
				ray_tree_lock.lock();
				if(NUM_BOUNCE_VIZ){
					RayTree::AddGeneralSegment(ray,0,hit.getT(), viz_color);
				}
				else if (DRAW_COLORED_PATHS){
					RayTree::AddGeneralSegment(ray,0,hit.getT(), photon_color_with_alpha);
				}
				else{
					RayTree::AddGeneralSegment(ray,0,hit.getT(), Vec4f(1.0,1.0,1.0,PHOTON_VISUALIZATION_ALPHA));
				}
				ray_tree_lock.unlock();
			}
			//Normal of the hit (white and fully opaque) -- make the time constant just to draw a line
			Ray hit_normal_ray = Ray(bounce_location, hit_normal);
			ray_tree_lock.lock();
			if(DRAW_NORMALS){
				if(DRAW_COLORED_NORMALS){
					RayTree::AddGeneralSegment(hit_normal_ray, 0, NORMAL_VISUALIZATION_LENGTH, photon_color_with_alpha);
				}
				else{
					RayTree::AddGeneralSegment(hit_normal_ray, 0, NORMAL_VISUALIZATION_LENGTH, Vec4f(1.0, 1.0, 1.0, 1.0));
				}
			}
			ray_tree_lock.unlock();
		}

		//Change the color again to catch the refractive case (so it doesn't continue to be orange)
		viz_color = colors[iter%colors.size()];
		iter++;

		if(hit.getIsBackfacing()){
			hit_normal = -1 * hit_normal;
		}

		//Reflectence for s and p polarized light respectively
		float rs = 0;
		float rp = 0;
		Vec3f refraction_direction = Vec3f(0,0,0); // Will be assigned later
		//Refraction calculations for Frensel equations
		Vec3f incoming_direction = ray.getDirection();
		float incoming_angle_with_surface = hit_normal.AngleBetweenRadians(-1 * incoming_direction);

		//reassigned when the photon either escapes or doesn't -- this is used for material 1 and material 2
		float next_n_val = nValueSellmeier(wavelength, material->getRefractiveB(), material->getRefractiveC());
		//float next_n_val = material->getRefractiveIndex();
		if(hit.getIsBackfacing()){
			//next_n_val = nValueSellmeier(wavelength, AIR_B, AIR_C);
			next_n_val = REFRACTIVE_INDEX_OF_AIR;
		}
		double n = current_n_val / next_n_val;
		double cosI = -1 * hit_normal.Dot3(incoming_direction);
		double sinT2 = n * n * (1.0 - cosI * cosI);
		if(sinT2 <= 1.0){
			double cosT = sqrt(1.0 - sinT2);
			refraction_direction = n * incoming_direction + (n * cosI - cosT) * hit_normal;
			float outgoing_angle_with_surface = refraction_direction.AngleBetweenRadians(-1 * hit_normal);
			//TODO: if we are leaving the material, send the refractive index of air instead of the material's refractive index
			//	to do this, check angle of the normal and the hit
			/*
			if(hit.getMaterial()->getName() == "diamond"){
				std::cout << "    Incoming angle: " << incoming_angle_with_surface << std::endl;
				std::cout << "    Outgoing angle: " << outgoing_angle_with_surface << std::endl;
				std::cout << "    Current N val: " << current_n_val << std::endl;
				std::cout << "    Next N val: "  << next_n_val << std::endl;
			}
			*/
			float rs_numerator = (current_n_val * cos(incoming_angle_with_surface)) - (next_n_val * cos(outgoing_angle_with_surface));
			float rs_denominator = (current_n_val * cos(incoming_angle_with_surface)) + (next_n_val * cos(outgoing_angle_with_surface));
			rs = pow(rs_numerator/rs_denominator, 2);
			
			float rp_numerator = (current_n_val * cos(outgoing_angle_with_surface)) - (next_n_val * cos(incoming_angle_with_surface));
			float rp_denominator = (current_n_val * cos(outgoing_angle_with_surface)) + (next_n_val * cos(incoming_angle_with_surface));
			rp = pow(rp_numerator/rp_denominator, 2);
		}
		//With total internal reflection, reflectance values are 1
		else{
			rs = 1;
			rp = 1;
		}

		//For unpolarized light, as I am assuming my light source is
		float reflectance = (rs + rp) / 2.0;
		/*
		if(hit.getMaterial()->getName() == "diamond"){
			std::cout << "    Reflectence: " << reflectance << std::endl;
			std::cout << "    rs: " << rs << std::endl;
			std::cout << "    rp: " << rp << std::endl;
		}
		*/
		if(hit.getFace() != NULL)
			assert(reflectance >= 0-eps && reflectance <= 1+eps);

		//Determine if it is reflection or refraction
		bool reflection = false;
		float refract_reflect_random = GLOBAL_mtrand.rand();
		if(refract_reflect_random < reflectance){
			reflection = true;
		}
		
		
		//REFLECTION
		if(reflection){
			Vec3f reflection_direction = ray.getDirection() - 2.0* (ray.getDirection().Dot3(hit_normal) * hit_normal);
			//TODO: double check that reflection always leads to the photon going into air
			//std::cout << "Iter " << iter << ": " << hit.getIsBackfacing() << " " << reflection_direction.AngleBetweenRadians(hit_normal) << std::endl << "    ";
			//printEscapingFacePhoton();
			//Count hits from the back (espcaping and interior bounces)
			
			if(hit.getFace() != NULL){	//for primitives
				if(hit.getIsBackfacing()){
					if(reflection_direction.AngleBetweenRadians(hit_normal) > M_PI/2.0){
						hit.getFace()->incrementNumRaysLeaving();
						hit.getFace()->addLightLeavingDirection(reflection_direction);
						next_n_val = REFRACTIVE_INDEX_OF_AIR;
					}
					else{
						hit.getFace()->incrementNumInteriorBounces();
						next_n_val = current_n_val;
					}
				}
				else{
					if(reflection_direction.AngleBetweenRadians(hit_normal) > M_PI/2.0){
						hit.getFace()->incrementNumRaysEntering();
						next_n_val = nValueSellmeier(wavelength, material->getRefractiveB(), material->getRefractiveC());
					}
					else{
						hit.getFace()->incrementNumRaysReflected();
						next_n_val = current_n_val;
					}
				}
			}
			TracePhoton(bounce_location,reflection_direction,wavelength,iter, viz_color, hit.getMaterial(), current_n_val, single_photon);
		}
		// REFRACTION
		else{

			//std::cout << "Iter " << iter << ": " << hit.getIsBackfacing() << " " << refraction_direction.AngleBetweenRadians(hit_normal) << std::endl << "    ";
			//printEscapingFacePhoton();
			//Count hits from the front (entering and reflecting bounces)

			if(hit.getFace() != NULL){
				if(hit.getIsBackfacing()){

					if(refraction_direction.AngleBetweenRadians(hit_normal) > M_PI/2.0){
						hit.getFace()->incrementNumRaysLeaving();
						hit.getFace()->addLightLeavingDirection(refraction_direction);
						next_n_val = REFRACTIVE_INDEX_OF_AIR;
					}
					else{
						hit.getFace()->incrementNumInteriorBounces();
						next_n_val = current_n_val;
					}
				}
				else{

					if(refraction_direction.AngleBetweenRadians(hit_normal) > M_PI/2.0){
						hit.getFace()->incrementNumRaysEntering();
						next_n_val = nValueSellmeier(wavelength, material->getRefractiveB(), material->getRefractiveC());
					}
					else{
						hit.getFace()->incrementNumRaysReflected();
						next_n_val = current_n_val;
					}
				}
			}
			TracePhoton(bounce_location,refraction_direction,wavelength,iter, Vec4f(1.0, 0.5, 0.0, PHOTON_VISUALIZATION_ALPHA), hit.getMaterial(), next_n_val, single_photon);
		}

	}
	//Visualize photons even when they don't hit anything.
	else {
		if(DRAW_ESCAPING_PHOTONS){
			ray_tree_lock.lock();
			//Ray bounce = Ray(bounce_location, bounce);
			RayTree::AddGeneralSegment(ray, 0, NORMAL_VISUALIZATION_LENGTH * 5, Vec4f(0.0,0.0,0.0,1.0));
			ray_tree_lock.unlock();
		}
		//This is a single photon escaping. Single_photon is only true when
		//'t' is pressed to send a single photon through the scene.
		if(single_photon){
			printEscapingFacePhoton();
			printOutputFile();
		}

		if(iter == 0){
			return false;
		}
	}
	return true;
		
}

// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {
	std::cout << "trace photons" << std::endl;

	RayTree::Activate();
	
	// first, throw away any existing photons
	delete kdtree;
	kdtree = NULL;
	photons.clear();
	
	unsigned long long start_time = 
		std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::system_clock::now().time_since_epoch()).count();
	
	
	int num_threads = args->num_shoot_threads;
	
	// photons emanate from the light sources
	const std::vector<Face*>& lights = mesh->getLights();

	// compute the total area of the lights
	double total_lights_area = 0;
	for (unsigned int i = 0; i < lights.size(); i++) {
		total_lights_area += lights[i]->getArea();
	}
	

	if(num_threads > 1)
	{
		std::thread** threads = new std::thread*[num_threads];
		for(int i = 0; i < num_threads; ++i)
		{
			threads[i] = new std::thread(
				&PhotonMapping::TracePhotonsWorker,
				this,
				std::ref(lights),
				total_lights_area,
				num_threads
			);
		}
		
		// Join the threads back in when they're done.
		for(int i = 0; i < num_threads; ++i)
		{
			threads[i]->join();
			delete threads[i];
			threads[i] = NULL;
		}
		
		delete [] threads;
		threads = NULL;
	}
	else
	{
		// Just 1 thread
		TracePhotonsWorker(lights, total_lights_area, num_threads);
	}
	
	start_time = std::chrono::duration_cast<std::chrono::microseconds>(
		std::chrono::system_clock::now().time_since_epoch()).count() 
		- start_time;

	std::cout << "Photon tracing completed in " << start_time/1000000.0 
		<< " seconds (" << (args->num_photons_to_shoot)/(start_time/1000000.0)
		<< " photons/second) using " << num_threads << " threads\n";
	
	// Construct the KDTree
	makeKDTree();
	
	printEscapingFacePhoton();
	printOutputFile();

	RayTree::Deactivate();
}


void PhotonMapping::TracePhotonsWorker(
	const std::vector<Face*>& lights, 
	double total_lights_area,
	int num_threads
)
{
	

	// shoot a constant number of photons per unit area of light source
	// (alternatively, this could be based on the total energy of each light)
	for (unsigned int i = 0; i < lights.size(); i++) {	
		// Calculate number to shoot per thread
		double my_area = lights[i]->getArea();
		int num = args->num_photons_to_shoot * my_area / total_lights_area / num_threads;
	
		// the initial energy for this photon 
		//Vec3f energy = my_area/double(num) * lights[i]->getMaterial()->getEmittedColor();
		//replace energy with photon color 
		Vec3f normal = lights[i]->computeNormal();
		for (int j = 0; j < num; j++) {

			if (j%50000 == 0){
				std::cout << "Thread: " << std::this_thread::get_id() << "is " << (double)j/(double)num *100 << "% done" << std:: endl;
			}

			bool success = false;

			while(!success){

				Vec3f start = lights[i]->RandomPoint();
				// the initial direction for this photon (for diffuse light sources)
				Vec3f direction = RandomDiffuseDirection(normal);
				if(DIRECTED_LIGHT){
					
					//slightly down (for rectagular prism)
					//direction = Vec3f(1,-0.25,0);

					//slightly up (for pyramid prism)
					//direction = Vec3f(1,0.25,0);

					//directed at diamond from light
					//direction = Vec3f(2.0, -8.0, 2.0);

					//straight up
					//direction = Vec3f(0.0, 1.0, 0.0);

					//left
					//direction = Vec3f(-1.0, 0.0, 0.0);

					// 10 degree rotation
					//direction = Vec3f(-0.17365, -0.984807, 0.0);

					// 15 degree rotation 
					//direction = Vec3f(-0.25882, -0.965926, 0.0);

					// 19.47 degree rotation 
					direction = Vec3f(-0.581396, -0.813621, 0.0);

					// 25 degree rotation (was this 20 or 25? -- I think it was 20)
					//direction = Vec3f(-0.342021, -0.939692, 0.0);

					// 30 degree rotation
					//direction = Vec3f(-0.499998, -0.866027, 0.0);

					// 60 degree rotation
					//direction = Vec3f(-0.866027, -0.499998, 0.0);

					// 120 degree rotation
					//direction = Vec3f(-0.866027, 0.499998, 0.0);

					// 150 degree rotation
					//direction = Vec3f(-0.499998, 0.866027, 0.0);

					//straight down 
					//direction = Vec3f(0.0, -1.0, 0.0);
				}
				else{
					// Floaty face filter

					// 0 degree rotation
					Vec3f a_vec = Vec3f(-1.0, 1.5, -1.0);
					Vec3f b_vec = Vec3f(1.0, 1.5, -1.0);
					Vec3f c_vec = Vec3f(1.0, 1.5, 1.0);
					Vec3f d_vec = Vec3f(-1.0, 1.5, 1.0);

					// 30 degree rotation
					//Vec3f a_vec = Vec3f(-0.116025, 1.79904, -1.0);
					//Vec3f b_vec = Vec3f(1.61603, 0.799038, -1.0);
					//Vec3f c_vec = Vec3f(1.61603, 0.799038, 1.0);
					//Vec3f d_vec = Vec3f(-0.116025, 1.79904, 1.0);

					// 60 degree rotation
					//Vec3f a_vec = Vec3f(0.799038, 1.61603, -1.0);
					//Vec3f b_vec = Vec3f(1.79904, -0.116025, -1.0);
					//Vec3f c_vec = Vec3f(1.79904, -0.116025, 1.0);
					//Vec3f d_vec = Vec3f(0.799038, 1.61603, 1.0);

					// 90 degree rotation
					//Vec3f a_vec = Vec3f(1.5, 1.0, -1.0);
					//Vec3f b_vec = Vec3f(1.5, -1.0, -1.0);
					//Vec3f c_vec = Vec3f(1.5, -1.0, 1.0);
					//Vec3f d_vec = Vec3f(1.5, 1.0, 1.0);

					// 120 degree rotation
					//Vec3f a_vec = Vec3f(1.79904, 0.116025, -1.0);
					//Vec3f b_vec = Vec3f(0.799038, -1.61603, -1.0);
					//Vec3f c_vec = Vec3f(0.799038, -1.61603, 1.0);
					//Vec3f d_vec = Vec3f(1.79904, 0.116025, 1.0);

					// 150 degree rotation
					//Vec3f a_vec = Vec3f(1.61603, -0.799038, -1.0);
					//Vec3f b_vec = Vec3f(-0.116025, -1.79904, -1.0);
					//Vec3f c_vec = Vec3f(-0.116025, -1.79904, 1.0);
					//Vec3f d_vec = Vec3f(1.61603, -0.799038, 1.0);

					// 180 degree rotation
					//Vec3f a_vec = Vec3f(1.0, -1.5, -1.0);
					//Vec3f b_vec = Vec3f(-1.0, -1.5, -1.0);
					//Vec3f c_vec = Vec3f(-1.0, -1.5, 1.0);
					//Vec3f d_vec = Vec3f(1.0, -1.5, 1.0);

					Vertex a_vert (10000, a_vec);
					Vertex b_vert (10001, b_vec);
					Vertex c_vert (10002, c_vec);
					Vertex d_vert (10003, d_vec);

					Vertex* a = &a_vert;
					Vertex* b = &b_vert;
					Vertex* c = &c_vert;
					Vertex* d = &d_vert;

					Material m;

					Face floaty_face(&m);

					Edge e1(a,b,&floaty_face);
					Edge e2(b,c,&floaty_face);
					Edge e3(c,d,&floaty_face);
					Edge e4(d,a,&floaty_face);

					e1.setNext(&e2);
					e2.setNext(&e3);
					e3.setNext(&e4);
					e4.setNext(&e1);

					floaty_face.setEdge(&e1);


					Ray r(start,direction);
					bool * backfacing_hit;
					bool backfacing = false;
					backfacing_hit = &backfacing;
					Hit h;
					if(!(floaty_face.intersect(r,h,true, backfacing_hit))){
						continue;
					}
				}
				//Vec4f photon_color = Vec4f(1.0, 0.0, 1.0, PHOTON_VISUALIZATION_ALPHA);
				Vec4f photon_color = Vec4f(1.0, 1.0, 1.0, PHOTON_VISUALIZATION_ALPHA);
				float wavelength = (GLOBAL_mtrand.rand() * 400) + 380; //random number between 380 and 780 (visible light)

				success = TracePhoton(
					start,
					direction,
					wavelength,
					0, 
					photon_color, 
					NULL, 
					REFRACTIVE_INDEX_OF_AIR, 
					false);
			
			}
		}
	}
}


// Returns true if the point is contained within the sphere
bool inSphere(const Vec3f & center, float radius_squared, const Vec3f & point)
{
	float x = center.x() - point.x();
	float y = center.y() - point.y();
	float z = center.z() - point.z();
	return ((x* x) + (y* y) + (z * z)) < radius_squared;
}


// ======================================================================
// During ray tracing, when a diffuse (or partially diffuse) object is
// hit, gather the nearby photons to approximate indirect illumination

Vec3f PhotonMapping::GatherIndirect(const Vec3f &point, const Vec3f &normal, const Vec3f &direction_from) {

	if (kdtree == NULL) { 
		std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
		return Vec3f(0,0,0); 
	}
	
	/*
	double radius = last_radius;
	
	std::vector<Photon> closest;
	std::vector<Photon> distance;
	
	// Ununsed (not sure why it was there)
	//std::vector<double> radii;
	int i = 1;
	//std::cout << "Collecting Photons for point: " << point << std::endl;
	while(distance.size() < args->num_photons_to_collect){
		//std::cout << "Round " << i << ", radius = " << radius;
		
		closest.clear();
		distance.clear();
		Vec3f bb1 = Vec3f(point.x() + radius, point.y() + radius, point.z() + radius);
		Vec3f bb2 = Vec3f(point.x() - radius, point.y() - radius, point.z() - radius);
		BoundingBox bb = BoundingBox(bb2,bb1);
		kdtree->CollectPhotonsInBox(bb, closest);

		for(unsigned int j = 0; j < closest.size(); j++){
			Photon tmp = closest[j];
			tmp.setPosition(Vec3f(
				tmp.getPosition().x() - point.x(), 
				tmp.getPosition().y() - point.y(), 
				tmp.getPosition().z() - point.z())
			);
			distance.push_back(tmp);
		}
		std::sort( distance.begin(), distance.end(), compareLengths );
		while(distance.size() != 0 
			&& distance[distance.size() - 1].getPosition().Length() > i)
		{
			distance.pop_back();
		}
		
		//std::cout << " and " << closest.size() << " photons found.\n";
		
		// Update radius
		if(radius > 0) 
		{
			radius *=2;
		}
		else 
		{
			radius = args->default_radius;
		}
		++i;
	}
	
	while(distance.size() > args->num_photons_to_collect){
		distance.pop_back();
	}
	if(distance.size() == 0)
		return Vec3f(0,0,0);
	
	radius = distance.back().getPosition().Length();
	
	// Update last radius for next iteration
	last_radius = 1.15 * radius;
	
	// Compute color
	Vec3f color = mixColors(distance);

	double area = PI * pow(radius, 2.0);
	color.Scale(1/(area * distance.size()));
	color.Scale(distance.size()/args->num_photons_to_collect);
	//std::cout << "\n";
	return color;
	*/
	
	
	int num_photons_to_collect = args->num_photons_to_collect;
	float radius = last_radius;
	float radius_squared = radius * radius;
	std::vector<Photon> photons;
	photons.reserve(num_photons_to_collect);
	std::vector<std::pair<Photon, float> > kept_photons;
	kept_photons.reserve(num_photons_to_collect);

	while (1)
	{
		// construct the bounding box
		Vec3f min(point.x() - radius, point.y() - radius, point.z() - radius);
		Vec3f max(point.x() + radius, point.y() + radius, point.z() + radius);
		BoundingBox bb(min, max);
		
		// Collect the photons
		kdtree->CollectPhotonsInBox(bb, photons);


		//std::cout << "Collected " << photons.size() << " photons\n";

		//std::cout << "Bounding box is from " << min << " to " << max << "\n";

		// Throw out the photons not inside the sphere
		for (unsigned int i = 0; i < photons.size();++i)
		{
			//std::cout << "Photon " << itr->getPosition();
			if (inSphere(point, radius_squared, (photons[i]).getPosition()))
			{
				//std::cout << " is outside the sphere.\n";
				//photons.erase(photons.begin() + i);
				
				// Pre-compute distance for use with sorting
				float dist = (photons[i].getPosition() - point).Length();
				kept_photons.push_back(std::make_pair(photons[i], dist));
			}
			// TODO: throw out wrong direction photons
		}

		//std::cout << "Kept " << kept_photons.size() << " photons\n";
	
		// If we don't have enough photons, double the radius and try again. 
		// Otherwise, break out and proceed to the next part.
		if ((int)kept_photons.size() < num_photons_to_collect)
		{
			radius *= 2;
			radius_squared = radius * radius;
			photons.clear();
			kept_photons.clear();
		}
		else
		{
			break;
		}

	}
	
	//std::cout << "The wide radius is " << radius << "\n";

	assert((int)kept_photons.size() >= num_photons_to_collect);

	// Sort the photons by their distance from the point
	std::sort(kept_photons.begin(), kept_photons.end(), [&point](const std::pair<Photon, float> & a, const std::pair<Photon, float> & b) -> bool{
		return a.second < b.second;
	});

	// Drop the extra photons - we don't need them
	// taking this out
	kept_photons.resize(num_photons_to_collect);


	//std::cout << "After dropping, now have " << kept_photons.size() << " photons\n";

	// Find out the radius actually required to capture these photons by
	// taking length of the farthest photon from our point.
	radius = (kept_photons.back().first.getPosition() - point).Length();

	//std::cout << "The actual radius is " << radius << "\n";

	// Update the "last radius" for the next call to this function. Slightly
	// overestimate, since it's better to sort a few more photons than get more
	// from the KD tree.
	last_radius = 1.15 * radius;

	// Compute color
	Vec3f color = mixColors(kept_photons);

	double area = PI * radius * radius;
	//color.Scale(1/(area));
	color.Scale(1/(area* kept_photons.size()));
	
	//color.Scale(kept_photons.size()/args->num_photons_to_collect);
	//std::cout << "\n";
	return color;
}

void PhotonMapping::makeKDTree()
{
	// If not waiting to construct a balanced tree
	if(!args->balanced_tree)
	{
		std::cout << "Creating unbalanced KDTree...\n";

		unsigned long long start_time = 
			std::chrono::duration_cast<std::chrono::microseconds>(
				std::chrono::system_clock::now().time_since_epoch()).count();
		
		
		// consruct a kdtree to store the photons
		BoundingBox *bb = mesh->getBoundingBox();
		Vec3f min = bb->getMin();
		Vec3f max = bb->getMax();
		Vec3f diff = max-min;
		min -= 0.001*diff;
		max += 0.001*diff;
		kdtree = new KDTree(BoundingBox(min,max));
		
		// Add all the photons, 1 by 1.
		for(unsigned int i = 0; i < photons.size(); ++i)
		{
			kdtree->AddPhoton(photons[i]);
		}
		start_time = std::chrono::duration_cast<std::chrono::microseconds>(
				std::chrono::system_clock::now().time_since_epoch()).count() 
				- start_time;
		std::cout << "Creation completed in " << start_time/1000000.0 << " seconds.\n";
	}
	else
	{
		std::cout << "Creating balanced KDTree...\n";
		
		unsigned long long start_time = 
			std::chrono::duration_cast<std::chrono::microseconds>(
				std::chrono::system_clock::now().time_since_epoch()).count();
		
		
		// Find the bounding box for the tree
		BoundingBox *bb = mesh->getBoundingBox();
		Vec3f min = bb->getMin();
		Vec3f max = bb->getMax();
		Vec3f diff = max-min;
		min -= 0.001*diff;
		max += 0.001*diff;
		
		// Construct the balanced tree using the special constructor
		kdtree = new KDTree(BoundingBox(min, max), 0, photons);
		
		// Throw out the vector representation of photons to save memory
		photons.clear();
		
		start_time = std::chrono::duration_cast<std::chrono::microseconds>(
				std::chrono::system_clock::now().time_since_epoch()).count() 
				- start_time;
		std::cout << "Creation completed in " << start_time/1000000.0 
			<< " seconds.\n";
	}
}

// Print the kdtree
void PhotonMapping::printKDTree() const
{
	std::vector<int> nums;
	kdtree->printLeaves(0, nums);
	std::sort(nums.begin(), nums.end());
	/*
	for(int i = 0; i < nums.size(); ++i)
	{
		std::cout << nums[i] << "\n";
	}
	*/
	int max_ = 0;
	int count = 0;
	std::vector<int> stats(13, 0);
	for(unsigned int i = 0; i < nums.size(); ++i)
	{
		if(nums[i] >= 0 && nums[i] <= 10)
		{
			stats[0]++;
		}
		else if(nums[i] > 10 && nums[i] <= 20)
		{
			stats[1]++;
		}
		else if(nums[i] > 20 && nums[i] <= 30)
		{
			stats[2]++;
		}
		else if(nums[i] > 30 && nums[i] <= 40)
		{
			stats[3]++;
		}
		else if(nums[i] > 40 && nums[i] <= 50)
		{
			stats[4]++;
		}
		else if(nums[i] > 50 && nums[i] <= 60)
		{
			stats[5]++;
		}
		else if(nums[i] > 60 && nums[i] <= 70)
		{
			stats[6]++;
		}
		else if(nums[i] > 70 && nums[i] <= 80)
		{
			stats[7]++;
		}
		else if(nums[i] > 80 && nums[i] <= 90)
		{
			stats[8]++;
		}
		else if(nums[i] > 90 && nums[i] <= 100)
		{
			stats[9]++;
		}
		else if(nums[i] > 100 && nums[i] <= 400)
		{
			stats[10]++;
		}
		else if(nums[i] > 400 && nums[i] <= 700)
		{
			stats[11]++;
		}
		else if(nums[i] > 700)
		{
			stats[12]++;
		}
		max_ = std::max(max_, nums[i]);
		count += nums[i];
	}
	
	std::cout << "There are " << stats[0] << " Nodes with 0 -> 10 photons\n";
	std::cout << "There are " << stats[1] << " Nodes with 11 -> 20 photons\n";
	std::cout << "There are " << stats[2] << " Nodes with 21 -> 30 photons\n";
	std::cout << "There are " << stats[3] << " Nodes with 31 -> 40 photons\n";
	std::cout << "There are " << stats[4] << " Nodes with 41 -> 50 photons\n";
	std::cout << "There are " << stats[5] << " Nodes with 51 -> 60 photons\n";
	std::cout << "There are " << stats[6] << " Nodes with 61 -> 70 photons\n";
	std::cout << "There are " << stats[7] << " Nodes with 71 -> 80 photons\n";
	std::cout << "There are " << stats[8] << " Nodes with 81 -> 90 photons\n";
	std::cout << "There are " << stats[9] << " Nodes with 91 -> 100 photons\n";
	std::cout << "There are " << stats[10] << " Nodes with 101 -> 400 photons\n";
	std::cout << "There are " << stats[11] << " Nodes with 401 -> 700 photons\n";
	std::cout << "There are " << stats[12] << " Nodes with 701+ photons\n";
	std::cout << "There are " << count << " total photons.\n";
	
	
	return;
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
			//TODO: check this with Barb
			//if we have a backfacing hit, the light is escaping via that face
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

	//Commented out because I already do this in the Trace Photon Function
	// Count the number of rays leaving the faces
	//if(h.getIsBackfacing() && h.getFace() != NULL){
	//	std::cout << "that happened" << std::endl;
	//	h.getFace()->incrementNumRaysLeaving();
	//}

	return answer;
}
