#include "glCanvas.h"

#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "glCanvas.h"
#include "sphere.h"
#include "raytree.h"
#include "raytracer.h"
#include "utils.h"

#define eps 0.1

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh *m, ArgParser *a) {
	mesh = m;
	args = a;
	num_faces = -1;	
	formfactors = NULL;
	area = NULL;
	undistributed = NULL;
	absorbed = NULL;
	radiance = NULL;
	max_undistributed_patch = -1;
	total_area = -1;
	Reset();
}

Radiosity::~Radiosity() {
	Cleanup();
	cleanupVBOs();
}

void Radiosity::Cleanup() {
	delete [] formfactors;
	delete [] area;
	delete [] undistributed;
	delete [] absorbed;
	delete [] radiance;
	num_faces = -1;
	formfactors = NULL;
	area = NULL;
	undistributed = NULL;
	absorbed = NULL;
	radiance = NULL;
	max_undistributed_patch = -1;
	total_area = -1;
}

void Radiosity::Reset() {
	delete [] area;
	delete [] undistributed;
	delete [] absorbed;
	delete [] radiance;

	// create and fill the data structures
	num_faces = mesh->numFaces();
	area = new double[num_faces];
	undistributed = new Vec3f[num_faces];
	absorbed = new Vec3f[num_faces];
	radiance = new Vec3f[num_faces];
	for (int i = 0; i < num_faces; i++) {
		Face *f = mesh->getFace(i);
		f->setRadiosityPatchIndex(i);
		setArea(i,f->getArea());
		Vec3f emit = f->getMaterial()->getEmittedColor();
		setUndistributed(i,emit);
		setAbsorbed(i,Vec3f(0,0,0));
		setRadiance(i,emit);
	}

	// find the patch with the most undistributed energy
	findMaxUndistributed();
}


// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
	// find the patch with the most undistributed energy 
	// don't forget that the patches may have different sizes!
	max_undistributed_patch = -1;
	total_undistributed = 0;
	total_area = 0;
	double max = -1;
	for (int i = 0; i < num_faces; i++) {
		double m = getUndistributed(i).Length() * getArea(i);
		total_undistributed += m;
		total_area += getArea(i);
		if (max < m) {
			max = m;
			max_undistributed_patch = i;
		}
	}
	assert (max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}


void Radiosity::ComputeFormFactors() {
	assert (formfactors == NULL);
	assert (num_faces > 0);
	formfactors = new double[num_faces*num_faces];
	//this is the integral method -> I should probably be using the ray casting method
	for(int face1 = 0; face1 < num_faces; face1++){
		double form_factor_sum = 0.0;
		for(int face2 = 0; face2 < num_faces; face2++){
			//a face shouldn't absorb any of it's own energy.
			if(face1 == face2){
				setFormFactor(face1, face2, 0.0);
			}
			//two different faces -> find form factor
			else{
	
				//find the angles between the normals and a vector connecting the center points
				Vec3f connecting_vector = mesh->getFace(face2)->computeCentroid() - mesh->getFace(face1)->computeCentroid();
				double r = connecting_vector.Length();
				connecting_vector.Normalize();
				double cos_angle_face_1 = mesh->getFace(face1)->computeNormal().Dot3(connecting_vector);
				
				connecting_vector *= -1.0; //make it go the other way for the other face
				double cos_angle_face_2 = mesh->getFace(face2)->computeNormal().Dot3(connecting_vector);
				
				double form_factor =  ( (cos_angle_face_1 * cos_angle_face_2)/(pow(r,2) * M_PI ) )  * mesh->getFace(face2)->getArea() ;

				//this means that the angles is > 90 and <270, so the patch is behind you
				if(cos_angle_face_1 < 0 || cos_angle_face_2 < 0){
					form_factor = 0;
				}
				if(form_factor > 1.0) {
					std::cout << "Form factor calculation" << std::endl;
					std::cout << "    Form factor greater than 1!!! " << form_factor << std::endl;
					std::cout << "    first cos is" << cos_angle_face_1 << std::endl;
					std::cout << "    second cos is" << cos_angle_face_2 << std::endl;
					std::cout << "    distance between is " << r << std::endl;
				}
				form_factor_sum += form_factor;



				//only detect occlusion of shadow samples > 0
				if(args->num_shadow_samples == 0){
					setFormFactor(face1, face2, form_factor);
				}
				else{
					connecting_vector *= -1.0; //put it back to going the right way
					Hit hit_for_shadow = Hit();
					Ray connecting_ray = Ray(mesh->getFace(face1)->computeCentroid(), connecting_vector);
					raytracer->CastRay(connecting_ray, hit_for_shadow, true);
					if( (hit_for_shadow.getT() * connecting_ray.getDirection()).Length() - r < eps && r - (hit_for_shadow.getT() * connecting_ray.getDirection()).Length() < eps ) {
						setFormFactor(face1, face2, form_factor);
			//			std::cout << "Hit the object, no occlusion" << std::endl;
					}
					else{
			//			std::cout << "occlusion!! r was: " << r << " and the ray length was: " << (hit_for_shadow.getT() * connecting_ray.getDirection()).Length() << std::endl;

						setFormFactor(face1, face2, 0.0);
					}
					


				}
			}
		}

		//NOTE: Barb has already normalized this somewhere else
		//loop back through to normalize the sum of the form factor to be 1 (divide each number by the sum)
		double normalized_form_factor_sum = 0.0;
		for(int face2 = 0; face2 < num_faces;  face2++){
			setFormFactor(face1, face2, getFormFactor(face1, face2)/form_factor_sum);
		//	if(face1 != face2)
		//		setFormFactor(face1, face2, 1/ (double)num_faces);
		//	else 
		//		setFormFactor(face1, face2, 0);
			normalized_form_factor_sum+= getFormFactor(face1, face2)/form_factor_sum;
		}
//		std::cout << "Sum of form factors for face " << face1 << " is " << form_factor_sum << std::endl;
	}

	

	// =====================================
	// COMPUTE THE FORM FACTORS
	// =====================================

	std::cout << "Form factors calculated!" << std::endl;

}

// ================================================================
// ================================================================

double Radiosity::Iterate() {
	if (formfactors == NULL) {
		ComputeFormFactors();
		//this is the first iteration, so for each face, get the emitted light
		for(int f = 0; f < num_faces; f++){
			setUndistributed(f, getRadiance(f)); //light starts as radiance
		}
	}
	assert (formfactors != NULL);




	// ==========================================
	// ASSIGNMENT:	IMPLEMENT RADIOSITY ALGORITHM
	// ==========================================

	
	//Progressive Refinement with Ambient Term
	bool ambient = false;
	Vec3f total_undistributed_light = Vec3f(0.0,0.0,0.0);
	Vec3f undistributed_per_face;
	if(args->progressive_ambient){
		//loop through the faces to find the total undistributed light
		
		for(int u = 0; u < num_faces; u++){
			total_undistributed_light += getUndistributed(u);		
		}
		//Go through and get rid of the total undistributed before next iteration starts
		undistributed_per_face = total_undistributed_light/num_faces;
		for(int u = 0; u < num_faces; u++){
			setRadiance(u, getRadiance(u) - undistributed_per_face);
			//this checks if it was the first iteration where all of the radiance is 0, but there is undistributed light
			if(getRadiance(u).r() < -0.01 || getRadiance(u).b() < -0.01 ||getRadiance(u).b() < -0.01){
				setRadiance(u, Vec3f(0.0,0.0,0.0));
			}
		}
	}

	//get the max undistributed amount of light and which patch it's in
	findMaxUndistributed();
//	std::cout << "New iteration" << std::endl;
//	std::cout << "    Index of patch with most light:  " << max_undistributed_patch << std::endl;
	
	Vec3f undistributed_at_max = getUndistributed(max_undistributed_patch) * getArea(max_undistributed_patch);

//	std::cout << "    Amount of light in that patch: (" << undistributed_at_max.r() << ", " << undistributed_at_max.g() << ", " << undistributed_at_max.g() << ")" << std::endl;  
	
	//go through and distribute the light to the other faces based on the form factors
	for(int u = 0; u < num_faces; u++){
		Vec3f reflected_color = mesh->getFace(u)->getMaterial()->getDiffuseColor();
//		std::cout << "         reflected color is " << reflected_color.r() << ", " << reflected_color.g() << ", " << reflected_color.g() << ")" << std::endl; 
		Vec3f absorbed_color = Vec3f(1.0,1.0,1.0) - reflected_color;
		setUndistributed(u, getUndistributed(u) + (reflected_color * (undistributed_at_max*getFormFactor(max_undistributed_patch, u))/getArea(u) ) );
		setRadiance(u, getRadiance(u) + (reflected_color * (undistributed_at_max*getFormFactor(max_undistributed_patch, u))/getArea(u) ) );
		setAbsorbed(u, getAbsorbed(u) + (absorbed_color * (undistributed_at_max*getFormFactor(max_undistributed_patch, u))/getArea(u) ) );
	}

	//zero out the undistributed light for this patch
	setUndistributed(max_undistributed_patch, Vec3f(0.0,0.0,0.0));
//	std::cout << "    After clearing the light, the max_undistributed_patch has " << getUndistributed(max_undistributed_patch).r() << ", " << getUndistributed(max_undistributed_patch).g() << ", " << getUndistributed(max_undistributed_patch).g() << ")" << std::endl;  

	//this updates the variables total_undistributed
	findMaxUndistributed();

	//Progressive Refinement with Ambient Term (the part where you add it in)
	if(args->progressive_ambient){
		total_undistributed_light = Vec3f(0.0,0.0,0.0);
		for(int u = 0; u < num_faces; u++){
			total_undistributed_light += getUndistributed(u);		
		}
		//Go through and get rid of the total undistributed before next iteration starts
		undistributed_per_face = total_undistributed_light/num_faces;
		for(int u = 0; u < num_faces; u++){
			setRadiance(u, getRadiance(u) + undistributed_per_face);
		}
	}
	// return the total light yet undistributed
	// (so we can decide when the solution has sufficiently converged)
//	std::cout << "    total undistributed light is " << total_undistributed << " at the end of an iteration" << std::endl;
	return total_undistributed;
}


// =======================================================================================
// VBO & DISPLAY FUNCTIONS
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex *have, Face *f, std::vector<Face*> &faces) {
	for (unsigned int i = 0; i < faces.size(); i++) {
		if (faces[i] == f) return;
	}
	if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3]) return;
	faces.push_back(f);
	for (int i = 0; i < 4; i++) {
		Edge *ea = f->getEdge()->getOpposite();
		Edge *eb = f->getEdge()->getNext()->getOpposite();
		Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
		Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
		if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),faces);
		if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),faces);
		if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),faces);
		if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),faces);
	}
}

// different visualization modes
Vec3f Radiosity::setupHelperForColor(Face *f, int i, int j) {
	assert (mesh->getFace(i) == f);
	assert (j >= 0 && j < 4);
	if (args->render_mode == RENDER_MATERIALS) {
		return f->getMaterial()->getDiffuseColor();
	} else if (args->render_mode == RENDER_RADIANCE && args->interpolate == true) {
		std::vector<Face*> faces;
		CollectFacesWithVertex((*f)[j],f,faces);
		double total = 0;
		Vec3f color = Vec3f(0,0,0);
		Vec3f normal = f->computeNormal();
		for (unsigned int i = 0; i < faces.size(); i++) {
			Vec3f normal2 = faces[i]->computeNormal();
			double area = faces[i]->getArea();
			if (normal.Dot3(normal2) < 0.5) continue;
			assert (area > 0);
			total += area;
			color += area * getRadiance(faces[i]->getRadiosityPatchIndex());
		}
		assert (total > 0);
		color /= total;
		return color;
	} else if (args->render_mode == RENDER_LIGHTS) {
		return f->getMaterial()->getEmittedColor();
	} else if (args->render_mode == RENDER_UNDISTRIBUTED) { 
		return getUndistributed(i);
	} else if (args->render_mode == RENDER_ABSORBED) {
		return getAbsorbed(i);
	} else if (args->render_mode == RENDER_RADIANCE) {
		return getRadiance(i);
	} else if (args->render_mode == RENDER_FORM_FACTORS) {
		if (formfactors == NULL) ComputeFormFactors();
		double scale = 0.2 * total_area/getArea(i);
		double factor = scale * getFormFactor(max_undistributed_patch,i);
		return Vec3f(factor,factor,factor);
	} else {
		assert(0);
	}
	exit(0);
}


void Radiosity::initializeVBOs() {
	// create a pointer for the vertex & index VBOs
	glGenBuffers(1, &mesh_quad_verts_VBO);
	glGenBuffers(1, &mesh_quad_indices_VBO);
	glGenBuffers(1, &mesh_textured_quad_indices_VBO);
	glGenBuffers(1, &mesh_interior_edge_indices_VBO);
	glGenBuffers(1, &mesh_border_edge_indices_VBO);
}


void Radiosity::setupVBOs() {
	mesh_quad_verts.clear();
	mesh_quad_indices.clear();
	mesh_textured_quad_indices.clear();
	mesh_border_edge_indices.clear();
	mesh_interior_edge_indices.clear();

	// initialize the data in each vector
	int num_faces = mesh->numFaces();
	assert (num_faces > 0);
	for (int i = 0; i < num_faces; i++) {
		Face *f = mesh->getFace(i);
		Edge *e = f->getEdge();
		for (int j = 0; j < 4; j++) {
			Vec3f pos = ((*f)[j])->get();
			Vec3f normal = f->computeNormal();
			Vec3f color = setupHelperForColor(f,i,j);
			color = Vec3f(linear_to_srgb(color.r()),
				linear_to_srgb(color.g()),
				linear_to_srgb(color.b()));
			mesh_quad_verts.push_back(VBOPosNormalColorTexture(pos,normal,color,(*f)[j]->get_s(),(*f)[j]->get_t()));
			if (e->getOpposite() == NULL) { 
	mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+j,i*4+(j+1)%4));
			} else if (e->getStartVertex()->getIndex() < e->getEndVertex()->getIndex()) {
	mesh_interior_edge_indices.push_back(VBOIndexedEdge(i*4+j,i*4+(j+1)%4));
			}
			e = e->getNext();
		}
		if (f->getMaterial()->hasTextureMap()) {
			mesh_textured_quad_indices.push_back(VBOIndexedQuad(i*4,i*4+1,i*4+2,i*4+3));
		} else {
			mesh_quad_indices.push_back(VBOIndexedQuad(i*4,i*4+1,i*4+2,i*4+3));
		}
		// also outline the max_undistributed patch
		if (args->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
			mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+0,i*4+1));
			mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+1,i*4+2));
			mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+2,i*4+3));
			mesh_border_edge_indices.push_back(VBOIndexedEdge(i*4+3,i*4+0));
		}
	}
	assert ((int)mesh_quad_verts.size() == num_faces*4);
	assert ((int)mesh_quad_indices.size() + (int)mesh_textured_quad_indices.size() == num_faces);

	// cleanup old buffer data (if any)
	cleanupVBOs();

	// copy the data to each VBO
	glBindBuffer(GL_ARRAY_BUFFER,mesh_quad_verts_VBO); 
	glBufferData(GL_ARRAY_BUFFER,
				 sizeof(VBOPosNormalColorTexture) * num_faces * 4,
				 &mesh_quad_verts[0],
				 GL_STATIC_DRAW); 
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_quad_indices_VBO); 
	glBufferData(GL_ELEMENT_ARRAY_BUFFER,
				 sizeof(VBOIndexedQuad) * mesh_quad_indices.size(),
				 &mesh_quad_indices[0], GL_STATIC_DRAW);
	if (mesh_textured_quad_indices.size() > 0) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_textured_quad_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedQuad) * mesh_textured_quad_indices.size(),
		 &mesh_textured_quad_indices[0], GL_STATIC_DRAW);
	}
	if (mesh_interior_edge_indices.size() > 0) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_interior_edge_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * mesh_interior_edge_indices.size(),
		 &mesh_interior_edge_indices[0], GL_STATIC_DRAW);
	}
	if (mesh_border_edge_indices.size() > 0) {
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,mesh_border_edge_indices_VBO); 
		glBufferData(GL_ELEMENT_ARRAY_BUFFER,
		 sizeof(VBOIndexedEdge) * mesh_border_edge_indices.size(),
		 &mesh_border_edge_indices[0], GL_STATIC_DRAW);
	}

	// WARNING: this naive VBO implementation only allows a single texture
	int num_textured_materials = 0;
	for (unsigned int mat = 0; mat < mesh->materials.size(); mat++) {
		Material *m = mesh->materials[mat];
		if (m->hasTextureMap()) {
			glBindTexture(GL_TEXTURE_2D,m->getTextureID());
			num_textured_materials++;
		}
	}
	assert (num_textured_materials <= 1);
}


void Radiosity::drawVBOs() {
	// =====================
	// DRAW ALL THE POLYGONS
	if (args->render_mode == RENDER_MATERIALS) {
		glEnable(GL_LIGHTING);
	} else {
		glDisable(GL_LIGHTING);
	}
	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.1,4.0);
	int num_faces = mesh->numFaces();
	assert ((int)mesh_quad_indices.size() + (int)mesh_textured_quad_indices.size() == num_faces);

	glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
	glEnableClientState(GL_VERTEX_ARRAY);
	glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
	glEnableClientState(GL_NORMAL_ARRAY);
	glNormalPointer(GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(12));
	glEnableClientState(GL_COLOR_ARRAY);
	glColorPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(24));

	// draw non textured faces
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_quad_indices_VBO);
	glDrawElements(GL_QUADS, 
		 mesh_quad_indices.size()*4,
		 GL_UNSIGNED_INT,
		 BUFFER_OFFSET(0));

	// draw textured faces
	if (args->render_mode == RENDER_MATERIALS) {
		glEnable(GL_TEXTURE_2D);
		glEnableClientState(GL_TEXTURE_COORD_ARRAY);
		glTexCoordPointer( 2, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(36));
	}
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_textured_quad_indices_VBO);
	glDrawElements(GL_QUADS, 
		 mesh_textured_quad_indices.size()*4,
		 GL_UNSIGNED_INT,
		 BUFFER_OFFSET(0));
	if (args->render_mode == RENDER_MATERIALS) {
		glDisableClientState(GL_TEXTURE_COORD_ARRAY);
		glDisable(GL_TEXTURE_2D);
	}

	glDisableClientState(GL_NORMAL_ARRAY);
	glDisableClientState(GL_COLOR_ARRAY);
	glDisableClientState(GL_VERTEX_ARRAY);

	glDisable(GL_POLYGON_OFFSET_FILL);	


	// =====================
	// DRAW WIREFRAME
	if (args->wireframe) {
		glDisable(GL_LIGHTING);
		if (mesh_interior_edge_indices.size() > 0) {
			glLineWidth(1);
			glColor3f(0,0,0);
			glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_interior_edge_indices_VBO);
			glDrawElements(GL_LINES, mesh_interior_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
		if (mesh_border_edge_indices.size() > 0) {
			glLineWidth(3);
			glColor3f(1,0,0);
			glBindBuffer(GL_ARRAY_BUFFER, mesh_quad_verts_VBO);
			glEnableClientState(GL_VERTEX_ARRAY);
			glVertexPointer(3, GL_FLOAT, sizeof(VBOPosNormalColorTexture), BUFFER_OFFSET(0));
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh_border_edge_indices_VBO);
			glDrawElements(GL_LINES, mesh_border_edge_indices.size()*2, GL_UNSIGNED_INT, 0);
			glDisableClientState(GL_VERTEX_ARRAY);
		}
	}
	HandleGLError(); 
}


void Radiosity::cleanupVBOs() {
	glDeleteBuffers(1, &mesh_quad_verts_VBO);
	glDeleteBuffers(1, &mesh_quad_indices_VBO);
	glDeleteBuffers(1, &mesh_textured_quad_indices_VBO);
	glDeleteBuffers(1, &mesh_interior_edge_indices_VBO);
	glDeleteBuffers(1, &mesh_border_edge_indices_VBO);
}

