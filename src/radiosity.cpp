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
	if(args->color_by_normal == true){
		
		Vec3f normal = f->computeNormal();
		Vec3f up(0.0f, 1.0f, 0.0f);

		Vec3f blue(0.1f, 0.1f, 1.0f);
		Vec3f orange(0.8f, 0.2f, 0.0f);
		Vec3f grey(0.25f, 0.25f, 0.25f);
		Vec3f yellow(1.0f, 1.0f, 0.0f);

		float angle = normal.AngleBetweenDegrees(up);

		//std::cout << "Normal: " << normal.r() << " " << normal.g() << " " << normal.b() << std::endl;
		//std::cout << "    angle: " << angle << std::endl;
		if(angle < 45){
			std::cout << "    color: blue" << std::endl;
			return blue;
		}
		else if (angle <90){
			std::cout << "    color: orange" << std::endl;
			return orange;
		}
		else if (angle <135){
			std::cout << "    color: grey" << std::endl;
			return grey;
		}
		else{
			std::cout << "    color: yellow" << std::endl;
			return yellow;
		}

		
	}
	else if (args->render_mode == RENDER_MATERIALS) {
		return f->getMaterial()->getDiffuseColor();
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

