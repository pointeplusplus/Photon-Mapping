#include "glCanvas.h"
#include "argparser.h"
#include "boundingbox.h"
#include "camera.h"
#include "radiosity.h"
#include "raytracer.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "raytree.h"
#include "utils.h"

#include <thread>

#ifndef REFRACTIVE_INDEX_OF_AIR
#define REFRACTIVE_INDEX_OF_AIR 1.000293
#endif

// ========================================================
// static variables of GLCanvas class

ArgParser* GLCanvas::args = NULL;
Mesh* GLCanvas::mesh = NULL;
Radiosity* GLCanvas::radiosity = NULL;
RayTracer* GLCanvas::raytracer = NULL;
PhotonMapping* GLCanvas::photon_mapping = NULL;

// State of the mouse cursor
int GLCanvas::mouseButton = 0;
int GLCanvas::mouseX = 0;
int GLCanvas::mouseY = 0;
bool GLCanvas::controlPressed = false;
bool GLCanvas::shiftPressed = false;
bool GLCanvas::altPressed = false;

// params for the raytracing animation
int GLCanvas::raytracing_x;
int GLCanvas::raytracing_y;
int GLCanvas::raytracing_skip;

// For timing animation
time_t GLCanvas::rendering_time;

// ========================================================
// Initialize all appropriate OpenGL variables, set
// callback functions, and start the main event loop.
// This function will not return but can be terminated
// by calling 'exit(0)'
// ========================================================

void GLCanvas::initialize(ArgParser *_args, Mesh *_mesh, 
				RayTracer *_raytracer, Radiosity *_radiosity, PhotonMapping *_photon_mapping) {

	args = _args;
	mesh = _mesh;
	raytracer = _raytracer;
	radiosity = _radiosity;
	photon_mapping = _photon_mapping;

	// setup glut stuff
	glutInitWindowSize(args->width, args->height);
	glutInitWindowPosition(100,100);
	glutInitDisplayMode(GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGB);
	glutCreateWindow("OpenGL Viewer");
	HandleGLError("in glcanvas initialize");

#ifdef _WIN32
	GLenum err = glewInit();
	if (err != GLEW_OK) {
			fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
			exit(1);
	}
#endif
	// basic rendering 
	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_NORMALIZE);
	glShadeModel(GL_SMOOTH);
	glLightModeli(GL_LIGHT_MODEL_LOCAL_VIEWER, GL_TRUE);
	GLfloat ambient[] = { 0.2, 0.2, 0.2, 1.0 };
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_FALSE);
	glCullFace(GL_BACK);
	glDisable(GL_CULL_FACE);

	// Initialize callback functions
	glutMouseFunc(mouse);
	glutMotionFunc(motion);
	glutDisplayFunc(display);
	glutReshapeFunc(reshape);
	glutKeyboardFunc(keyboard);
	glutIdleFunc(idle);

	HandleGLError("finished glcanvas initialize");

	RayTree::initializeVBOs();
	if (radiosity) radiosity->initializeVBOs();
	if (photon_mapping) photon_mapping->initializeVBOs();

	RayTree::setupVBOs();
	if (radiosity) radiosity->setupVBOs();
	if (photon_mapping) photon_mapping->setupVBOs();

	HandleGLError("finished glcanvas initialize");

	// Enter the main rendering loop
	glutMainLoop();
}


// ========================================================

void GLCanvas::InitLight() {
	// Set the last component of the position to 0 to indicate
	// a directional light source

	GLfloat position[4] = { 30,30,100, 1};
	GLfloat diffuse[4] = { 0.75,0.75,0.75,1};
	GLfloat specular[4] = { 0,0,0,1};
	GLfloat ambient[4] = { 0.2, 0.2, 0.2, 1.0 };

	GLfloat zero[4] = {0,0,0,0};
	glLightfv(GL_LIGHT1, GL_POSITION, position);
	glLightfv(GL_LIGHT1, GL_DIFFUSE, diffuse);
	glLightfv(GL_LIGHT1, GL_SPECULAR, specular);
	glLightfv(GL_LIGHT1, GL_AMBIENT, zero);
	glEnable(GL_LIGHT1);
	glLightModeli(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE);
	glEnable(GL_COLOR_MATERIAL);
	glLightModelfv(GL_LIGHT_MODEL_AMBIENT, ambient);

	GLfloat spec_mat[4] = {1,1,1,1};
	float glexponent = 30;
	glMaterialfv(GL_FRONT, GL_SHININESS, &glexponent);
	glMaterialfv(GL_FRONT, GL_SPECULAR, spec_mat);

	glColorMaterial(GL_FRONT, GL_AMBIENT_AND_DIFFUSE);
	float back_color[] = { 0.2,0.8,0.8,1};
	glMaterialfv(GL_BACK, GL_AMBIENT_AND_DIFFUSE, back_color);
	glEnable(GL_LIGHT1);
}


void GLCanvas::display(void) {
	glDrawBuffer(GL_BACK);

	Vec3f bg = mesh->background_color;
	// Clear the display buffer, set it to the background color
	glClearColor(bg.r(),bg.g(),bg.b(),0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// Set the camera parameters
	mesh->camera->glInit(args->width, args->height);
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	mesh->camera->glPlaceCamera();
	InitLight(); // light will be a headlamp!

	if (args->intersect_backfacing)
		glDisable(GL_CULL_FACE);
	else
		glEnable(GL_CULL_FACE);

	glEnable(GL_LIGHTING);
	glEnable(GL_DEPTH_TEST);
	
	//	glCallList(display_list_index);
	HandleGLError(); 

	radiosity->drawVBOs();
	photon_mapping->drawVBOs();
	RayTree::drawVBOs();
	 
	// Swap the back buffer with the front buffer to display
	// the scene
	glutSwapBuffers();
}

// ========================================================
// Callback function for window resize
// ========================================================

void GLCanvas::reshape(int w, int h) {
	args->width = w;
	args->height = h;

	// Set the OpenGL viewport to fill the entire window
	glViewport(0, 0, (GLsizei)args->width, (GLsizei)args->height);

	// Set the camera parameters to reflect the changes
	mesh->camera->glInit(args->width, args->height);
}

// ========================================================
// Callback function for mouse click or release
// ========================================================

void GLCanvas::mouse(int button, int /*state*/, int x, int y) {
	args->raytracing_animation = false;
	// Save the current state of the mouse.	This will be
	// used by the 'motion' function
	mouseButton = button;
	mouseX = x;
	mouseY = y;

	shiftPressed = (glutGetModifiers() & GLUT_ACTIVE_SHIFT) != 0;
	controlPressed = (glutGetModifiers() & GLUT_ACTIVE_CTRL) != 0;
	altPressed = (glutGetModifiers() & GLUT_ACTIVE_ALT) != 0;
}

// ========================================================
// Callback function for mouse drag
// ========================================================

void GLCanvas::motion(int x, int y) {
	// Left button = rotation
	// (rotate camera around the up and horizontal vectors)
	if (mouseButton == GLUT_LEFT_BUTTON) {
		mesh->camera->rotateCamera(0.005*(mouseX-x), 0.005*(mouseY-y));
		mouseX = x;
		mouseY = y;
	}
	// Middle button = translation
	// (move camera perpendicular to the direction vector)
	else if (mouseButton == GLUT_MIDDLE_BUTTON) {
		mesh->camera->truckCamera((mouseX-x)*0.5, (y-mouseY)*0.5);
		mouseX = x;
		mouseY = y;
	}
	// Right button = dolly or zoom
	// (move camera along the direction vector)
	else if (mouseButton == GLUT_RIGHT_BUTTON) {
		if (controlPressed) {
			mesh->camera->zoomCamera(mouseY-y);
		} else {
			mesh->camera->dollyCamera(mouseY-y);
		}
		mouseX = x;
		mouseY = y;
	}

	// Redraw the scene with the new camera parameters
	glutPostRedisplay();
}

// ========================================================
// Callback function for keyboard events
// ========================================================

void GLCanvas::keyboard(unsigned char key, int x, int y) {
	args->raytracing_animation = false;
	switch (key) {
		// RAYTRACING STUFF
	case 'r':	case 'R':
		// animate raytracing of the scene
		args->gather_indirect=false;
		args->raytracing_animation = !args->raytracing_animation;
		if (args->raytracing_animation) {
			// time the animation
			rendering_time = time(NULL);
			
			raytracing_skip = my_max(args->width,args->height) / 10;
			if (raytracing_skip % 2 == 0) raytracing_skip++;
			assert (raytracing_skip >= 1);
			raytracing_x = raytracing_skip/2;
			raytracing_y = raytracing_skip/2;
			display(); // clear out any old rendering
			printf ("raytracing animation started, press 'R' to stop\n");
		} else {
			printf ("raytracing animation stopped, press 'R' to start\n");		
			// Print time to render
			rendering_time =  time(NULL) - rendering_time;
			std::cout << "Rendering completed in " << rendering_time 
					<< " seconds." << std::endl;
		}
		break;
	case 't':	case 'T': {
		// visualize the ray tree for the pixel at the current mouse position
		int i = x;
		int j = glutGet(GLUT_WINDOW_HEIGHT)-y;
		RayTree::Activate();
		raytracing_skip = 1;
		//TraceRay(i,j);
		//photon_mapping->setupVBOs();
		//photon_mapping->initializeVBOs();
		std::cout << "about to trace photons " << std::endl;
		TracePhoton(i,j);
		//photon_mapping->drawVBOs();
		//RayTree::drawVBOs();
		RayTree::Deactivate();
		// redraw
		std::cout << "about to call setup VBOs after trace photons" << std::endl;
		RayTree::setupVBOs();
		radiosity->setupVBOs();
		photon_mapping->setupVBOs();


		glutPostRedisplay();
		break; 
	}

	case 'l':	case 'L': { 
		// toggle photon rendering
		args->render_photons = !args->render_photons;
		glutPostRedisplay();
		break; }
	case 'k':	case 'K': { 
		// toggle photon rendering
		args->render_kdtree = !args->render_kdtree;
		glutPostRedisplay();
		break; }
	case 'p':	case 'P': { 
		// toggle photon rendering
		rendering_time = time(NULL);
		photon_mapping->TracePhotons();
		rendering_time = time(NULL) - rendering_time;
		std::cout << "Photon tracing completed in " << rendering_time 
			<< " seconds (" << (args->num_photons_to_shoot)/((float)rendering_time)
			<< " photons/second)\n";
		photon_mapping->setupVBOs();
		RayTree::setupVBOs();
		glutPostRedisplay();
		break; }
	case 'g':	case 'G': { 
		args->gather_indirect = true;
		args->raytracing_animation = !args->raytracing_animation;
		if (args->raytracing_animation) {
			// time the animation
			rendering_time = time(NULL);
			
			raytracing_skip = 1;//my_max(args->width,args->height) / 10;
			if (raytracing_skip % 2 == 0) raytracing_skip++;
			assert (raytracing_skip >= 1);
			raytracing_x = raytracing_skip/2;
			raytracing_y = raytracing_skip/2;
			display(); // clear out any old rendering
			printf ("photon mapping animation started, press 'G' to stop\n");
		} else {
			printf ("photon mapping animation stopped, press 'G' to start\n");		
			// Print time to render
			rendering_time =  time(NULL) - rendering_time;
			std::cout << "Rendering completed in " << rendering_time 
					<< " seconds." << std::endl;
		}
		break; 
	}
		
	case 's': case 'S':
		// subdivide the mesh for radiosity
		radiosity->Cleanup();
		radiosity->getMesh()->Subdivision();
		radiosity->Reset();
		radiosity->setupVBOs();
		glutPostRedisplay();
		break;

		// VISUALIZATIONS
	case 'w':	case 'W':
		// render wireframe mode
		args->wireframe = !args->wireframe;
		glutPostRedisplay();
		break;
	case 'v': case 'V':
		std::cout << "RENDER_MATERIALS is the only mode\n";
		break;

	case 'q':	case 'Q':
		// quit
		delete GLCanvas::photon_mapping;
		delete GLCanvas::raytracer;
		delete GLCanvas::radiosity;
		delete GLCanvas::mesh;
		exit(0);
		break;
	default:
		printf("UNKNOWN KEYBOARD INPUT	'%c'\n", key);
	}
}


// trace a ray through pixel (i,j) of the image an return the color
Vec3f GLCanvas::TraceRay(double i, double j) {
	// compute and set the pixel color
	int max_d = my_max(args->width,args->height);
	Vec3f color = Vec3f(0,0,0);
	double x = (i+0.5-args->width/2.0)/double(max_d) + 0.5;
	double y = (j+0.5-args->height/2.0)/double(max_d) + 0.5;
	Ray r = mesh->camera->generateRay(x,y);
	Hit hit = Hit();


	// ==================================
	// IMPLEMENT ANTIALIASING
	// ==================================
	double sqrt_num_antialias_samples = sqrt((double)args->num_antialias_samples);
	int stratified_num_antialais_samples = ceil(sqrt_num_antialias_samples);
	if(!args->stratified_antialiasing){
		for(int a = 0; a < args->num_antialias_samples; a++){
			double random_x = GLOBAL_mtrand.rand();
			double random_y = GLOBAL_mtrand.rand();
		//	std::cout << "x is " << random_x << " and y is " << random_y << std::endl;
			x = (i+random_x-args->width/2.0)/double(max_d) + 0.5;
			y = (j+random_y-args->height/2.0)/double(max_d) + 0.5;

			r = mesh->camera->generateRay(x,y); 
			Hit hit = Hit();
			color += raytracer->TraceRay(r,hit,args->num_bounces);
	//		std::cout << "color in for loop " << color <<std::endl;

		}
	}
	else{
		for(int x_iter = 0; x_iter < stratified_num_antialais_samples; x_iter++){
			for(int y_iter = 0; y_iter < stratified_num_antialais_samples; y_iter++){
				double random_x = GLOBAL_mtrand.rand();
				double random_y = GLOBAL_mtrand.rand();

//				std:: cout << "    Oriniginal random x:  " << random_x << std::endl;
//				std:: cout << "    Oriniginal random y:  " << random_y << std::endl;
				//scale/translate random_x and random_y to get them into the right stratified bucket
				random_x = (random_x + x)/stratified_num_antialais_samples;
				random_y = (random_y + y)/stratified_num_antialais_samples;


//				std:: cout << "    new random x:  " << random_x << std::endl;
//				std:: cout << "    new random y:  " << random_y << std::endl;
				x = (i+random_x-args->width/2.0)/double(max_d) + 0.5;
				y = (j+random_y-args->height/2.0)/double(max_d) + 0.5;
				
				r = mesh->camera->generateRay(x,y); 
				Hit hit = Hit();
				color += raytracer->TraceRay(r,hit,args->num_bounces);


			}
		}

	}
	if(!args->stratified_antialiasing){
		color /= (double)args->num_antialias_samples;
	}
	else{
		color /= (double)pow(stratified_num_antialais_samples,2); 
	}
//	std::cout << "color after loop " << color <<std::endl;

	
	// Here's what we do with a single sample per pixel:
	// construct & trace a ray through the center of the pixel
	
	
	// add that ray for visualization
	RayTree::AddMainSegment(r,0,hit.getT());



	// return the color
	return color;
}

// trace a ray through pixel (i,j) of the image an return the color
void GLCanvas::TracePhoton(double i, double j) {
	// compute and set the pixel color
	int max_d = my_max(args->width,args->height);
	
	// DO NOT DELETE (currently unused)
	// Vec3f color = Vec3f(0,0,0);
	
	double x = (i+0.5-args->width/2.0)/double(max_d) + 0.5;
	double y = (j+0.5-args->height/2.0)/double(max_d) + 0.5;
	Ray r = mesh->camera->generateRay(x,y);
	Hit hit = Hit();

	// ==================================
	// IMPLEMENT ANTIALIASING
	// ==================================
	double sqrt_num_antialias_samples = sqrt((double)args->num_antialias_samples);
	int stratified_num_antialais_samples = ceil(sqrt_num_antialias_samples);
	/*
	if(!args->stratified_antialiasing){
		for(int a = 0; a < args->num_antialias_samples; a++){
			double random_x = GLOBAL_mtrand.rand();
			double random_y = GLOBAL_mtrand.rand();

			x = (i+random_x-args->width/2.0)/double(max_d) + 0.5;
			y = (j+random_y-args->height/2.0)/double(max_d) + 0.5;

			r = mesh->camera->generateRay(x,y); 
			Hit hit = Hit();
			color += raytracer->TraceRay(r,hit,args->num_bounces);
		}
	}
*/

	for(int x_iter = 0; x_iter < stratified_num_antialais_samples; x_iter++){
		for(int y_iter = 0; y_iter < stratified_num_antialais_samples; y_iter++){
			
			double random_x = GLOBAL_mtrand.rand();
			double random_y = GLOBAL_mtrand.rand();

			random_x = (random_x + x)/stratified_num_antialais_samples;
			random_y = (random_y + y)/stratified_num_antialais_samples;

			x = (i+random_x-args->width/2.0)/double(max_d) + 0.5;
			y = (j+random_y-args->height/2.0)/double(max_d) + 0.5;
			
			r = mesh->camera->generateRay(x,y); 
			Hit hit = Hit();
			//TODO: determine what color antialiasing photons should be, if any
			photon_mapping->TracePhoton(r.getOrigin(), r.getDirection(), 580, 0, Vec4f(1.0, 1.0, 1.0, 0.7), NULL, REFRACTIVE_INDEX_OF_AIR, true);
			//TracePhoton(const Vec3f &position, const Vec3f &direction, const Vec3f &energy, int iter, Vec4f viz_color);

		}
	}

	// Here's what we do with a single sample per pixel:
	// construct & trace a ray through the center of the pixel
	
	
	// add that ray for visualization
	//RayTree::AddMainSegment(r,0,hit.getT());

}

// Scan through the image from the lower left corner across each row
// and then up to the top right.	Initially the image is sampled very
// coarsely.	Increment the static variables that track the progress
// through the scans
int GLCanvas::DrawPixel() {
	if (raytracing_x > args->width) {
		raytracing_x = raytracing_skip/2;
		raytracing_y += raytracing_skip;
	}
	if (raytracing_y > args->height) {
		if (raytracing_skip == 1)
		{
			// stop rendering, matches resolution of current camera
			rendering_time =  time(NULL) - rendering_time;
			std::cout << "Rendering completed in " << rendering_time 
					<< " seconds." << std::endl;
			return 0;
		}
		raytracing_skip = raytracing_skip / 2;
		if (raytracing_skip % 2 == 0) raytracing_skip++;
		assert (raytracing_skip >= 1);
		raytracing_x = raytracing_skip/2;
		raytracing_y = raytracing_skip/2;
		glEnd();
		glPointSize(raytracing_skip);
		glBegin(GL_POINTS);
	}

	// compute the color and position of intersection
	Vec3f color= TraceRay(raytracing_x, raytracing_y);
	double r = linear_to_srgb(color.x());
	double g = linear_to_srgb(color.y());
	double b = linear_to_srgb(color.z());
	glColor3f(r,g,b);
	//	glColor3f(1,0,0);
	double x = 2 * (raytracing_x/double(args->width)) - 1;
	double y = 2 * (raytracing_y/double(args->height)) - 1;
	glVertex3f(x,y,-1);
	raytracing_x += raytracing_skip;
	return 1;
}


void GLCanvas::idle() {
	if (args->raytracing_animation) {
		// draw 100 pixels and then refresh the screen and handle any user input
		glDisable(GL_LIGHTING);
		glDrawBuffer(GL_FRONT);
		glDisable(GL_DEPTH_TEST);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		glPointSize(raytracing_skip);
		glBegin(GL_POINTS);
		for (int i = 0; i < 400; i++) {
			if (!DrawPixel()) {
				args->raytracing_animation = false;
				
				//TODO: remove this
				//exit(0);
	break;
			}
		}
		glEnd();
		glFlush();
	}
}

// ========================================================
// ========================================================

int HandleGLError(const std::string &message) {
	GLenum error;
	int i = 0;
	while ((error = glGetError()) != GL_NO_ERROR) {
		if (message != "") {
			std::cout << "[" << message << "] ";
		}
		std::cout << "GL ERROR(" << i << ") " << gluErrorString(error) << std::endl;
		i++;
	}
	if (i == 0) return 1;
	return 0;
}

// ========================================================
// ========================================================
