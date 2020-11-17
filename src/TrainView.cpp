/************************************************************************
     File:        TrainView.cpp

     Author:     
                  Michael Gleicher, gleicher@cs.wisc.edu

     Modifier
                  Yu-Chi Lai, yu-chi@cs.wisc.edu
     
     Comment:     
						The TrainView is the window that actually shows the 
						train. Its a
						GL display canvas (Fl_Gl_Window).  It is held within 
						a TrainWindow
						that is the outer window with all the widgets. 
						The TrainView needs 
						to be aware of the window - since it might need to 
						check the widgets to see how to draw

	  Note:        we need to have pointers to this, but maybe not know 
						about it (beware circular references)

     Platform:    Visio Studio.Net 2003/2005

*************************************************************************/

#include <Fl/fl.h>

// we will need OpenGL, and OpenGL needs windows.h
// #include <windows.h>
//#include "GL/gl.h"
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <FL/glu.h>

#include "TrainView.H"
#include "TrainWindow.H"
#include "Utilities/3DUtils.H"


#ifdef EXAMPLE_SOLUTION
#	include "TrainExample/TrainExample.H"
#endif

//************************************************************************
//
// * Constructor to set up the GL window
//========================================================================
TrainView::
TrainView(int x, int y, int w, int h, const char* l) 
	: Fl_Gl_Window(x,y,w,h,l)
//========================================================================
{
	mode( FL_RGB|FL_ALPHA|FL_DOUBLE | FL_STENCIL );

	resetArcball();
	
	trainParams = new Pnt3f[3];
	for (int i = 0; i < 3; ++i) {
		trainParams[i] = Pnt3f();
	}
	
	trainPosition = new TrainPosition();

	bool res = loadOBJ(TRAIN_OBJ_NAME.c_str(), trainVertices, trainUvs, trainNormals);
	if (!res) {
		throw std::runtime_error("Could not open train obj file");
	} else {
		std::cout << "load train obj file finished\n";
	}
	res = loadOBJ(CAR_OBJ_NAME.c_str(), carVertices, carUvs, carNormals);
	if (!res) {
		throw std::runtime_error("Could not open car obj file");
	} else {
		std::cout << "load car obj file finished\n";
	}
}

//************************************************************************
//
// * Reset the camera to look at the world
//========================================================================
void TrainView::
resetArcball()
//========================================================================
{
	// Set up the camera to look at the world
	// these parameters might seem magical, and they kindof are
	// a little trial and error goes a long way
	arcball.setup(this, 40, 250, .2f, .4f, 0);
}

//************************************************************************
//
// * FlTk Event handler for the window
//########################################################################
// TODO: 
//       if you want to make the train respond to other events 
//       (like key presses), you might want to hack this.
//########################################################################
//========================================================================
int TrainView::handle(int event)
{
	// see if the ArcBall will handle the event - if it does, 
	// then we're done
	// note: the arcball only gets the event if we're in world view
	if (tw->worldCam->value())
		if (arcball.handle(event)) 
			return 1;

	// remember what button was used
	static int last_push;

	switch(event) {
		// Mouse button being pushed event
		case FL_PUSH:
			last_push = Fl::event_button();
			// if the left button be pushed is left mouse button
			if (last_push == FL_LEFT_MOUSE  ) {
				doPick();
				damage(1);
				return 1;
			};
			break;

	   // Mouse button release event
		case FL_RELEASE: // button release
			damage(1);
			last_push = 0;
			return 1;

		// Mouse button drag event
		case FL_DRAG:

			// Compute the new control point position
			if ((last_push == FL_LEFT_MOUSE) && (selectedCube >= 0)) {
				ControlPoint* cp = &m_pTrack->points[selectedCube];

				double r1x, r1y, r1z, r2x, r2y, r2z;
				getMouseLine(r1x, r1y, r1z, r2x, r2y, r2z);

				double rx, ry, rz;
				mousePoleGo(r1x, r1y, r1z, r2x, r2y, r2z, 
								static_cast<double>(cp->pos.x), 
								static_cast<double>(cp->pos.y),
								static_cast<double>(cp->pos.z),
								rx, ry, rz,
								(Fl::event_state() & FL_CTRL) != 0);

				cp->pos.x = (float) rx;
				cp->pos.y = (float) ry;
				cp->pos.z = (float) rz;
				damage(1);
			}
			break;

		// in order to get keyboard events, we need to accept focus
		case FL_FOCUS:
			return 1;

		// every time the mouse enters this window, aggressively take focus
		case FL_ENTER:	
			focus(this);
			break;

		case FL_KEYBOARD:
		 		int k = Fl::event_key();
				int ks = Fl::event_state();
				if (k == 'p') {
					// Print out the selected control point information
					if (selectedCube >= 0) 
						printf("Selected(%d) (%g %g %g) (%g %g %g)\n",
								 selectedCube,
								 m_pTrack->points[selectedCube].pos.x,
								 m_pTrack->points[selectedCube].pos.y,
								 m_pTrack->points[selectedCube].pos.z,
								 m_pTrack->points[selectedCube].orient.x,
								 m_pTrack->points[selectedCube].orient.y,
								 m_pTrack->points[selectedCube].orient.z);
					else
						printf("Nothing Selected\n");

					return 1;
				};
				break;
	}

	return Fl_Gl_Window::handle(event);
}

//************************************************************************
//
// * this is the code that actually draws the window
//   it puts a lot of the work into other routines to simplify things
//========================================================================
void TrainView::draw()
{

	//*********************************************************************
	//
	// * Set up basic opengl informaiton
	//
	//**********************************************************************
	//initialized glad
	if (gladLoadGL())
	{
		//initiailize VAO, VBO, Shader...
		if (!inited) {
			glGenBuffers(1, &trainVertexBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, trainVertexBuffer);
			glBufferData(GL_ARRAY_BUFFER, trainVertices.size() * sizeof(glm::vec3), &trainVertices[0], GL_STATIC_DRAW);

			glGenBuffers(1, &trainUvBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, trainUvBuffer);
			glBufferData(GL_ARRAY_BUFFER, trainUvs.size() * sizeof(glm::vec2), &trainUvs[0], GL_STATIC_DRAW);

			glGenBuffers(1, &trainNormalBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, trainNormalBuffer);
			glBufferData(GL_ARRAY_BUFFER, trainNormals.size() * sizeof(glm::vec3), &trainNormals[0], GL_STATIC_DRAW);
			std::cout << (int)trainVertexBuffer << "," << (int)trainUvBuffer << "," << (int)trainNormalBuffer << std::endl;
			
			glGenBuffers(1, &carVertexBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, carVertexBuffer);
			glBufferData(GL_ARRAY_BUFFER, carVertices.size() * sizeof(glm::vec3), &carVertices[0], GL_STATIC_DRAW);

			glGenBuffers(1, &carUvBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, carUvBuffer);
			glBufferData(GL_ARRAY_BUFFER, carUvs.size() * sizeof(glm::vec2), &carUvs[0], GL_STATIC_DRAW);

			glGenBuffers(1, &carNormalBuffer);
			glBindBuffer(GL_ARRAY_BUFFER, carNormalBuffer);
			glBufferData(GL_ARRAY_BUFFER, carNormals.size() * sizeof(glm::vec3), &carNormals[0], GL_STATIC_DRAW);
			std::cout << (int)carVertexBuffer << "," << (int)carUvBuffer << "," << (int)carNormalBuffer << std::endl;

			inited = true;
		}
	}
	else
		throw std::runtime_error("Could not initialize GLAD!");

	// Set up the view port
	#ifdef RETINA
	glViewport(0,0,2*w(),2*h());
	#else
	glViewport(0,0,w(),h());
	#endif

	// clear the window, be sure to clear the Z-Buffer too
	glClearColor(0,0,.3f,0);		// background should be blue

	// we need to clear out the stencil buffer since we'll use
	// it for shadows
	glClearStencil(0);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT | GL_STENCIL_BUFFER_BIT);
	glEnable(GL_DEPTH);

	// Blayne prefers GL_DIFFUSE
    glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

	// prepare for projection
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity();
	setProjection();		// put the code to set up matrices here

	//######################################################################
	// TODO: 
	// you might want to set the lighting up differently. if you do, 
	// we need to set up the lights AFTER setting up the projection
	//######################################################################
	// enable the lighting
	glEnable(GL_COLOR_MATERIAL);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_LIGHTING);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);
	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);
	glDisable(GL_LIGHT5);

	// top view only needs one light
	if (tw->topCam->value()) {
		glDisable(GL_LIGHT0);
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHT2);
		glDisable(GL_LIGHT4);
		glDisable(GL_LIGHT5);
	} else {
		// glEnable(GL_LIGHT3);
		// glEnable(GL_LIGHT4);
	}



	//*********************************************************************
	//
	// * set the light parameters
	//
	//**********************************************************************
	GLfloat lightPosition1[]	= {0,1,1,0}; // {50, 200.0, 50, 1.0};
	GLfloat lightPosition2[]	= {1, 0, 0, 0};
	GLfloat lightPosition3[]	= {0, -1, 0, 0};
	GLfloat yellowLight[]		= {0.5f, 0.5f, .1f, 1.0};
	GLfloat whiteLight[]		= {1.0f, 1.0f, 1.0f, 1.0};
	GLfloat blueLight[]			= {.1f,.1f,.3f,1.0};
	GLfloat grayLight[]			= {.3f, .3f, .3f, 1.0};

	if (tw->baseLight->value()) {
		glEnable(GL_LIGHT3);
		glLightfv(GL_LIGHT3, GL_DIFFUSE, whiteLight);
		glLightfv(GL_LIGHT3, GL_AMBIENT, blueLight);
		glLightfv(GL_LIGHT3, GL_POSITION, lightPosition1);
	}

	// glLightfv(GL_LIGHT1, GL_POSITION, lightPosition2);
	// glLightfv(GL_LIGHT1, GL_DIFFUSE, yellowLight);

	// glLightfv(GL_LIGHT2, GL_POSITION, lightPosition3);
	// glLightfv(GL_LIGHT2, GL_DIFFUSE, blueLight);

	if (tw->dirLight->value()) {
		// * Directional Light
		glEnable(GL_LIGHT0);
		float light3NoAmbient[4] = {.1, .1, .1, 1.0};
		float light3WhiteDiffuse[4] = {1.0, .0, 0.0, 1.0}; // yellow
		float light3Position[4] = {1.0, 0.0, 0.0, 0.0};
		glLightfv(GL_LIGHT0, GL_AMBIENT, light3NoAmbient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light3WhiteDiffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, light3Position);
	}
	else if (tw->posLight->value()) {
		// * Point Light
		glEnable(GL_LIGHT1);
		float light4Ambient[4] = {.1, .1, .1, 1.0};
		float light4Diffuse[4] = {1.0, .0, .0, 1.0};
		float light4Position[4] = {0.0f, 10.0f, 0.0f, 1.0f};
		glLightfv(GL_LIGHT1, GL_AMBIENT, light4Ambient);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light4Diffuse);
		glLightfv(GL_LIGHT1, GL_POSITION, light4Position);
	}
	else if (tw->spotLight->value()) {
		// * Spot Light
		glEnable(GL_LIGHT2);
		float ligth5NoAmbient[4] = {.1, .1, .1, 1.0};
		float light5Diffuse[4] = {1.0, 1.0, .0, 1.0};
			// train position
		float light5Position[4] = {frameTable[trainPosition->frame].trainParams[0].x, frameTable[trainPosition->frame].trainParams[0].y + 20, frameTable[trainPosition->frame].trainParams[0].z, 1.0f};
		glLightfv(GL_LIGHT2, GL_AMBIENT, ligth5NoAmbient);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, light5Diffuse);
		glLightfv(GL_LIGHT2, GL_POSITION, light5Position);

			// train target - target position
		Pnt3f forward = frameTable[trainPosition->frame].trainParams[1] - frameTable[trainPosition->frame].trainParams[0];
		forward.normalize();
		float light5Direction[4] = {forward.x, forward.y, forward.z, 0.0};
		glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light5Direction);
		glLightf(GL_LIGHT2, GL_SPOT_CUTOFF, 60.0f);

		glLightf(GL_LIGHT2, GL_SPOT_EXPONENT, 15.0f);

		glLightf(GL_LIGHT2, GL_CONSTANT_ATTENUATION, 0.2f);
		glLightf(GL_LIGHT2, GL_LINEAR_ATTENUATION, .0f);
		glLightf(GL_LIGHT2, GL_QUADRATIC_ATTENUATION, .0f);
	}

	//*********************************************************************
	// now draw the ground plane
	//*********************************************************************
	// set to opengl fixed pipeline(use opengl 1.x draw function)
	glUseProgram(0);

	setupFloor();
	if (!tw->spotLight->value())
		glDisable(GL_LIGHTING);
	drawFloor(200,10);


	//*********************************************************************
	// now draw the object and we need to do it twice
	// once for real, and then once for shadows
	//*********************************************************************
	glEnable(GL_LIGHTING);
	setupObjects();

	drawStuff();

	// this time drawing is for shadows (except for top view)
	if (!tw->topCam->value()) {
		setupShadows();
		drawStuff(true);
		unsetupShadows();
	}
}

//************************************************************************
//
// * This sets up both the Projection and the ModelView matrices
//   HOWEVER: it doesn't clear the projection first (the caller handles
//   that) - its important for picking
//========================================================================
void TrainView::
setProjection()
//========================================================================
{
	// Compute the aspect ratio (we'll need it)
	float aspect = static_cast<float>(w()) / static_cast<float>(h());

	// Check whether we use the world camp
	if (tw->worldCam->value())
		arcball.setProjection(false);
	// Or we use the top cam
	else if (tw->topCam->value()) {
		float wi, he;
		if (aspect >= 1) {
			wi = 110;
			he = wi / aspect;
		} 
		else {
			he = 110;
			wi = he * aspect;
		}

		// Set up the top camera drop mode to be orthogonal and set
		// up proper projection matrix
		glMatrixMode(GL_PROJECTION);
		glOrtho(-wi, wi, -he, he, 200, -200);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		glRotatef(-90,1,0,0);
	} 
	// Or do the train view or other view here
	//####################################################################
	// TODO: 
	// put code for train view projection here!	
	//####################################################################
	else {
		glMatrixMode(GL_PROJECTION);
		glLoadIdentity();
		gluPerspective(120, aspect, 1, 200);
		glMatrixMode(GL_MODELVIEW);
		glLoadIdentity();
		gluLookAt(frameTable[trainPosition->frame].trainParams[0].x, frameTable[trainPosition->frame].trainParams[0].y + 5, frameTable[trainPosition->frame].trainParams[0].z,
			frameTable[trainPosition->frame].trainParams[1].x, frameTable[trainPosition->frame].trainParams[1].y + 5, frameTable[trainPosition->frame].trainParams[1].z,
			frameTable[trainPosition->frame].trainParams[2].x, frameTable[trainPosition->frame].trainParams[2].y, frameTable[trainPosition->frame].trainParams[2].z);

#ifdef EXAMPLE_SOLUTION
		trainCamView(this,aspect);
#endif
	}
}

void TrainView::getCurve(int i, const int& size, float* g, float* m, float* t, float time, Pnt3f& qt, int type) {
	for (int j = -1; j < 3; ++j) {
		int index = i + j;
		if (index < 0) index += size;
		else if (index >= size) index -= size;

		if (type == 0) {
			g[j + 1] = m_pTrack->points[index].pos.x;
			g[4 + (j + 1)] = m_pTrack->points[index].pos.y;
			g[8 + (j + 1)] = m_pTrack->points[index].pos.z;
		} else if (type == 1) {
			g[j + 1] = m_pTrack->points[index].orient.x;
			g[4 + (j + 1)] = m_pTrack->points[index].orient.y;
			g[8 + (j + 1)] = m_pTrack->points[index].orient.z;
		}

		t[j + 1] = pow(time, 3 - (1 + j));
	}

	// multiplication
	int mSize;
	float* gm = multiply(g, m, 3, 4, 4, 4, &mSize);
	float* gmt = multiply(gm, t, 4, 4, 4, 1, &mSize);
	qt.x = gmt[0];
	qt.y = gmt[1];
	qt.z = gmt[2];
	free(gm);
	free(gmt);
}

size_t TrainView::getFrameNum(double location) const {
	for (size_t i = 0; i < frameTable.size(); ++i) {
		if (location <= frameTable[i].placements) {
			return i;
		}
	}
	return frameTable.size() - 1;
}

float* TrainView::multiply(float* m, float* n, int m1, int m2, int n1, int n2, int* size) {
	float* result = (float*) malloc(sizeof(float) * (m1 * n2));
	for (int i = 0; i < m1; ++i) {
		for (int j = 0; j < n2; ++j) {
			float sum = 0.0;
			for (int k = 0; k < m2; ++k) {
				sum += m[(i * m2) + k] * n[(n2 * k) + j];
			}
			result[i * n2 + j] = sum;
		}
	}

	*size = m1 * n2;
	return result;
}

double TrainView::getDistance(const Pnt3f& a, const Pnt3f& b) {
	return pow(pow(a.x - b.x, 2) + pow(a.y - b.y, 2) + pow(a.z - b.z, 2), 0.5);
}

//************************************************************************
//
// * this draws all of the stuff in the world
//
//	NOTE: if you're drawing shadows, DO NOT set colors (otherwise, you get 
//       colored shadows). this gets called twice per draw 
//       -- once for the objects, once for the shadows
//########################################################################
// TODO: 
// if you have other objects in the world, make sure to draw them
//########################################################################
//========================================================================
void TrainView::drawStuff(bool doingShadows)
{
	// Draw the control points
	// don't draw the control points if you're driving 
	// (otherwise you get sea-sick as you drive through them)
	if (!tw->trainCam->value()) {
		for(size_t i=0; i<m_pTrack->points.size(); ++i) {
			if (!doingShadows) {
				if ( ((int) i) != selectedCube)
					glColor3ub(240, 60, 60);
				else
					glColor3ub(240, 240, 30);
			}
			m_pTrack->points[i].draw();
		}
	}

	// draw the track
	//####################################################################
	// TODO: 
	// call your own track drawing code
	//####################################################################

#ifdef EXAMPLE_SOLUTION
	drawTrack(this, doingShadows);
#endif

	int splineType = tw->splineBrowser->value() - 1;

	static float g[12];
	static float mt[4];
	static float mBSpline[16] = {
		-1/6.0, 3/6.0, -3/6.0, 1/6.0,
		3/6.0, -6/6.0, 0.0, 4/6.0,
		-3/6.0, 3/6.0, 3/6.0, 1/6.0,
		1/6.0, 0.0, 0.0, 0.0
	};
	static float mCardinal[16] = {
		-1/2.0, 2/2.0, -1/2.0, 0.0,
		3/2.0, -5/2.0, 0.0, 2/2.0,
		-3/2.0, 4/2.0, 1/2.0, 0.0,
		1/2.0, -1/2.0, 0.0, 0.0
	};

	Pnt3f qt, qt0, qt1;
	int size = m_pTrack->points.size();
	double sumOfPlacements = 0;
	frameTable = std::vector<TrackTrainInfo>(size * DIVIDE_LINES);
	for (size_t i = 0; i < size; ++i) {
		// pos
		Pnt3f cpPosP1 = m_pTrack->points[i].pos;
		Pnt3f cpPosP2 = m_pTrack->points[(i+1) % (m_pTrack->points.size())].pos;
		
		// orient
		Pnt3f cpOrientP1 = m_pTrack->points[i].orient;
		Pnt3f cpOrientP2 = m_pTrack->points[(i+1) % (m_pTrack->points.size())].orient;

		float percent = 1.0f / DIVIDE_LINES;
		float t = 0;

		if (splineType == SplineLinear) {
			qt = (1 - t) * cpPosP1 + t * cpPosP2;
		} else if (splineType == SplineCardinalCubic) {
			getCurve(i, size, g, mCardinal, mt, t, qt, 0);
		} else if (splineType == SplineCubicBSpline) {
			getCurve(i, size, g, mBSpline, mt, t, qt, 0);
		}

		for (size_t j = 0; j < DIVIDE_LINES; ++j) {
			// line from (qt0 -> qt1)
			qt0 = qt;
			t += percent;
			if (splineType == SplineLinear) {
				qt = (1 - t) * cpPosP1 + t * cpPosP2;
			} else if (splineType == SplineCardinalCubic) {
				getCurve(i, size, g, mCardinal, mt, t, qt, 0);
			} else if (splineType == SplineCubicBSpline) {
				getCurve(i, size, g, mBSpline, mt, t, qt, 0);
			}
			qt1 = qt;

			// cross
			Pnt3f orientT;
			if (splineType == SplineLinear) {
				orientT = (1 - t) * cpOrientP1 + t * cpOrientP2;
			} else if (splineType == SplineCardinalCubic) {
				getCurve(i, size, g, mCardinal, mt, t, orientT, 1);
			} else if (splineType == SplineCubicBSpline) {
				getCurve(i, size, g, mBSpline, mt, t, orientT, 1);
			}
			orientT.normalize();
			Pnt3f crossT = (qt1 - qt0) * orientT;
			crossT.normalize();
			crossT *= 2.5f;

			// * lines with cross
			// draw
			glLineWidth(3);
			glBegin(GL_LINES);
			if (!doingShadows) {
				glColor3ub(32, 32, 64);
			}
			// line1
			glVertex3f(qt0.x + crossT.x, qt0.y + crossT.y, qt0.z + crossT.z);
			glVertex3f(qt1.x + crossT.x, qt1.y + crossT.y, qt1.z + crossT.z);
			// line2
			glVertex3f(qt0.x - crossT.x, qt0.y - crossT.y, qt0.z - crossT.z);
			glVertex3f(qt1.x - crossT.x, qt1.y - crossT.y, qt1.z - crossT.z);
			glEnd();
			glLineWidth(1);

			if ((j / SLEEPS_WIDTH) % 2) {
				// * sleeps
				crossT *= 2.0f;
				glBegin(GL_QUADS);
				if (!doingShadows) {
					glColor3ub(106, 73, 64);
				}

				glNormal3f(0.0f, 0.0f, 1.0f);
				// line1
				glVertex3f(qt0.x + crossT.x, qt0.y + crossT.y, qt0.z + crossT.z);
				glVertex3f(qt1.x + crossT.x, qt1.y + crossT.y, qt1.z + crossT.z);
				// line2
				glVertex3f(qt1.x - crossT.x, qt1.y - crossT.y, qt1.z - crossT.z);
				glVertex3f(qt0.x - crossT.x, qt0.y - crossT.y, qt0.z - crossT.z);
				glEnd();
			}

			frameTable[i * DIVIDE_LINES + j].placements = sumOfPlacements;
			frameTable[i * DIVIDE_LINES + j].frameNum = i * DIVIDE_LINES + j;
			frameTable[i * DIVIDE_LINES + j].trainParams[0] = qt0;
			frameTable[i * DIVIDE_LINES + j].trainParams[1] = qt1;
			frameTable[i * DIVIDE_LINES + j].trainParams[2] = orientT;
			frameTable[i * DIVIDE_LINES + j].trainParams[3] = qt;

			sumOfPlacements += getDistance(qt1, qt0);
		}
	}

#ifdef _DEBUG
	if (!doingShadows) {
		for (size_t i = 0; i < frameTable.size(); ++i) {
			std::cout << frameTable[i] << " ";
		}
		std::cout << "\n\n";
	}
#endif

	if (!doingShadows && !tw->trainCam->value()) {
		Pnt3f u = (frameTable[trainPosition->frame].trainParams[1] - frameTable[trainPosition->frame].trainParams[0]); u.normalize();
		Pnt3f w = u * frameTable[trainPosition->frame].trainParams[2]; w.normalize();
		Pnt3f v = w * u; v.normalize();

		float rotation[16] = {
			u.x, u.y, u.z, 0.0,
			v.x, v.y, v.z, 0.0,
			w.x, w.y, w.z, 0.0,
			0.0, 0.0, 0.0, 1.0
		};

		glPushMatrix();
		glTranslatef(frameTable[trainPosition->frame].trainParams[3].x, frameTable[trainPosition->frame].trainParams[3].y, frameTable[trainPosition->frame].trainParams[3].z);
		glMultMatrixf(rotation);
		glScalef(5.0f, 5.0f, 5.0f);
		glTranslatef(0.0f, 0.0f, 0.0f);
		drawCube(doingShadows);
		glPopMatrix();

		for (size_t i = 0; i < tw->carNum; ++i) {
			double carLocation = trainPosition->location - (30 + (i) * 30);
			if (carLocation < 0) carLocation += frameTable.back().placements;
			size_t carFrame = getFrameNum(carLocation);

			u = (frameTable[carFrame].trainParams[1] - frameTable[carFrame].trainParams[0]); u.normalize();
			w = u * frameTable[carFrame].trainParams[2]; w.normalize();
			v = w * u; v.normalize();

			rotation[0] = u.x; rotation[1] = u.y; rotation[2] = u.z; rotation[3] = 0.0f;
			rotation[4] = v.x; rotation[5] = v.y; rotation[6] = v.z; rotation[7] = 0.0f;
			rotation[8] = w.x; rotation[9] = w.y; rotation[10] = w.z; rotation[11] = 0.0f;
			rotation[12] = 0.0f; rotation[13] = 0.0f; rotation[14] = 0.0f; rotation[15] = 1.0f;

			glPushMatrix();
			glTranslatef(frameTable[carFrame].trainParams[3].x, frameTable[carFrame].trainParams[3].y, frameTable[carFrame].trainParams[3].z);
			glMultMatrixf(rotation);
			glScalef(5.0f, 5.0f, 5.0f);
			glTranslatef(0.0f, 0.0f, 0.0f);
			drawCar(doingShadows);
			glPopMatrix();
		}
	}

	if (!doingShadows && tw->runButton->value()) {
		if (tw->arcLength->value()) {
			trainPosition->location += tw->speed->value() / 2.5;
			if (trainPosition->location > frameTable.back().placements) {
				trainPosition->location = 0;
			}
			trainPosition->frame = getFrameNum(trainPosition->location);
		} else {
			trainPosition->frame += (int)tw->speed->value() <= 0 ? 1 : (int)tw->speed->value();
			if (trainPosition->frame >= m_pTrack->points.size() * DIVIDE_LINES) trainPosition->frame -= m_pTrack->points.size() * DIVIDE_LINES;
			trainPosition->location = frameTable[trainPosition->frame].placements;
		}
	}

	// draw the train
	//####################################################################
	// TODO: 
	//	call your own train drawing code
	//####################################################################
#ifdef EXAMPLE_SOLUTION
	// don't draw the train if you're looking out the front window
	if (!tw->trainCam->value())
		drawTrain(this, doingShadows);
#endif
}

void TrainView::drawCube(bool doingShadows) {
	if (!inited) return;
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < trainVertices.size(); ++i) {
		glNormal3f(trainNormals[i].x, trainNormals[i].y, trainNormals[i].z);
		glVertex3f(trainVertices[i].x, trainVertices[i].y, trainVertices[i].z);
	}
	glEnd();

	// glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
	// // Top face (y = 1.0f)
	// // Define vertices in counter-clockwise (CCW) order with normal pointing out
	// glNormal3f(0.0f, 1.0f, 0.0f);
	// if (!doingShadows)
	// 	glColor3ub(0, 255, 0);     // Green
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(1.0f, 1.0f, -1.0f);
	// glVertex3f(-1.0f, 1.0f, -1.0f);
	// glVertex3f(-1.0f, 1.0f, 1.0f);
	// glVertex3f(1.0f, 1.0f, 1.0f);

	// // Bottom face (y = -1.0f)
	// glNormal3f(0.0f, -1.0f, 0.0f);
	// if (!doingShadows)
	// 	glColor3ub(255, 128, 0);     // Orange
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(1.0f, -1.0f, 1.0f);
	// glVertex3f(-1.0f, -1.0f, 1.0f);
	// glVertex3f(-1.0f, -1.0f, -1.0f);
	// glVertex3f(1.0f, -1.0f, -1.0f);

	// // Front face  (z = 1.0f)
	// glNormal3f(0.0f, 0.0f, 1.0f);
	// if (!doingShadows)
	// 	glColor3ub(255, 0, 0);     // Red
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(1.0f, 1.0f, 1.0f);
	// glVertex3f(-1.0f, 1.0f, 1.0f);
	// glVertex3f(-1.0f, -1.0f, 1.0f);
	// glVertex3f(1.0f, -1.0f, 1.0f);

	// // Back face (z = -1.0f)
	// glNormal3f(0.0f, 0.0f, -1.0f);
	// if (!doingShadows)
	// 	glColor3ub(255, 255, 0);     // Yellow
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(1.0f, -1.0f, -1.0f);
	// glVertex3f(-1.0f, -1.0f, -1.0f);
	// glVertex3f(-1.0f, 1.0f, -1.0f);
	// glVertex3f(1.0f, 1.0f, -1.0f);

	// // Left face (x = -1.0f)
	// glNormal3f(-1.0f, 0.0f, 0.0f);
	// if (!doingShadows)
	// 	glColor3ub(0, 0, 255);     // Blue
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(-1.0f, 1.0f, 1.0f);
	// glVertex3f(-1.0f, 1.0f, -1.0f);
	// glVertex3f(-1.0f, -1.0f, -1.0f);
	// glVertex3f(-1.0f, -1.0f, 1.0f);

	// // Right face (x = 1.0f)
	// glNormal3f(1.0f, 0.0f, 0.0f);
	// if (!doingShadows)
	// 	glColor3ub(255, 0, 255);     // Magenta
	// 	// glColor3ub(255, 255, 255);
	// glVertex3f(1.0f, 1.0f, -1.0f);
	// glVertex3f(1.0f, 1.0f, 1.0f);
	// glVertex3f(1.0f, -1.0f, 1.0f);
	// glVertex3f(1.0f, -1.0f, -1.0f);
	// glEnd();  // End of drawing color-cube
}

void TrainView::drawCar(bool doingShawdows) {
	if (!inited) return;
	glBegin(GL_TRIANGLES);
	for (size_t i = 0; i < carVertices.size(); ++i) {
		glNormal3f(carNormals[i].x, carNormals[i].y, carNormals[i].z);
		glVertex3f(carVertices[i].x, carVertices[i].y, carVertices[i].z);
	}
	glEnd();
}

// 
//************************************************************************
//
// * this tries to see which control point is under the mouse
//	  (for when the mouse is clicked)
//		it uses OpenGL picking - which is always a trick
//########################################################################
// TODO: 
//		if you want to pick things other than control points, or you
//		changed how control points are drawn, you might need to change this
//########################################################################
//========================================================================
void TrainView::
doPick()
//========================================================================
{
	// since we'll need to do some GL stuff so we make this window as 
	// active window
	make_current();		

	// where is the mouse?
	#ifdef RETINA
	int mx = 2 * Fl::event_x(); 
	int my = 2 * Fl::event_y();
	#else
	int mx = Fl::event_x(); 
	int my = Fl::event_y();
	#endif

	// get the viewport - most reliable way to turn mouse coords into GL coords
	int viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);

	// Set up the pick matrix on the stack - remember, FlTk is
	// upside down!
	glMatrixMode(GL_PROJECTION);
	glLoadIdentity ();
	gluPickMatrix((double)mx, (double)(viewport[3]-my), 
						5, 5, viewport);

	// now set up the projection
	setProjection();

	// now draw the objects - but really only see what we hit
	GLuint buf[100];
	glSelectBuffer(100,buf);
	glRenderMode(GL_SELECT);
	glInitNames();
	glPushName(0);

	// draw the cubes, loading the names as we go
	for(size_t i=0; i<m_pTrack->points.size(); ++i) {
		glLoadName((GLuint) (i+1));
		m_pTrack->points[i].draw();
	}

	// go back to drawing mode, and see how picking did
	int hits = glRenderMode(GL_RENDER);
	if (hits) {
		// warning; this just grabs the first object hit - if there
		// are multiple objects, you really want to pick the closest
		// one - see the OpenGL manual 
		// remember: we load names that are one more than the index
		selectedCube = buf[3]-1;
	} else // nothing hit, nothing selected
		selectedCube = -1;

	// printf("Selected Cube %d\n",selectedCube);
}

std::ostream& operator<<(std::ostream& os, const TrackTrainInfo& t) {
	os << "#" << t.frameNum << ", " << t.placements << "\n";
	os << t.trainParams[0].print() << "\n"
		<< t.trainParams[1].print() << "\n"
		<< t.trainParams[2].print() << "\n"
		<< t.trainParams[3].print() << "\n";
	return os;
}
