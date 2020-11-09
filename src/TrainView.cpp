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

#define _DEBUG 0


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

	// g = new float[12];
	// memset(g, 0.0f, 12);
	// mt = new float[4];
	// memset(mt, 0.0f, 4);
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
	}
	else
		throw std::runtime_error("Could not initialize GLAD!");

	// Set up the view port
	glViewport(0,0,2*w(),2*h());

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
	// glEnable(GL_LIGHT0);
	glDisable(GL_LIGHT0);
	glDisable(GL_LIGHT1);
	glDisable(GL_LIGHT2);

	// top view only needs one light
	if (tw->topCam->value()) {
		glDisable(GL_LIGHT1);
		glDisable(GL_LIGHT2);
		glDisable(GL_LIGHT3);
		glDisable(GL_LIGHT4);
		glDisable(GL_LIGHT5);
	} else {
		// glEnable(GL_LIGHT3);
		// glEnable(GL_LIGHT4);
	}


	glDisable(GL_LIGHT3);
	glDisable(GL_LIGHT4);
	glDisable(GL_LIGHT5);

	//*********************************************************************
	//
	// * set the light parameters
	//
	//**********************************************************************

	// glLightfv(GL_LIGHT0, GL_POSITION, lightPosition1);
	// glLightfv(GL_LIGHT0, GL_DIFFUSE, whiteLight);
	// glLightfv(GL_LIGHT0, GL_AMBIENT, grayLight);

	// glLightfv(GL_LIGHT1, GL_POSITION, lightPosition2);
	// glLightfv(GL_LIGHT1, GL_DIFFUSE, yellowLight);

	// glLightfv(GL_LIGHT2, GL_POSITION, lightPosition3);
	// glLightfv(GL_LIGHT2, GL_DIFFUSE, blueLight);

	if (tw->dirLight->value()) {
		// * Directional Light
		float light3NoAmbient[4] = {.0, .0, .0, 1.0};
		float light3WhiteDiffuse[4] = {1.0, 0.0, 0.0, 1.0};
		float light3Position[4] = {1.0, 1.0, 0.0, 0.0};
		glLightfv(GL_LIGHT0, GL_AMBIENT, light3NoAmbient);
		glLightfv(GL_LIGHT0, GL_DIFFUSE, light3WhiteDiffuse);
		glLightfv(GL_LIGHT0, GL_POSITION, light3Position);
		glEnable(GL_LIGHT0);
	}
	if (tw->posLight->value()) {
		// * Point Light
		glEnable(GL_LIGHT1);
		float light4Ambient[4] = {0.1, 0.1, .0, 1.0};
		float light4Diffuse[4] = {1.0, 1.0, .0, 1.0};
		float light4Position[4] = {0.0f, 0.0f, 0.0f, 1.0f};
		glLightfv(GL_LIGHT1, GL_AMBIENT, light4Ambient);
		glLightfv(GL_LIGHT1, GL_DIFFUSE, light4Diffuse);
		glLightfv(GL_LIGHT1, GL_POSITION, light4Position);
	}
	if (tw->spotLight->value()) {
		// * Spot Light
		glEnable(GL_LIGHT2);
		float ligth5NoAmbient[4] = {1.0, 1.0, .2, 1.0};
		float light5Diffuse[4] = {1.0, 1.0, .0, 1.0};
			// train position
		float light5Position[4] = {trainParams[0].x, trainParams[0].y + 20, trainParams[0].z, 1.0f};
		glLightfv(GL_LIGHT2, GL_AMBIENT, ligth5NoAmbient);
		glLightfv(GL_LIGHT2, GL_DIFFUSE, light5Diffuse);
		glLightfv(GL_LIGHT2, GL_POSITION, light5Position);

			// train target - target position
		Pnt3f forward = trainParams[1] - trainParams[0];
		forward.normalize();
		float light5Direction[4] = {forward.x, forward.y, forward.z, 0.0};
		glLightfv(GL_LIGHT2, GL_SPOT_DIRECTION, light5Direction);
		glLightf(GL_LIGHT2, GL_SPOT_CUTOFF, 60.0f);

		glLightf(GL_LIGHT2, GL_SPOT_EXPONENT, 15.0f);

		glLightf(GL_LIGHT2, GL_CONSTANT_ATTENUATION, 0.4f);
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
		gluLookAt(trainParams[0].x, trainParams[0].y + 5, trainParams[0].z,
			trainParams[1].x, trainParams[1].y + 5, trainParams[1].z,
			trainParams[2].x, trainParams[2].y, trainParams[2].z);

#ifdef EXAMPLE_SOLUTION
		trainCamView(this,aspect);
#endif
	}
}

void TrainView::getCurve(int i, const int& size, float* g, float* m, float* t, float time, Pnt3f& qt) {
	for (int j = -1; j < 3; ++j) {
		int index = i + j;
		if (index < 0) index += size;
		else if (index >= size) index -= size;

		g[j + 1] = m_pTrack->points[index].pos.x;
		g[4 + (j + 1)] = m_pTrack->points[index].pos.y;
		g[8 + (j + 1)] = m_pTrack->points[index].pos.z;

		m[j + 1] = pow(time, 3 - (1 + j));
	}

	// multiplication
	int mSize;
	float* gm = multiply(g, m, 3, 4, 4, 4, &mSize);
	float* gmt = multiply(gm, t, 4, 4, 4, 1, &mSize);

	#ifdef _DEBUG
	std::cout << "g \n";
	for (int j = 0; j < sizeof(g)/sizeof(float); ++j) {
		std::cout << g[j] << ",";
	}
	std::cout << "end g\n";
	std::cout << "m \n";
	for (int j = 0; j < sizeof(m)/sizeof(float); ++j) {
		std::cout << m[j] << ",";
	}
	std::cout << "end m\n";
	std::cout << "t \n";
	for (int j = 0; j < sizeof(t)/sizeof(float); ++j) {
		std::cout << t[j] << ",";
	}
	std::cout << "end t\n";
	std::cout << "gmt \n";
	for (int j = 0; j < sizeof(gmt)/sizeof(float); ++j) {
		std::cout << gmt[j] << ",";
	}
	std::cout << "end gmt\n";
	#endif
	qt.x = gmt[0];
	qt.y = gmt[1];
	qt.z = gmt[2];
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
	static size_t frame = 0;
	if (!doingShadows) {
		frame += (int)tw->speed->value() <= 0 ? 1 : (int)tw->speed->value();
		if (frame >= m_pTrack->points.size() * DIVIDE_LINES) frame -= m_pTrack->points.size() * DIVIDE_LINES;
	}

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
			for (int j = -1; j < 3; ++j) {
				int index = i + j;
				if (index < 0) index += size;
				else if (index >= size) index -= size;

				g[j + 1] = m_pTrack->points[index].pos.x;
				g[4 + (j + 1)] = m_pTrack->points[index].pos.y;
				g[8 + (j + 1)] = m_pTrack->points[index].pos.z;

				mt[j + 1] = pow(t, 3 - (1 + j));
			}

			// multiplication
			int mSize;
			float* gm = multiply(g, mCardinal, 3, 4, 4, 4, &mSize);
			float* gmt = multiply(gm, mt, 4, 4, 4, 1, &mSize);
			qt.x = gmt[0];
			qt.y = gmt[1];
			qt.z = gmt[2];
		} else if (splineType == SplineCubicBSpline) {
			// getCurve((int)i, size, g, m, mt, t, qt);
			for (int j = -1; j < 3; ++j) {
				int index = i + j;
				if (index < 0) index += size;
				else if (index >= size) index -= size;

				g[j + 1] = m_pTrack->points[index].pos.x;
				g[4 + (j + 1)] = m_pTrack->points[index].pos.y;
				g[8 + (j + 1)] = m_pTrack->points[index].pos.z;

				mt[j + 1] = pow(t, 3 - (1 + j));
			}

			// multiplication
			int mSize;
			float* gm = multiply(g, mBSpline, 3, 4, 4, 4, &mSize);
			float* gmt = multiply(gm, mt, 4, 4, 4, 1, &mSize);
			qt.x = gmt[0];
			qt.y = gmt[1];
			qt.z = gmt[2];
		}

		for (size_t j = 0; j < DIVIDE_LINES; ++j) {
			// line from (qt0 -> qt1)
			qt0 = qt;
			t += percent;
			if (splineType == SplineLinear) {
				qt = (1 - t) * cpPosP1 + t * cpPosP2;
			} else if (splineType == SplineCardinalCubic) {
				for (int j = -1; j < 3; ++j) {
					int index = i + j;
					if (index < 0) index += size;
					else if (index >= size) index -= size;

					g[j + 1] = m_pTrack->points[index].pos.x;
					g[4 + (j + 1)] = m_pTrack->points[index].pos.y;
					g[8 + (j + 1)] = m_pTrack->points[index].pos.z;

					mt[j + 1] = pow(t, 3 - (1 + j));
				}

				// multiplication
				int mSize;
				float* gm = multiply(g, mCardinal, 3, 4, 4, 4, &mSize);
				float* gmt = multiply(gm, mt, 4, 4, 4, 1, &mSize);
				qt.x = gmt[0];
				qt.y = gmt[1];
				qt.z = gmt[2];
			} else if (splineType == SplineCubicBSpline) {
				for (int k = -1; k < 3; ++k) {
					int index = i + k;
					if (index < 0) index += size;
					else if (index >= size) index -= size;

					g[(k + 1)] = m_pTrack->points[index].pos.x;
					g[(k + 1) + 4] = m_pTrack->points[index].pos.y;
					g[(k + 1) + 8] = m_pTrack->points[index].pos.z;

					mt[k + 1] = pow(t, 3 - (1 + k));
				}

				// multiplication
				int gmSize;
				float* gm = multiply(g, mBSpline, 3, 4, 4, 4, &gmSize);
				int gmtSize;
				float* gmt = multiply(gm, mt, 4, 4, 4, 1, &gmtSize);

				qt.x = gmt[0];
				qt.y = gmt[1];
				qt.z = gmt[2];
			}
			qt1 = qt;

			// // draw
			// * single middle line
			// * not used anymore
			// glLineWidth(3);
			// glBegin(GL_LINES);
			// if (!doingShadows) {
			// 	glColor3ub(32, 32, 64);
			// }
			// glVertex3f(qt0.x, qt0.y, qt0.z);
			// glVertex3f(qt1.x, qt1.y, qt1.z);
			// glEnd();
			// glLineWidth(1);

			// cross
			Pnt3f orientT;
			if (splineType == SplineLinear) {
				orientT = (1 - t) * cpOrientP1 + t * cpOrientP2;
			} else if (splineType == SplineCardinalCubic) {
				for (int k = -1; k < 3; ++k) {
					int index = i + k;
					if (index < 0) index += size;
					else if (index >= size) index -= size;

					g[(k + 1)] = m_pTrack->points[index].orient.x;
					g[(k + 1) + 4] = m_pTrack->points[index].orient.y;
					g[(k + 1) + 8] = m_pTrack->points[index].orient.z;

					mt[k + 1] = pow(t, 3 - (1 + k));
				}

				// multiplication
				int gmSize;
				float* gm = multiply(g, mCardinal, 3, 4, 4, 4, &gmSize);
				int gmtSize;
				float* gmt = multiply(gm, mt, 4, 4, 4, 1, &gmtSize);

				orientT.x = gmt[0];
				orientT.y = gmt[1];
				orientT.z = gmt[2];
			} else if (splineType == SplineCubicBSpline) {
				for (int k = -1; k < 3; ++k) {
					int index = i + k;
					if (index < 0) index += size;
					else if (index >= size) index -= size;

					g[(k + 1)] = m_pTrack->points[index].orient.x;
					g[(k + 1) + 4] = m_pTrack->points[index].orient.y;
					g[(k + 1) + 8] = m_pTrack->points[index].orient.z;

					mt[k + 1] = pow(t, 3 - (1 + k));
				}

				// multiplication
				int gmSize;
				float* gm = multiply(g, mBSpline, 3, 4, 4, 4, &gmSize);
				int gmtSize;
				float* gmt = multiply(gm, mt, 4, 4, 4, 1, &gmtSize);

				orientT.x = gmt[0];
				orientT.y = gmt[1];
				orientT.z = gmt[2];
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

			if ((j / 10) % 2) {
				// * sleeps
				crossT *= 2.0f;
				glBegin(GL_QUADS);
				if (!doingShadows) {
					glColor3ub(255, 255, 255);
				}
				// line1
				glVertex3f(qt0.x + crossT.x, qt0.y + crossT.y, qt0.z + crossT.z);
				glVertex3f(qt1.x + crossT.x, qt1.y + crossT.y, qt1.z + crossT.z);

				// line2
				glVertex3f(qt1.x - crossT.x, qt1.y - crossT.y, qt1.z - crossT.z);
				glVertex3f(qt0.x - crossT.x, qt0.y - crossT.y, qt0.z - crossT.z);
				glEnd();
			}

			if (frame == (i * DIVIDE_LINES) + j) {
				if (!tw->trainCam->value()) {
					// glBegin(GL_QUADS);
					// if (!doingShadows) {
					// 	glColor3ub(255, 255, 255);	
					// }
					// glTexCoord2f(0.0f, 0.0f);
					// glVertex3f(qt.x - 5, qt.y - 5, qt.z - 5);
					// glTexCoord2f(1.0f, 0.0f);
					// glVertex3f(qt.x + 5, qt.y - 5, qt.z - 5);
					// glTexCoord2f(1.0f, 1.0f);
					// glVertex3f(qt.x + 5, qt.y + 5, qt.z - 5);
					// glTexCoord2f(0.0f, 1.0f);
					// glVertex3f(qt.x - 5, qt.y + 5, qt.z - 5);
					// glEnd();

					Pnt3f u = (qt1 - qt0); u.normalize();
					Pnt3f w = u * orientT; w.normalize();
					Pnt3f v = w * u; v.normalize();

					float rotation[16] = {
						u.x, u.y, u.z, 0.0,
						v.x, v.y, v.z, 0.0,
						w.x, w.y, w.z, 0.0,
						0.0, 0.0, 0.0, 1.0
					};

					glPushMatrix();
					glTranslatef(qt.x, qt.y, qt.z);
					glMultMatrixf(rotation);
					glScalef(5.0f, 5.0f, 5.0f);
					glTranslatef(0.0f, 1.0f, 0.0f);
					draw_cube(doingShadows);
					glPopMatrix();
				}
				// setup train params
				trainParams[0] = qt0; trainParams[1] = qt1; trainParams[2] = orientT;
			}
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

void TrainView::draw_cube(bool doingShadows) {
	glBegin(GL_QUADS);                // Begin drawing the color cube with 6 quads
	// Top face (y = 1.0f)
	// Define vertices in counter-clockwise (CCW) order with normal pointing out
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(0, 255, 0);     // Green
	glNormal3f(0.0f, 1.0f, 0.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);

	// Bottom face (y = -1.0f)
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(255, 128, 0);     // Orange
	glNormal3f(0.0f, -1.0f, 0.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);

	// Front face  (z = 1.0f)
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(255, 0, 0);     // Red
	glNormal3f(0.0f, 0.0f, 1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);

	// Back face (z = -1.0f)
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(255, 255, 0);     // Yellow
	glNormal3f(0.0f, 0.0f, -1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);

	// Left face (x = -1.0f)
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(0, 0, 255);     // Blue
	glNormal3f(-1.0f, 0.0f, 0.0f);
	glVertex3f(-1.0f, 1.0f, 1.0f);
	glVertex3f(-1.0f, 1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, -1.0f);
	glVertex3f(-1.0f, -1.0f, 1.0f);

	// Right face (x = 1.0f)
	if (!doingShadows)
		glColor3ub(255, 255, 255);
		// glColor3ub(255, 0, 255);     // Magenta
	glNormal3f(1.0f, 0.0f, 0.0f);
	glVertex3f(1.0f, 1.0f, -1.0f);
	glVertex3f(1.0f, 1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, 1.0f);
	glVertex3f(1.0f, -1.0f, -1.0f);
	glEnd();  // End of drawing color-cube
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
	int mx = 2 * Fl::event_x(); 
	int my = 2 * Fl::event_y();

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

	printf("Selected Cube %d\n",selectedCube);
}