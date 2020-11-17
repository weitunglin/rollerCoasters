/************************************************************************
     File:        TrainView.H

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

#pragma once
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include "objloader.hpp"

// Preclarify for preventing the compiler error
class TrainWindow;
class CTrack;
class Pnt3f;


//#######################################################################
// TODO
// You might change the TrainView in order to add different objects to
// be drawn, or adjust the lighting, or ...
//#######################################################################


// The TrainView also needs to handle its user interface events.
// since this is an FLTk Window, we need to include FlTk Headers
// the pragmas turn off the warnings from FlTk - our code should have no warnings,
// but their code doesn't meet that standard
#pragma warning(push)
#pragma warning(disable:4312)
#pragma warning(disable:4311)
#include <Fl/Fl_Gl_Window.h>
#pragma warning(pop)

// this uses the old ArcBall Code
#include "Utilities/ArcBallCam.H"

#define DIVIDE_LINES 300
#define SLEEPS_WIDTH 15

typedef enum {
	SplineLinear = 0,
	SplineCardinalCubic = 1,
	SplineCubicBSpline = 2
} SplineT;

class TrainPosition {
	public:
		size_t frame;
		double location;

		TrainPosition() : frame(0), location(0) {

		}
};

class TrackTrainInfo {
	public:
		std::vector<Pnt3f> trainParams;
		double placements;
		size_t frameNum;

		TrackTrainInfo() : trainParams(std::vector<Pnt3f>(4)), placements(0), frameNum(0) {

		}

		friend std::ostream& operator<<(std::ostream& os, const TrackTrainInfo& t);
};

class TrainView : public Fl_Gl_Window
{
	public:
		// note that we keep the "standard widget" constructor arguments
		TrainView(int x, int y, int w, int h, const char* l = 0);
		~TrainView() {
			delete[] trainParams;
		}

		// overrides of important window things
		virtual int handle(int);
		virtual void draw();

		// all of the actual drawing happens in this routine
		// it has to be encapsulated, since we draw differently if
		// we're drawing shadows (no colors, for example)
		void drawStuff(bool doingShadows=false);

		// setup the projection - assuming that the projection stack has been
		// cleared for you
		void setProjection();

		// Reset the Arc ball control
		void resetArcball();

		// pick a point (for when the mouse goes down)
		void doPick();

		void drawCube(bool doingShadows);

		void drawCar(bool doingShadows);

		float* multiply(float* m, float* n, int m1, int m2, int n1, int n2, int* size);

		void getCurve(int i, const int& size, float* g, float* m, float* t, float time, Pnt3f& qt, int type);

		double getDistance(const Pnt3f& a, const Pnt3f& b);

		size_t getFrameNum(double location) const;

	public:
		ArcBallCam arcball;			// keep an ArcBall for the UI
		int selectedCube;  			// simple - just remember which cube is selected

		TrainWindow* tw;			// The parent of this display window
		CTrack*	m_pTrack;			// The track of the entire scene

		Pnt3f* trainParams;

		std::vector<TrackTrainInfo> frameTable;

		TrainPosition* trainPosition;

		std::vector<glm::vec3> trainVertices;
		std::vector<glm::vec2> trainUvs;
		std::vector<glm::vec3> trainNormals;

		std::vector<glm::vec3> carVertices;
		std::vector<glm::vec2> carUvs;
		std::vector<glm::vec3> carNormals;

		const std::string TRAIN_OBJ_NAME = "model.obj";
		const std::string CAR_OBJ_NAME = "car.obj";

		uint32_t trainVertexBuffer, trainUvBuffer, trainNormalBuffer;
		uint32_t carVertexBuffer, carUvBuffer, carNormalBuffer;
		bool inited = false;
};
