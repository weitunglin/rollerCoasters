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

// Preclarify for preventing the compiler error
class TrainWindow;
class CTrack;


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

#define DIVIDE_LINES 200

typedef enum {
	SplineLinear = 0,
	SplineCardinalCubic = 1,
	SplineCubicBSpline = 2
} SplineT;

class TrainView : public Fl_Gl_Window
{
	public:
		// note that we keep the "standard widget" constructor arguments
		TrainView(int x, int y, int w, int h, const char* l = 0);

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

		float* multiply(float* m, float* n, int m1, int m2, int n1, int n2, int* size) {
			float* result = (float*) malloc(sizeof(float) * (m1 * n2));

			// for (int i = 0; i < m1; ++i) {
			// 	for (int j = 0; j < m2; ++j) {
			// 		std::cout << m[i*m1 + j] << " ";
			// 	}
			// 	std::cout << "\n";
			// }
			// std::cout << " X " << "\n";
			// for (int i = 0; i < n1; ++i) {
			// 	for (int j = 0; j < n2; ++j) {
			// 		std::cout << n[i*n1 + j] << " ";
			// 	}
			// 	std::cout << "\n";
			// }
			// std::cout << " = " << "\n";

			for (int i = 0; i < m1; ++i) {
				for (int j = 0; j < n2; ++j) {
					float sum = 0.0;
					for (int k = 0; k < m2; ++k) {
						sum += m[(i * m2) + k] * n[(n2 * k) + j];
					}
					result[i * n2 + j] = sum;
					// std::cout << result[i * n2 + j] << " ";
				}
				// std::cout << "\n";
			}
			// std::cout << "\n";

			*size = m1 * n2;
			return result;
		}

	public:
		ArcBallCam		arcball;			// keep an ArcBall for the UI
		int				selectedCube;  // simple - just remember which cube is selected

		TrainWindow*	tw;				// The parent of this display window
		CTrack*			m_pTrack;		// The track of the entire scene

		Pnt3f lookatParam[3];
};
