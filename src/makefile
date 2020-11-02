CC = gcc
CXX = g++
DEBUG = -g
CXXFLAGS = $(shell fltk-config --use-gl --use-images --cxxflags ) -I. -I/usr/local/Cellar/glad/include
LDFLAGS = $(shell fltk-config --use-gl --use-glut --use-images --ldflags ) 
LDSTATIC = $(shell fltk-config --use-gl --use-glut --use-images --ldstaticflags ) 

TARGET = app
SRCS = CallBacks.cpp ControlPoint.cpp Track.cpp TrainView.cpp TrainWindow.cpp main.cpp
OBJS = $(SRCS:.cpp=.o)

DIR = ./Utilities
UTILITYSRCS = $(DIR)/3DUtils.cpp $(DIR)/ArcBallCam.cpp $(DIR)/Pnt3f.cpp
UTILITYOBJS=$(UTILITYSRCS:.cpp=.o)

%.o: %.cpp
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c $< -o $@

all: GLAD $(UTILITYOBJS) $(OBJS)
	$(CXX) -o $(TARGET) glad.o $(UTILITYOBJS) $(OBJS) $(LDSTATIC)

GLAD:
	$(CXX) $(CXXFLAGS) $(LDFLAGS) -c /usr/local/Cellar/glad/src/glad.c -o glad.o

.PHONY: clean
clean:
	rm -f *.o 2> /dev/null
	rm -f Utilities/*.o 2> /dev/null
	rm -f $(TARGET) 2> /dev/null
