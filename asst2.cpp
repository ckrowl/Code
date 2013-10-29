////////////////////////////////////////////////////////////////////////
//
//   Harvard University
//   CS175 : Computer Graphics
//   Professor Steven Gortler
//
////////////////////////////////////////////////////////////////////////
/****************************************************** 
* Project:         CS 116A Homework #2 
* File:            asst2.cpp 
* Purpose:         Translations, rotations, and scalings in OpenGL
* Start date:      9/28/13 
* Programmer:      Zane Melcho 
* 
****************************************************** 
*/

#include <vector>
#include <string>
#include <memory>
#include <stdexcept>
#include <stack>
#if __GNUG__
#   include <tr1/memory>
#endif

#include <GL/glew.h>
#ifdef __MAC__
#   include <GLUT/glut.h>
#else
#   include <glut.h>
#endif


#include "cvec.h"
#include "matrix4.h"
#include "geometrymaker.h"
#include "ppm.h"
#include "glsupport.h"
#include "rigtform.h"

#define M_PI 3.1415926535897932384626433832795;


using namespace std;      // for string, vector, iostream, and other standard C++ stuff
using namespace tr1; // for shared_ptr

// G L O B A L S ///////////////////////////////////////////////////

// --------- IMPORTANT --------------------------------------------------------
// Before you start working on this assignment, set the following variable
// properly to indicate whether you want to use OpenGL 2.x with GLSL 1.0 or
// OpenGL 3.x+ with GLSL 1.3.
//
// Set g_Gl2Compatible = true to use GLSL 1.0 and g_Gl2Compatible = false to
// use GLSL 1.3. Make sure that your machine supports the version of GLSL you
// are using. In particular, on Mac OS X currently there is no way of using
// OpenGL 3.x with GLSL 1.3 when GLUT is used.
//
// If g_Gl2Compatible=true, shaders with -gl2 suffix will be loaded.
// If g_Gl2Compatible=false, shaders with -gl3 suffix will be loaded.
// To complete the assignment you only need to edit the shader files that get
// loaded
// ----------------------------------------------------------------------------
static const bool g_Gl2Compatible = false;

struct RigidBody;

static const float g_frustMinFov = 60.0;  // A minimal of 60 degree field of view
static float g_frustFovY = g_frustMinFov; // FOV in y direction (updated by updateFrustFovY)

static const float g_frustNear = -0.1;    // near plane
static const float g_frustFar = -50.0;    // far plane
static const float g_groundY = 0.0;      // y coordinate of the ground
static const float g_groundSize = 15.0;   // half the ground length

static int g_windowWidth = 512;
static int g_windowHeight = 512;
static bool g_mouseClickDown = false;    // is the mouse button pressed
static bool g_mouseLClickButton, g_mouseRClickButton, g_mouseMClickButton;
static int g_mouseClickX, g_mouseClickY; // coordinates for mouse click event
static int g_activeShader = 0;
static const int g_numOfObjects = 2; //Number of cube objects to be drawn

struct ShaderState {
  GlProgram program;

  // Handles to uniform variables
  GLint h_uLight, h_uLight2;
  GLint h_uProjMatrix;
  GLint h_uModelViewMatrix;
  GLint h_uNormalMatrix;
  GLint h_uColor;

  // Handles to vertex attributes
  GLint h_aPosition;
  GLint h_aNormal;

  ShaderState(const char* vsfn, const char* fsfn) {
    readAndCompileShader(program, vsfn, fsfn);

    const GLuint h = program; // short hand

    // Retrieve handles to uniform variables
    h_uLight = safe_glGetUniformLocation(h, "uLight");
    h_uLight2 = safe_glGetUniformLocation(h, "uLight2");
    h_uProjMatrix = safe_glGetUniformLocation(h, "uProjMatrix");
    h_uModelViewMatrix = safe_glGetUniformLocation(h, "uModelViewMatrix");
    h_uNormalMatrix = safe_glGetUniformLocation(h, "uNormalMatrix");
    h_uColor = safe_glGetUniformLocation(h, "uColor");

    // Retrieve handles to vertex attributes
    h_aPosition = safe_glGetAttribLocation(h, "aPosition");
    h_aNormal = safe_glGetAttribLocation(h, "aNormal");

    if (!g_Gl2Compatible)
      glBindFragDataLocation(h, 0, "fragColor");
    checkGlErrors();
  }

};

static const int g_numShaders = 2;
static const char * const g_shaderFiles[g_numShaders][2] = {
  {"./shaders/basic-gl3.vshader", "./shaders/diffuse-gl3.fshader"},
  {"./shaders/basic-gl3.vshader", "./shaders/solid-gl3.fshader"}
};
static const char * const g_shaderFilesGl2[g_numShaders][2] = {
  {"./shaders/basic-gl2.vshader", "./shaders/diffuse-gl2.fshader"},
  {"./shaders/basic-gl2.vshader", "./shaders/solid-gl2.fshader"}
};
static vector<shared_ptr<ShaderState> > g_shaderStates; // our global shader states

// --------- Geometry

// Macro used to obtain relative offset of a field within a struct
#define FIELD_OFFSET(StructType, field) &(((StructType *)0)->field)

// A vertex with floating point position and normal
struct VertexPN {
  Cvec3f p, n;

  VertexPN() {}
  VertexPN(float x, float y, float z,
           float nx, float ny, float nz)
    : p(x,y,z), n(nx, ny, nz)
  {}

  // Define copy constructor and assignment operator from GenericVertex so we can
  // use make* functions from geometrymaker.h
  VertexPN(const GenericVertex& v) {
    *this = v;
  }

  VertexPN& operator = (const GenericVertex& v) {
    p = v.pos;
    n = v.normal;
    return *this;
  }
};

struct Geometry {
  GlBufferObject vbo, ibo;
  int vboLen, iboLen;

  Geometry(VertexPN *vtx, unsigned short *idx, int vboLen, int iboLen) {
    this->vboLen = vboLen;
    this->iboLen = iboLen;

    // Now create the VBO and IBO
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, sizeof(VertexPN) * vboLen, vtx, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned short) * iboLen, idx, GL_STATIC_DRAW);
  }

  void draw(const ShaderState& curSS) {
    // Enable the attributes used by our shader
    safe_glEnableVertexAttribArray(curSS.h_aPosition);
    safe_glEnableVertexAttribArray(curSS.h_aNormal);

    // bind vbo
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    safe_glVertexAttribPointer(curSS.h_aPosition, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, p));
    safe_glVertexAttribPointer(curSS.h_aNormal, 3, GL_FLOAT, GL_FALSE, sizeof(VertexPN), FIELD_OFFSET(VertexPN, n));

    // bind ibo
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);

    // draw!
    glDrawElements(GL_TRIANGLES, iboLen, GL_UNSIGNED_SHORT, 0);

    // Disable the attributes used by our shader
    safe_glDisableVertexAttribArray(curSS.h_aPosition);
    safe_glDisableVertexAttribArray(curSS.h_aNormal);
  }

	void draw(const ShaderState& curSS, Matrix4 MVM)
	{
		Matrix4 NMVM = normalMatrix(MVM);

		GLfloat glmatrix[16];
		MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
		safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

		NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
		safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);

		draw(curSS);
	}
};
/*-----------------------------------------------*/
#include "RigidBody.h"
/*-----------------------------------------------*/
// Vertex buffer and index buffer associated with the ground and cube geometry
static shared_ptr<Geometry> g_ground, g_cube, g_cylinder;

// --------- Scene

static const Cvec3 g_light1(2.0, 3.0, 14.0), g_light2(-2, -3.0, -5.0);  // define two lights positions in world space
static RigTForm g_skyRbt = RigTForm(Cvec3(0.0, 3, 10.0)); // Default camera
static RigTForm g_eyeRbt = g_skyRbt; //Set the g_eyeRbt frame to be default as the sky frame
//static Matrix4 g_objectRbt[g_numOfObjects];// = {Matrix4::makeTranslation(Cvec3(0,0,0))};  // currently only 1 obj is defined
//Different color for each body part for debug.
//static Cvec3f g_objectColors[g_numOfObjects] = {Cvec3f(0,0,0), Cvec3f(0, 1, 0), Cvec3f(0,0,0), Cvec3f(1,1,1), Cvec3f(0,1,0), Cvec3f(0,0,0), Cvec3f(0,0,0), Cvec3f(.5,.5,.5), Cvec3f(0,0,0), Cvec3f(1,0,0)}; // Array of colors to fill cubes with
static RigidBody g_rigidBodies[g_numOfObjects]; // Array that holds each cube(body part) of the robot

///////////////// END OF G L O B A L S //////////////////////////////////////////////////
/*-----------------------------------------------*/
static Matrix4 lookAt(Cvec3f eyePosition, Cvec3f lookAtPosition, Cvec3f upVector)
{
	Cvec3f x, y, z, w;
	double m[16];

	//Different from the book but works correctly
	z = normalize(eyePosition - lookAtPosition);
	x = normalize(cross(upVector,z));
	y = cross(z,x);	

	int k = 0;

	for (int i = 0; i < 3; i++)
	{
		m[k] = x[i];
		k++;
		m[k] = y[i];
		k++;
		m[k] = z[i];
		k++;
		m[k] = eyePosition[i];
		k++;
	}

	m[12] = 0;
	m[13] = 0;
	m[14] = 0;
	m[15] = 1;

	//return Matrix4();

	return Matrix4(m, true);
}
/*-----------------------------------------------*/
static void initCamera()
{
	Cvec3f eye = Cvec3f(0.0, 3.0, 10.0);
	Cvec3f at = Cvec3f(0.0, 0.0, 0.0);
	Cvec3f up = Cvec3f(0.0,1.0,0.0);
	//g_skyRbt = lookAt(eye, at, up); // Default camera
	//g_eyeRbt = g_skyRbt;
}
/*-----------------------------------------------*/
static void initGround() {
  // A x-z plane at y = g_groundY of dimension [-g_groundSize, g_groundSize]^2
  VertexPN vtx[4] = {
    VertexPN(-g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
    VertexPN(-g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY,  g_groundSize, 0, 1, 0),
    VertexPN( g_groundSize, g_groundY, -g_groundSize, 0, 1, 0),
  };
  unsigned short idx[] = {0, 1, 2, 0, 2, 3};
  g_ground.reset(new Geometry(&vtx[0], &idx[0], 4, 6));
}

static Geometry* initCubes() 
{
  int ibLen, vbLen;
  getCubeVbIbLen(vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCube(1, vtx.begin(), idx.begin());
  return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}

static Geometry* initCylinder() {
  int ibLen, vbLen;
  int slices = 12;
  int stacks = 12;
  float radius = 10.0;
  float height = 1;
  getCylinderVbIbLen(slices, vbLen, ibLen);

  // Temporary storage for cube geometry
  vector<VertexPN> vtx(vbLen);
  vector<unsigned short> idx(ibLen);

  makeCylinder(radius, slices, height, vtx.begin(), idx.begin());
  return new Geometry(&vtx[0], &idx[0], vbLen, ibLen);
}
/*-----------------------------------------------*/
void initRobot()
{
	/* PURPOSE:		Initializes each body part of the robot to it's position, scale, and rotation  
	*/

	//0-Robot - Used as a container for the robot as a whole
	RigTForm rigTemp = RigTForm(Cvec3(0,4.388,0), Quat().makeYRotation(90));
	Matrix4 scaleTemp = Matrix4::makeScale(Cvec3(2,2,2));//1.9
	
	RigidBody *robot = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,0,1));
	robot->name = "robot";
	robot->isVisible = false;

	//torso->num_children = 4;
   //torso->children = new RigidBody*[torso->num_children];

	//1-Head
	rigTemp = RigTForm(Cvec3(0.0,0.0,0.0));//1.85
	scaleTemp = Matrix4::makeScale(Cvec3(1,1,1));
	RigidBody *head = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,1,0));
	head->name = "head";

	//2-Body
	rigTemp = RigTForm(Cvec3(0.0,-1.6,0.0));//-1.6
	scaleTemp = Matrix4::makeScale(Cvec3(1.6,2.4,1.2));
	RigidBody *body = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,0,0));
	body->name = "body";

	//3-Left Arm
	rigTemp = RigTForm(Cvec3(-0.95,0.65,0.6), Quat().makeXRotation(90));
	scaleTemp = Matrix4::makeScale(Cvec3(0.25,0.6,0.25));
	RigidBody *leftArm = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,1,0));
	leftArm->name = "leftArm";

	//4-Right Arm
	rigTemp = RigTForm(Cvec3(.95,0.65,0.6), Quat().makeXRotation(90));
	scaleTemp = Matrix4::makeScale(Cvec3(0.25,0.6,.25));
	RigidBody *rightArm = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,1,0));
	rightArm->name = "rightArm";

	//5-Left Leg
	rigTemp = RigTForm(Cvec3(-0.25,-0.85,0.0));
	scaleTemp = Matrix4::makeScale(Cvec3(0.25,0.8,0.25));
	RigidBody *leftLeg = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,0,0));
	leftLeg->name = "leftLeg";

	//6-Right Leg
	rigTemp = RigTForm(Cvec3(0.25,-0.85,0.0));
	scaleTemp = Matrix4::makeScale(Cvec3(0.25,0.8,0.25));
	RigidBody *rightLeg = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,0,0));
	rightLeg->name = "rightLeg";

	//7-Bolt
	rigTemp = RigTForm(Cvec3(0.0,0.0,0.0));
	scaleTemp = Matrix4::makeScale(Cvec3(1.2,.1,0.1));
	RigidBody *bolt = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0.5,0.5,0.5));
	bolt->name = "bolt";

	//8-Hair
	rigTemp = RigTForm(Cvec3(0.0,.5,0.0));
	scaleTemp = Matrix4::makeScale(Cvec3(1.2,.2,1.2));
	RigidBody *hair = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(0,0,0));
	hair->name = "hair";

	//9-ScaleReference
	rigTemp = RigTForm(Cvec3(0.0,1.26,0.0));
	scaleTemp = Matrix4::makeScale(Cvec3(0.1,5,0.1));
	RigidBody *scaleReference = new RigidBody(rigTemp, scaleTemp, NULL, initCubes(), Cvec3(1,0,0));
	scaleReference->name = "scaleReference";

	//10-Cylinder
	//rigTemp = RigTForm(Cvec3(-2.0,3.0,0.0));
	//scaleTemp = Matrix4::makeScale(Cvec3(.1,.3,.1));
	//RigidBody *cylinder = RigidBody(rigTemp, scaleTemp, &RigidBody(), initCylinder(), Cvec3(1,1,1));
	//cylinder->name = "cylinder";

	//Setup Children
	robot->numOfChildren = 1;
	head->numOfChildren = 3;
	body->numOfChildren = 4;

	robot->children = new RigidBody*[robot->numOfChildren];
	robot->children[0] = head;

	head->children = new RigidBody*[head->numOfChildren];
	head->children[0] = hair;
	head->children[1] = bolt;
	head->children[2] = body;

	body->children = new RigidBody*[body->numOfChildren];
	body->children[0]  = leftArm;
	body->children[1]  = rightArm;
	body->children[2]  = leftLeg;
	body->children[3]  = rightLeg;
	
	//Add to global RigidBody array
	g_rigidBodies[0] = *robot;
	g_rigidBodies[1] = *scaleReference;
	//g_rigidBodies[2] = cylinder;

	glutPostRedisplay();
}
/*-----------------------------------------------*/
// takes a projection matrix and send to the the shaders
static void sendProjectionMatrix(const ShaderState& curSS, const Matrix4& projMatrix) {
  GLfloat glmatrix[16];
  projMatrix.writeToColumnMajorMatrix(glmatrix); // send projection matrix
  safe_glUniformMatrix4fv(curSS.h_uProjMatrix, glmatrix);
}

// takes MVM and its normal matrix to the shaders
static void sendModelViewNormalMatrix(const ShaderState& curSS, const Matrix4& MVM, const Matrix4& NMVM) {
  GLfloat glmatrix[16];
  MVM.writeToColumnMajorMatrix(glmatrix); // send MVM
  safe_glUniformMatrix4fv(curSS.h_uModelViewMatrix, glmatrix);

  NMVM.writeToColumnMajorMatrix(glmatrix); // send NMVM
  safe_glUniformMatrix4fv(curSS.h_uNormalMatrix, glmatrix);
}

// update g_frustFovY from g_frustMinFov, g_windowWidth, and g_windowHeight
static void updateFrustFovY() {
  if (g_windowWidth >= g_windowHeight)
    g_frustFovY = g_frustMinFov;
  else {
    const double RAD_PER_DEG = 0.5 * CS175_PI/180;
    g_frustFovY = atan2(sin(g_frustMinFov * RAD_PER_DEG) * g_windowHeight / g_windowWidth, cos(g_frustMinFov * RAD_PER_DEG)) / RAD_PER_DEG;
  }
}

static Matrix4 makeProjectionMatrix() {
  return Matrix4::makeProjection(
           g_frustFovY, g_windowWidth / static_cast <double> (g_windowHeight),
           g_frustNear, g_frustFar);
}
/*-----------------------------------------------*/
static void drawStuff() 
{
	/* PURPOSE:		Draws objects in relative 3d space  
	*/

	// short hand for current shader state
	const ShaderState& curSS = *g_shaderStates[g_activeShader];

	// build & send proj. matrix to vshader
	const Matrix4 projmat = makeProjectionMatrix();
	sendProjectionMatrix(curSS, projmat);

	// Use the g_eyeRbt as the eyeRbt;
	const RigTForm eyeRbt = g_eyeRbt;
	const RigTForm invEyeRbt = inv(eyeRbt);

	const Cvec3 eyeLight1 = Cvec3(invEyeRbt * Cvec4(g_light1, 1)); // g_light1 position in eye coordinates
	const Cvec3 eyeLight2 = Cvec3(invEyeRbt * Cvec4(g_light2, 1)); // g_light2 position in eye coordinates
	safe_glUniform3f(curSS.h_uLight, eyeLight1[0], eyeLight1[1], eyeLight1[2]);
	safe_glUniform3f(curSS.h_uLight2, eyeLight2[0], eyeLight2[1], eyeLight2[2]);

	// draw ground
	// ===========
	//
	const RigTForm groundRbt = RigTForm();  // identity
	Matrix4 MVM = RigTForm::makeTRmatrix(invEyeRbt * groundRbt,Matrix4::makeScale(Cvec3(1,1,1)));
	Matrix4 NMVM = normalMatrix(MVM);
	sendModelViewNormalMatrix(curSS, MVM, NMVM);
	safe_glUniform3f(curSS.h_uColor, 0.1, 0.95, 0.1); // set color
	g_ground->draw(curSS);

	// Draw all Rigid body objects
	for (int i = 0; i < g_numOfObjects; i++)
		g_rigidBodies[i].drawRigidBody(curSS, invEyeRbt);
}
/*-----------------------------------------------*/
static void display() {
  glUseProgram(g_shaderStates[g_activeShader]->program);
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);                   // clear framebuffer color&depth

  drawStuff();

  glutSwapBuffers();                                    // show the back buffer (where we rendered stuff)

  checkGlErrors();
}

static void reshape(const int w, const int h) {
  g_windowWidth = w;
  g_windowHeight = h;
  glViewport(0, 0, w, h);
  cerr << "Size of window is now " << w << "x" << h << endl;
  updateFrustFovY();
  glutPostRedisplay();
}

static void motion(const int x, const int y) {
	const double dx = x - g_mouseClickX;
	const double dy = g_windowHeight - y - 1 - g_mouseClickY;

	RigTForm m;
	if (g_mouseLClickButton && !g_mouseRClickButton) { // left button down?
		m = g_rigidBodies[0].rtf * RigTForm(Quat().makeXRotation(-dy)) * RigTForm(Quat().makeYRotation(dx)) * inv(g_rigidBodies[0].rtf);
	}
  else if (g_mouseRClickButton && !g_mouseLClickButton) { // right button down?
   m = g_eyeRbt * RigTForm(Cvec3(dx, 0, -dy) * 0.01) * inv(g_eyeRbt); //Update based on Eye Frame
	//m = Matrix4::makeTranslation(Cvec3(dx, dy, 0) * 0.01);
  }
  else if (g_mouseMClickButton || (g_mouseLClickButton && g_mouseRClickButton)) {  // middle or (left and right) button down?
    //m = g_eyeRbt * Matrix4::makeTranslation(Cvec3(0, 0, -dy) * 0.01) * inv(g_eyeRbt); //Update based on Eye Frame
	  m = RigTForm(Cvec3(0, -dy, 0) * 0.01); //Update based on Eye Frame
  }

  if (g_mouseClickDown) {
//	  g_objectRbt[0] *= m; // Simply right-multiply is WRONG
	  g_rigidBodies[0].rtf = g_rigidBodies[0].rtf * m; //Update Robot Container
	  glutPostRedisplay(); // we always redraw if we changed the scene
  }

  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;
}


static void mouse(const int button, const int state, const int x, const int y) {
  g_mouseClickX = x;
  g_mouseClickY = g_windowHeight - y - 1;  // conversion from GLUT window-coordinate-system to OpenGL window-coordinate-system

  g_mouseLClickButton |= (button == GLUT_LEFT_BUTTON && state == GLUT_DOWN);
  g_mouseRClickButton |= (button == GLUT_RIGHT_BUTTON && state == GLUT_DOWN);
  g_mouseMClickButton |= (button == GLUT_MIDDLE_BUTTON && state == GLUT_DOWN);

  g_mouseLClickButton &= !(button == GLUT_LEFT_BUTTON && state == GLUT_UP);
  g_mouseRClickButton &= !(button == GLUT_RIGHT_BUTTON && state == GLUT_UP);
  g_mouseMClickButton &= !(button == GLUT_MIDDLE_BUTTON && state == GLUT_UP);

  g_mouseClickDown = g_mouseLClickButton || g_mouseRClickButton || g_mouseMClickButton;
}
/*-----------------------------------------------*/
static void keyUp (const unsigned char key, const int x, const int y)
{
	/* PURPOSE:		An OpenGL callback when a keyboard key is released 
		RECEIVES:	unsigned char key - Key released
						int x - location of the mouse x_position when released
						int y - location of the mouse y_position when released
	*/

	if (key == '1' || key == '2' || key == '3')
	{
		//Reset to the default Eye view
		g_eyeRbt = g_skyRbt;
	}
	//Reset Scaling or Translation of Robot back to original from hop or scale
	if (key == 's' || key == 'h' || key == 'g')
	{
		//g_rigidBodies[0].originalPosition();
	}
	//Reset Torso back to original rotation from Twist
	if (key == 't' || key == 'y')
	{
		//g_rigidBodies[2].originalPosition();
	}
	//Reset Robot, arms, and legs back to original position from Jumping Jack
	if (key == 'j')
	{
		//g_rigidBodies[0].originalPosition();
		//for (int i = 3; i < g_numOfObjects; i++)
			//g_rigidBodies[i].originalPosition();
	}

	glutPostRedisplay();
}
/*-----------------------------------------------*/
static void keyboard(const unsigned char key, const int x, const int y) 
{
	/* PURPOSE:		OpenGL Callback for Keyboard presses 
		RECEIVES:	unsigned char key - Key pressed
						int x - Mouse x_position when key pressed
						int y - Mouse y_position when key pressed
		REMARKS:		Handles robot modifications based on key presses and then requests a redisplay
	*/

	switch (key) 
	{
		case 27:
			exit(0);                                  // ESC
		case 'i':
			cout << " ============== H E L P ==============\n\n"
			<< "i\t\thelp menu\n"
			<< "s\t\tsave screenshot\n"
			<< "f\t\tToggle flat shading on/off.\n"
			<< "o\t\tCycle object to edit\n"
			<< "v\t\tCycle view\n"
			<< "drag left mouse to rotate\n" << endl;
			break;
		case 's':
			glFlush();
			writePpmScreenshot(g_windowWidth, g_windowHeight, "out.ppm");
			break;
		case 'f':
			g_activeShader ^= 1;
			break;
  }

	if (key == '1')
	{
		/*
		//Set eye view to look down at Jumping Jack from above
		g_eyeRbt = Matrix4::makeTranslation(Cvec3(0.0, 10.0, 0.0));
		g_eyeRbt = g_eyeRbt * Matrix4::makeXRotation(-90);
		*/
	}
	else if (key == '2')
	{
		/*
		//Set Eye view to look at Jumping Jacks Right Side
		g_eyeRbt = Matrix4::makeTranslation(Cvec3(10.0, 2.5, 0.0));
		g_eyeRbt = g_eyeRbt * Matrix4::makeYRotation(90);
		*/
	}
	else if (key == '3')
	{
		/*
		//Set Eye view to look at Jumping Jacks Left Side
		g_eyeRbt = Matrix4::makeTranslation(Cvec3(-10.0, 2.5, 0.0));
		g_eyeRbt = g_eyeRbt * Matrix4::makeYRotation(-90);
		*/
	}
	
	if (key == 'j')
	{
		/*
		//Jumping jack with arms and legs outstretched

		//Robot
		g_rigidBodies[0].resetTransRbt();
		g_rigidBodies[0].setTransRbt(Matrix4::makeTranslation(Cvec3(0,-0.2,0)));
		//Left Arm
		g_rigidBodies[3].resetRotateRbt();
		g_rigidBodies[3].setRotateRbt(Matrix4::makeZRotation(-135));
		g_rigidBodies[3].resetTransRbt();
		g_rigidBodies[3].setTransRbt(Matrix4::makeTranslation(Cvec3(-0.75,0.5,0)));
		//Right Arm
		g_rigidBodies[4].resetRotateRbt();
		g_rigidBodies[4].setRotateRbt(Matrix4::makeZRotation(135));
		g_rigidBodies[4].resetTransRbt();
		g_rigidBodies[4].setTransRbt(Matrix4::makeTranslation(Cvec3(0.75,0.5,0)));
		//Left Leg
		g_rigidBodies[5].resetRotateRbt();
		g_rigidBodies[5].setRotateRbt(Matrix4::makeZRotation(-45));
		g_rigidBodies[5].resetTransRbt();
		Matrix4 lArmTransRbt = g_rigidBodies[5].getTransRbt();
		g_rigidBodies[5].setTransRbt(lArmTransRbt * Matrix4::makeTranslation(Cvec3(-0.2,0.1,0)));
		//Right Leg
		g_rigidBodies[6].resetRotateRbt();
		g_rigidBodies[6].setRotateRbt(Matrix4::makeZRotation(45));
		g_rigidBodies[6].resetTransRbt();
		Matrix4 rArmTransRbt = g_rigidBodies[6].getTransRbt();
		g_rigidBodies[6].setTransRbt(rArmTransRbt * Matrix4::makeTranslation(Cvec3(0.2,0.1,0)));
		*/
	}
	else if (key == 't')
	{
		/*
		//Torso twist by 45 degrees 
		g_rigidBodies[2].resetRotateRbt();
		g_rigidBodies[2].setRotateRbt(Matrix4::makeYRotation(45));
		*/
	}
	else if (key == 'y')
	{
		/*
		//Torso twist by -45 degrees 
		g_rigidBodies[2].resetRotateRbt();
		g_rigidBodies[2].setRotateRbt(Matrix4::makeYRotation(-45));
		*/
	}
	else if (key == 's')
	{
		/*
		//Supersize to twice the normal size
		g_rigidBodies[0].resetScaleRbt();
		g_rigidBodies[0].setScaleRbt(Matrix4::makeScale(Cvec3(2,2,2)));
		*/
	}
	else if (key == 'h')
	{
		/*
		//Hop a fixed displacement to the right side
		g_rigidBodies[0].resetTransRbt();
		g_rigidBodies[0].setTransRbt(Matrix4::makeTranslation(Cvec3(1,0,0)));
		*/
	}
	else if (key == 'g')
	{
		/*
		//Hop a fixed displacement to the left side
		g_rigidBodies[0].resetTransRbt();
		g_rigidBodies[0].setTransRbt(Matrix4::makeTranslation(Cvec3(-1,0,0)));
		*/
	}
	glutPostRedisplay();
}
/*-----------------------------------------------*/
static void initGlutState(int argc, char * argv[]) {
  glutInit(&argc, argv);                                  // initialize Glut based on cmd-line args
  glutInitDisplayMode(GLUT_RGBA|GLUT_DOUBLE|GLUT_DEPTH);  //  RGBA pixel channels and double buffering
  glutInitWindowSize(g_windowWidth, g_windowHeight);      // create a window
  glutCreateWindow("Assignment 2");                       // title the window

  glutDisplayFunc(display);                               // display rendering callback
  glutReshapeFunc(reshape);                               // window reshape callback
  glutMotionFunc(motion);                                 // mouse movement callback
  glutMouseFunc(mouse);                                   // mouse click callback
  glutKeyboardFunc(keyboard);
  glutKeyboardUpFunc(keyUp);										// Keyboard key released
}

static void initGLState() {
  glClearColor(128./255., 200./255., 255./255., 0.);
  glClearDepth(0.);
  glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
  glPixelStorei(GL_PACK_ALIGNMENT, 1);
  glCullFace(GL_BACK);
  glEnable(GL_CULL_FACE);
  glEnable(GL_DEPTH_TEST);
  glDepthFunc(GL_GREATER);
  glReadBuffer(GL_BACK);
  if (!g_Gl2Compatible)
    glEnable(GL_FRAMEBUFFER_SRGB);
}

static void initShaders() {
  g_shaderStates.resize(g_numShaders);
  for (int i = 0; i < g_numShaders; ++i) {
    if (g_Gl2Compatible)
      g_shaderStates[i].reset(new ShaderState(g_shaderFilesGl2[i][0], g_shaderFilesGl2[i][1]));
    else
      g_shaderStates[i].reset(new ShaderState(g_shaderFiles[i][0], g_shaderFiles[i][1]));
  }
}

static void initGeometry() 
{
	//Initialize Object Matrix array
	//initObjectRbt();
	initRobot();
	initGround();
	//initCubes();
	//initCylinder();
}

int main(int argc, char * argv[]) {
  try {
		initGlutState(argc,argv);

		glewInit(); // load the OpenGL extensions

		cout << (g_Gl2Compatible ? "Will use OpenGL 2.x / GLSL 1.0" : "Will use OpenGL 3.x / GLSL 1.3") << endl;
		if ((!g_Gl2Compatible) && !GLEW_VERSION_3_0)
		throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.3");
		else if (g_Gl2Compatible && !GLEW_VERSION_2_0)
		throw runtime_error("Error: card/driver does not support OpenGL Shading Language v1.0");

		initGLState();
		initShaders();
		initGeometry();
		initCamera();

/*		
		//Debug stuff
		cout << "\n";
		Matrix4::print(g_skyRbt);
		cout << "\n";
*/
		glutMainLoop();
		return 0;
  }
  catch (const runtime_error& e) {
    cout << "Exception caught: " << e.what() << endl;
    return -1;
  }
}
