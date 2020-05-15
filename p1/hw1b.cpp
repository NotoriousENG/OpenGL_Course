// Thank YOU to my amazing Computer Graphics Prof and/or TA for this link
// https://www.cise.ufl.edu/class/cap4730sp20/pdf/cap4730/1bfaq.html
// Iwould not have gotten anywhere without it
// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <sstream>
#include <math.h>
#include <iostream>
#include <algorithm> // std::swap
// Include GLEW
#include <GL/glew.h>
// Include GLFW
#include <glfw/glfw3.h>
// Include GLM
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
using namespace glm;
// Include AntTweakBar
#include <AntTweakBar.h>

#include <common/shader.hpp>
#include <common/controls.hpp>
#include <common/objloader.hpp>
#include <common/vboindexer.hpp>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#define MAXSIZE 320 
#define MAXLEVEL 5 
#define NUM_CASTELJAU_POINTS 15
#define NUM_INITIAL_POINTS 10
#define FAR_PLANE 100.0
#define NEAR_PLANE -10.0
typedef struct Vertex {
	float XYZW[4];
	float RGBA[4];
	void SetCoords(float* coords) {
		XYZW[0] = coords[0];
		XYZW[1] = coords[1];
		XYZW[2] = coords[2];
		XYZW[3] = coords[3];
	}
	void SetColor(float* color) {
		RGBA[0] = color[0];
		RGBA[1] = color[1];
		RGBA[2] = color[2];
		RGBA[3] = color[3];
	}
};

//typedef struct Vector2 {
//	float x;
//	float y;
//};
class Vector2
{
public:
	float x;
	float y;
	Vector2(float xValue, float yValue)
	{
		x = xValue;
		y = yValue;
	}
};

// ATTN: USE POINT STRUCTS FOR EASIER COMPUTATIONS
typedef struct point {
	float x, y, z;
	point(const float x = 0, const float y = 0, const float z = 0) : x(x), y(y), z(z){};
	point(float *coords) : x(coords[0]), y(coords[1]), z(coords[2]){};
	point operator -(const point& a)const {
		return point(x - a.x, y - a.y, z - a.z);
	}
	point operator +(const point& a)const {
		return point(x + a.x, y + a.y, z + a.z);
	}
	point operator *(const float& a)const {
		return point(x*a, y*a, z*a);
	}
	point operator /(const float& a)const {
		return point(x / a, y / a, z / a);
	}
	float* toArray() {
		float array[] = { x, y, z, 1.0f };
		return array;
	}
};

enum LineType { Subdivide, Bezier, Catmul };

// function prototypes
int initWindow(void);
void initOpenGL(void);
void createVAOs(Vertex[], unsigned short[], size_t, size_t, int);
void createObjects(void);
void pickVertex(void);
void moveVertex(void);
void drawScene(void);
void cleanup(void);
static void mouseCallback(GLFWwindow*, int, int, int);
static void keyCallback(GLFWwindow*, int, int, int, int);

// GLOBAL VARIABLES
GLFWwindow* window;
const GLuint window_width = 1024, window_height = 768;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;
glm::mat4 gViewMatrix2;

GLuint gPickedIndex;
std::string gMessage;

GLuint programID;
GLuint pickingProgramID;

// ATTN: INCREASE THIS NUMBER AS YOU CREATE NEW OBJECTS
const GLuint NumObjects = 4;	// number of different "objects" to be drawn
GLuint VertexArrayId[NumObjects] = { 0, 1, 2, 3 };
GLuint VertexBufferId[NumObjects] = { 0, 1, 2, 3 };
GLuint IndexBufferId[NumObjects] = { 0, 1, 2, 3 };
size_t NumVert[NumObjects] = { 0, 0, 0, 0 };

GLuint MatrixID;
GLuint ViewMatrixID;
GLuint ModelMatrixID;
GLuint PickingMatrixID;
GLuint pickingColorArrayID;
GLuint pickingColorID;
GLuint LightID;

// Store time
// double lastTime = glfwGetTime();
bool moveOnZAxis = false;
bool doubleView = false;
bool movingPointActive = false;

float color[4];

float red[4] = { 1,0,0,1 };
float green[4] = { 0,1,0,1 };
float yellow[4] = { 1,1,0,1 };

// Subdivision and Restraints
LineType lineType = LineType::Subdivide;
int numSubdivisions = 0;

// our vip global arrays for picking/etc
float pickingColor[NUM_INITIAL_POINTS];
unsigned short Indices[NUM_INITIAL_POINTS];
Vertex Vertices[NUM_INITIAL_POINTS];

// Control Points
const int MaxControlPoints = 4 * NUM_INITIAL_POINTS;
unsigned short controlInds[MaxControlPoints];
Vertex controlVerts[MaxControlPoints];
GLsizei NumControlPoints = 0;

// Line Vertices
unsigned short lineInds[MAXSIZE];
Vertex lineVerts[MAXSIZE];
GLsizei NumLinePoints = 0;

// looper dot
Vertex looperVertex;
unsigned short looperIndex = 0;
float looperPos = 0;

/// <summary>
/// convert hsv values to an RGB vector
/// </summary>
std::vector<float> HSV_to_RGB(float h, float s, float v)
{
	float c = v * s; // c = v * s
	float x = c * (1 - abs(fmod((h / 60), 2) - 1));
	float m = v - c; // m = v - c

	std::vector<float> RGB;

	if (h < 60)
		RGB = { c, x, 0 };
	else if (h >= 60 && h < 120)
		RGB = { x, c, 0 };
	else if (h >= 120 && h < 180)
		RGB = { 0, c, x };
	else if (h >= 180 && h < 240)
		RGB = { 0, x, c };
	else if (h >= 240 && h < 300)
		RGB = { x, 0, c };
	else if (h >= 300 && h < 360)
		RGB = { c, 0, x };

	RGB[0] += m;
	RGB[1] += m;
	RGB[2] += m;

	return RGB;
}


/// <summary>
/// moves the point to the nearest available value
/// taken from: https://stackoverflow.com/a/21657330
/// </summary>
double clamp(double x, double upper, double lower)
{
	return min(upper, max(x, lower));
}

/// <summary>
/// moves the point to the nearest available value
/// taken from: https://stackoverflow.com/a/51018529/2994229
/// </summary>
int mod(int x, int m)
{
	int r = x % m;
	return r < 0 ? r + m : r;
}

/// <summary>
/// swaps rbg0 with rgb1 (both values remain stored)
/// </summary>
void swapColor(float rgb0[], float rgb1[])
{
	float temp[4] = { 0,0,0,0 };
	for (int i = 0; i < 4; i++)
	{
		temp[i] = rgb0[i];
		rgb0[i] = rgb1[i];
		rgb1[i] = temp[i];
	}
}


/// <summary>
/// Set colors evenly spaced by hues given RGBA[], 
/// index and total number of vertices
/// Algorithm: https://www.rapidtables.com/convert/color/hsv-to-rgb.html
/// </summary>
void setColor(float RGBA[], int index, int numPts)
{
	// Allow hues from red to magenta
	// Saturation and Value are set to 100%
	float h = ((float)((float)index / (numPts - 1)) * 300.0f);
	float s = 1.0f;
	float v = 1.0f;

	std::vector<float> RGB = HSV_to_RGB(h, s, v);

	// set rgba values
	RGBA[0] = RGB[0];
	RGBA[1] = RGB[1];
	RGBA[2] = RGB[2];
	RGBA[3] = 1;

	/*std::cout << "H: " << h << std::endl;*/
}

/// <summary>
/// Arrange (numPts) vertices in a circle of a specified radius
/// </summary>
void setCircleVertices(Vertex arr[], float radius, int startPoint, int numPts, float yOffset, int stopPoint)
{
	// 2 * PI radians in a circle
	// get the delta for angling the pts
	const float PI = 3.14159f;
	float delta_angle = (yOffset * 2 * PI) / (numPts);
	float angle = -1 * yOffset *PI/2;

	for (int i = startPoint; i < stopPoint; i++)
	{
		Vertex v = { {radius * cos(angle), radius * sin(angle) + yOffset, 0.0f, 1.0f} , { 255, 182.143, 0, 1 } };
		setColor(v.RGBA, i, 10);
		arr[i] = v;
		angle += delta_angle;
	}
}

void setFigure8Vertices(Vertex arr[], float circleRadius, int numPts, float topOffset, float botOffset)
{
	int halfPoint = numPts / 2;
	setCircleVertices(arr, circleRadius, 0, numPts/2, topOffset, halfPoint);
	setCircleVertices(arr, circleRadius, halfPoint, numPts/2, botOffset, numPts);

	for (int i = 0; i < NUM_INITIAL_POINTS; i++) {
		Indices[i] = i;
		pickingColor[i] = i / 255.0f;
	}
	for (int i = 0; i < MAXSIZE; i++) {
		lineInds[i] = i;
	}
	for (int i = 0; i < MaxControlPoints; i++) {
		controlInds[i] = i;
	}
}

void storePointsToArr(point points[], int numPoints, Vertex VerticesArr[], float color[4])
{
	for (int i = 0; i < numPoints; i++)
	{
		VerticesArr[i].SetCoords(points[i].toArray());
		VerticesArr[i].SetColor(color);
	}
}

void initCasteljau(point casteljau[][4], int n, point curve[])
{
	// input curve
	for (int x = 0; x < n; x++)
	{
		casteljau[x][n - x - 1] = curve[x];
	}
}

void evalCasteljau(point casteljau[][4], int n, point points[])
{
	for (int u = 0; u < NUM_CASTELJAU_POINTS; u++)
	{
		double t = (double)u / NUM_CASTELJAU_POINTS;
		for (int j = 1; j < n; j++)
		{
			for (int i = 0; i < n - j; i++)
			{
				casteljau[i][n - j - i - 1] = casteljau[i][n - j - i] * (1 - t) + casteljau[i + 1][n - j - i - 1] * t;
			}
		}
		points[u] = casteljau[0][0];
	}
}
/// <summary>
/// Adapted from: 
/// https://math.stackexchange.com/questions/43947/casteljaus-algorithm-practical-example
/// </summary>
void computeCasteljau(point points[], point curve[])
{
	// new 2d point arr of degree 4
	const int n = 4;
	point casteljau[n][n];

	initCasteljau(casteljau, n, curve);
	evalCasteljau(casteljau, n, points);
	
}

/// <summary> 
/// set the picking color
/// </summary>
void setPickingColor(float arr[], int indexCount)
{
	for (int i = 0; i < indexCount; i++)
	{
		arr[i] = (i / 255.0f);
	}
}

void doSubdivide(point (&current)[MAXSIZE], point(&controlPoints)[MaxControlPoints])
{
	int numSubdividePts = NUM_INITIAL_POINTS;
	point previous[MAXSIZE];
	
	// We use mod to stop the points from becoming unlinked
	// if we are at point 9, next is point 0

	for (int k = 0; k < numSubdivisions; k++)
	{
		// http://www.cplusplus.com/reference/algorithm/swap/
		std::swap(previous, current);
		for (int i = 0; i < numSubdividePts; i++)
		{
			//   |k |   //
			//  P|  |   //
			//   |2i|   //
			current[2 * i] = 
				(
					previous[i] * 4 + 
						previous[mod(i - 1, numSubdividePts)] * 4
				) / 8;
			//   |  k   |   //
			//  P|      |   //
			//   |2i + 1|   //
			current[2 * i + 1] = 
				(
					previous[mod(i + 1, numSubdividePts)]
						+ previous[mod(i - 1, numSubdividePts)] 
							+ previous[i] * 6
				) / 8;
		}
		numSubdividePts = numSubdividePts * 2;
	}
	NumLinePoints = numSubdividePts;
}

void doBezier(point(&points)[MAXSIZE])
{
	point bez[NUM_INITIAL_POINTS][4];

	// interior points
	for (int i = 0; i < NUM_INITIAL_POINTS; i++) 
	{
		// C[i][1]
		// We use mod to stop the points from becoming unlinked
		// if we are at point 9, next is point 0
		bez[i][1] = 
			(
				points[mod(i + 1, NUM_INITIAL_POINTS)] 
					+ points[i] * 2
			) / 3;
		bez[i][2] = 
			(
				points[mod(i + 1, NUM_INITIAL_POINTS)] * 2 
					+ points[i]
			) / 3;
	}

	// exterior points
	// we can calculate these by interpolating 
	// between the previously calculated values (AVERAGE)
	// We use mod to stop the points from becoming unlinked
	// if we are at point 0, next is point 9
	for (int i = 0; i < NUM_INITIAL_POINTS; i++) {
		bez[i][0] = 
			(
				bez[mod(i - 1, NUM_INITIAL_POINTS)][2] 
					+ bez[i][1]
			) / 2;
		bez[i][3] = 
			(
				bez[mod(i + 1, NUM_INITIAL_POINTS)][1] 
					+ bez[i][2]
			) / 2;
	}

	// evaluate using De Casteljau's Algorithm 
	for (int i = 0; i < NUM_INITIAL_POINTS; i++) {
		computeCasteljau(&points[i * NUM_CASTELJAU_POINTS], bez[i]);
	}

	// NumControlPoints = NUM_INITIAL_POINTS * 4;
	NumLinePoints = NUM_INITIAL_POINTS * NUM_CASTELJAU_POINTS;
	// std::cout << "nCtrlPts: " << NumControlPoints << ", nLinePts: "<< NumLinePoints"<< std::endl;
}

void doCatmul(point(&points)[MAXSIZE])
{
	// We use mod to stop the points from becoming unlinked
	// algorithm applied from http://algorithmist.net/docs/catmullrom.pdf

	point cat[NUM_INITIAL_POINTS][4];

	// alpha is set to 1/6 because alpha must be between 0 and 1 ... 
	// we have 2 points and we are using degree 3
	double alpha = (1.0 / 6.0);
	

	// exterior points (also control points)
	for (int i = 0; i < NUM_INITIAL_POINTS; i++) 
	{
		cat[i][0] = points[i];
		cat[i][3] = points[mod(i + 1, NUM_INITIAL_POINTS)];
	
	// get derivative (implicitly deffined tangent)
		point deltaL = 
			(
				points[mod(i + 1, NUM_INITIAL_POINTS)] 
					- points[mod(i - 1, NUM_INITIAL_POINTS)]
			) * alpha;
		point deltaR = 
			(
				points[i] 
					- points[mod(i + 2, NUM_INITIAL_POINTS)]
			) * alpha;
	
	// use derivative to interpolate interior points along tangents
		cat[i][1] = deltaL + cat[i][0];
		cat[i][2] = deltaR+ cat[i][3];
	}

	// evaluate curve
	// in it's own loop because it depends on all cat[] values
	for (int i = 0; i < NUM_INITIAL_POINTS; i++) {
		computeCasteljau(&points[i * NUM_CASTELJAU_POINTS], cat[i]);
	}

	// NumControlPoints = NUM_INITIAL_POINTS * 4;
	NumLinePoints = NUM_INITIAL_POINTS * NUM_CASTELJAU_POINTS;

	// std::cout << "nCtrlPts: " << NumControlPoints << ", nLinePts: "<< NumLinePoints"<< std::endl;
}

void evalObjects(point(&points)[MAXSIZE], point(&controlPoints)[MaxControlPoints])
{
	if (lineType == LineType::Subdivide){
		doSubdivide(points, controlPoints);
		std::copy(std::begin(red), std::end(red), std::begin(color));
	} else if (lineType == LineType::Bezier){
		doBezier(points);
		std::copy(std::begin(yellow), std::end(yellow), std::begin(color));
	} else if (lineType == LineType::Catmul){
		doCatmul(points);
		std::copy(std::begin(green), std::end(green), std::begin(color));
	}
}

void doLoopingBehaviour(point points[])
{
	// change looping vertex position

	float vPos = looperPos * NumLinePoints;
	int looperVert = (int) vPos;

	float n = vPos - looperVert;
	float c = 1 - n;

	point currP = points[mod(looperVert, NumLinePoints)];
	point nextP = points[mod(1 + looperVert, NumLinePoints)];

	point midP = (nextP * n) + (currP * c);

	looperVertex.SetCoords(midP.toArray());

	looperPos += 0.0005;
	looperPos -= (int) looperPos;
}
void createObjects(void)
{
	point points[MAXSIZE];
	point controlPoints[MaxControlPoints];

	// get points from the initial vertices
	for (int i = 0; i < NUM_INITIAL_POINTS; i++)
		points[i] = point(Vertices[i].XYZW);

	evalObjects(points, controlPoints);

	storePointsToArr(points, NumLinePoints, lineVerts, color);

	doLoopingBehaviour(points);
}

void drawLine()
{
	glBindVertexArray(VertexArrayId[1]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[1]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(lineVerts), lineVerts);
	glDrawElements(GL_LINE_LOOP, NumLinePoints, GL_UNSIGNED_SHORT, (void*)0);
}

void drawControlPoints()
{
	glBindVertexArray(VertexArrayId[2]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[2]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(controlVerts), controlVerts);
}

void drawLoopingPoint()
{
	if (!movingPointActive)
	{
		return;
	}
	glBindVertexArray(VertexArrayId[3]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[3]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(looperVertex), &looperVertex);
	glDrawElements(GL_POINTS, 1, GL_UNSIGNED_SHORT, (void*)0);
}
void drawViewport(glm::mat4 ModelMatrix)
{
	glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

	// Send our transformation to the currently bound shader, 
	// in the "MVP" uniform
	glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVP[0][0]);
	glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
	glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
	glm::vec3 lightPos = glm::vec3(4, 4, 4);
	glUniform3f(LightID, lightPos.x, lightPos.y, lightPos.z);

	glEnable(GL_PROGRAM_POINT_SIZE);

	glBindVertexArray(VertexArrayId[0]);	// draw Vertices
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]);
	glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertices), Vertices);				// update buffer data
	//glDrawElements(GL_LINE_LOOP, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);
	glDrawElements(GL_POINTS, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);

	// ATTN: OTHER BINDING AND DRAWING COMMANDS GO HERE, one set per object:
	drawLine();
	drawControlPoints();
	drawLoopingPoint();

	glBindVertexArray(0);
}
void drawScene(void)
{
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(programID);
	{
		const float PI = 3.14159f;

		glm::mat4 ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		
		if (!doubleView)
		{
			drawViewport(ModelMatrix);
		}
		
		if (doubleView)
		{
			auto ModelMatrix1 = ModelMatrix;
			vec4 transformVector(2, 2, 2, 0);
			ModelMatrix1 = glm::translate(ModelMatrix1, glm::vec3(0, 1, 0));
			ModelMatrix1 = glm::scale(ModelMatrix1, glm::vec3(0.5, 0.5, 0.5));
			drawViewport(ModelMatrix1);
			// ModelMatrix1 *= transformVector;

			auto ModelMatrix2 = ModelMatrix;
			vec3 EulerAngles(0, -PI / 2.0, 0);
			quat rot = quat(EulerAngles);
			mat4 RotationMatrix = glm::tmat4x4<float, glm::highp>(rot);
			ModelMatrix2 = glm::translate(ModelMatrix2, glm::vec3(0, -1, 0));
			ModelMatrix2 *= RotationMatrix;
			ModelMatrix2 = glm::scale(ModelMatrix2, glm::vec3(0.5, 0.5, 0.5));
			drawViewport(ModelMatrix2);
		}
	}
	glUseProgram(0);
	// Draw GUI
	TwDraw();

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
}

void pickVertex(void)
{
	// Clear the screen in white
	glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(pickingProgramID);
	{
		glm::mat4 ModelMatrix = glm::mat4(1.0); // TranslationMatrix * RotationMatrix;
		glm::mat4 MVP = gProjectionMatrix * gViewMatrix * ModelMatrix;

		// Send our transformation to the currently bound shader, in the "MVP" uniform
		glUniformMatrix4fv(PickingMatrixID, 1, GL_FALSE, &MVP[0][0]);
		glUniform1fv(pickingColorArrayID, NumVert[0], pickingColor);	// here we pass in the picking marker array

		// Draw the ponts
		glEnable(GL_PROGRAM_POINT_SIZE);
		glBindVertexArray(VertexArrayId[0]);
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[0]);
		glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(Vertices), Vertices);	// update buffer data
		glDrawElements(GL_POINTS, NumVert[0], GL_UNSIGNED_SHORT, (void*)0);
		glBindVertexArray(0);
	}
	glUseProgram(0);
	// Wait until all the pending drawing commands are really done.
	// Ultra-mega-over slow ! 
	// There are usually a long time between glDrawElements() and
	// all the fragments completely rasterized.
	glFlush();
	glFinish();

	glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

	// Read the pixel at the center of the screen.
	// You can also use glfwGetMousePos().
	// Ultra-mega-over slow too, even for 1 pixel, 
	// because the framebuffer is on the GPU.
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
	unsigned char data[4];
	glReadPixels(xpos, window_height - ypos, 1, 1, GL_RGBA, GL_UNSIGNED_BYTE, data); // OpenGL renders with (0,0) on bottom, mouse reports with (0,0) on top

	int pickedIndex = int(data[0]);
	// getaround values over 128 for border clicking
	// (WARNING) limits to 128 objects
	// (JUSTIFICATION) We only need 10 today :)
	if (pickedIndex < 128 || pickedIndex == 255) 
	{
		gPickedIndex = pickedIndex;
	}
}

// fill this function in!
void moveVertex(void)
{
	glm::mat4 ModelMatrix = glm::mat4(1.0);
	GLint viewport[4];
	glGetIntegerv(GL_VIEWPORT, viewport);
	glm::vec4 vp = glm::vec4(viewport[0], viewport[1], viewport[2], viewport[3]);
	
	// retrieve your cursor position
	double xpos, ypos;
	
	glfwGetCursorPos(window, &xpos, &ypos);

	double zpos = ypos;

	xpos = clamp(xpos, window_width - 10, 10);
	ypos = clamp(ypos, window_height - 10, 10);


	// get your world coordinates
	glm::vec3 newPos;

	// move points

	if (gPickedIndex == 255){ // Full white, must be the background !
		gMessage = "background";
	}
	else {
		newPos = glm::unProject(glm::vec3(xpos, viewport[3] - ypos, 0.0), ModelMatrix * gViewMatrix, gProjectionMatrix, vp);
		if (moveOnZAxis)
		{
			int scale = 1;
			ypos -= 10;
			zpos = ypos / (window_height - 20);
			if (zpos <= .5)
			{
				zpos *= 2;
				zpos *= 7.8;
				zpos -= 7.8;
				zpos = abs(zpos);
				zpos = clamp(zpos, 7.8, 0);
			}
			else
			{
				zpos *= 2;
				zpos -= 1;
				zpos *= -7.8;
			}
			newPos.z = zpos; //- newPos.y;
			// newPos.z = clamp(newPos.z, 7.8, -7.8);

			newPos.x = Vertices[gPickedIndex].XYZW[0];
			newPos.y = Vertices[gPickedIndex].XYZW[1];
		}
		else
		{
			newPos.z = Vertices[gPickedIndex].XYZW[2];
		}
		float newPosArr[4] = { newPos.x, newPos.y, newPos.z, 1.0 };
		Vertices[gPickedIndex].SetCoords(newPosArr);
		// Vertices[gPickedIndex].SetCoords(newPosArr);
		for (auto& point : Vertices[gPickedIndex].XYZW)
		{
			std::cout << point << ",";
		}
		std::cout << std::endl;
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();
	}
}

int initWindow(void)
{
	// Initialise GLFW
	if (!glfwInit()) {
		fprintf(stderr, "Failed to initialize GLFW\n");
		return -1;
	}

	glfwWindowHint(GLFW_SAMPLES, 4);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
	glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // FOR MAC

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "O'Connell, Michael (9990-6691)", NULL, NULL);
	if (window == NULL) {
		fprintf(stderr, "Failed to open GLFW window. If you have an Intel GPU, they are not 3.3 compatible. Try the 2.1 version of the tutorials.\n");
		glfwTerminate();
		return -1;
	}
	glfwMakeContextCurrent(window);

	// Initialize GLEW
	glewExperimental = true; // Needed for core profile
	if (glewInit() != GLEW_OK) {
		fprintf(stderr, "Failed to initialize GLEW\n");
		return -1;
	}

	// Initialize the GUI
	TwInit(TW_OPENGL_CORE, NULL);
	TwWindowSize(window_width, window_height);
	TwBar * GUI = TwNewBar("Picking");
	TwSetParam(GUI, NULL, "refresh", TW_PARAM_CSTRING, 1, "0.1");
	TwAddVarRW(GUI, "Last picked object", TW_TYPE_STDSTRING, &gMessage, NULL);

	// Set up inputs
	glfwSetInputMode(window, GLFW_STICKY_KEYS, GL_FALSE);
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetMouseButtonCallback(window, mouseCallback);
	glfwSetKeyCallback(window, keyCallback);

	return 0;
}

void initOpenGL(void)
{
	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.4f, 0.0f);

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	// Projection matrix : 45� Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	// glm::mat4 ProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f
	// gProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, (float) NEAR_PLANE, (float) FAR_PLANE);
	// Or, for an ortho camera :
	gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, (float) NEAR_PLANE, (float) FAR_PLANE); // In world coordinates

	// Camera matrix
	gViewMatrix = glm::lookAt(
		glm::vec3(0, 0, -5), // Camera is at (4,3,3), in World Space
		glm::vec3(0, 0, 0), // and looks at the origin
		glm::vec3(0, 1, 0)  // Head is up (set to 0,-1,0 to look upside-down)
		);

	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("hw1bShade.vertexshader", "hw1bShade.fragmentshader");  //programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("hw1bPick.vertexshader", "hw1bPick.fragmentshader"); //pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	// Get a handle for our "pickingColorID" uniform
	pickingColorArrayID = glGetUniformLocation(pickingProgramID, "PickingColorArray");
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");
	// Get a handle for our "LightPosition" uniform
	LightID = glGetUniformLocation(programID, "LightPosition_worldspace");

	createVAOs(Vertices, Indices, sizeof(Vertices), sizeof(Indices), 0);
	createObjects();

	// ATTN: create VAOs for each of the newly created objects here:
	// createVAOs(<fill this appropriately>);

	createVAOs(lineVerts, lineInds, sizeof(lineVerts), sizeof(lineInds), 1);
	createVAOs(controlVerts, controlInds, sizeof(controlVerts), sizeof(controlInds), 2);

	createVAOs(&looperVertex, &looperIndex, sizeof(looperVertex), sizeof(looperIndex), 3);
}

void createVAOs(Vertex Vertices[], unsigned short Indices[], size_t BufferSize, size_t IdxBufferSize, int ObjectId) {

	NumVert[ObjectId] = IdxBufferSize / (sizeof GLubyte);

	GLenum ErrorCheckValue = glGetError();
	size_t VertexSize = sizeof(Vertices[0]);
	size_t RgbOffset = sizeof(Vertices[0].XYZW);

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);
	glBindVertexArray(VertexArrayId[ObjectId]);

	// Create Buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, BufferSize, Vertices, GL_STATIC_DRAW);

	// Create Buffer for indices
	glGenBuffers(1, &IndexBufferId[ObjectId]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
	glBufferData(GL_ELEMENT_ARRAY_BUFFER, IdxBufferSize, Indices, GL_STATIC_DRAW);

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color

	// Disable our Vertex Buffer Object 
	glBindVertexArray(0);

	ErrorCheckValue = glGetError();
	if (ErrorCheckValue != GL_NO_ERROR)
	{
		fprintf(
			stderr,
			"ERROR: Could not create a VBO: %s \n",
			gluErrorString(ErrorCheckValue)
			);
	}
}

void cleanup(void)
{
	// Cleanup VBO and shader
	for (int i = 0; i < NumObjects; i++) {
		glDeleteBuffers(1, &VertexBufferId[i]);
		glDeleteBuffers(1, &IndexBufferId[i]);
		glDeleteVertexArrays(1, &VertexArrayId[i]);
	}
	glDeleteProgram(programID);
	glDeleteProgram(pickingProgramID);

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		pickVertex();
	}
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	if (action == GLFW_PRESS)
	{
		// std::cout << "pressed: "<< key << std::endl;
		if (key == GLFW_KEY_1)
		{
			if (lineType == LineType::Subdivide)
			{
				numSubdivisions = (numSubdivisions + 1) % (MAXLEVEL + 1);
			}
			lineType = LineType::Subdivide;
		}
		else if (key == GLFW_KEY_2)
		{
			lineType = LineType::Bezier;
		}
		else if (key == GLFW_KEY_3)
		{
			lineType = LineType::Catmul;
		}
		else if (key == GLFW_KEY_RIGHT_SHIFT || key == GLFW_KEY_LEFT_SHIFT)
		{
			moveOnZAxis = !moveOnZAxis;
		}
		else if (key == GLFW_KEY_4)
		{
			doubleView = !doubleView;
		}
		else if (key == GLFW_KEY_5)
		{
			movingPointActive = !movingPointActive;
		}

	}
}

int main(void)
{
	// initialize vertices into figure 8
	
	// pi/2 , -pi/2
	// +y, -y

	// Vector2 vec0 = Vector2(0, 1);
	// Vector2 vec1 = Vector2(0, -1);
	// initialize window
	setFigure8Vertices(Vertices, 1, 10, 1, -1);
	looperVertex.SetColor(yellow);
	
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	// initialize OpenGL pipeline
	initOpenGL();

	// For speed computation
	double lastTime = glfwGetTime();
	int nbFrames = 0;
	do {
		// Measure speed
		double currentTime = glfwGetTime();
		nbFrames++;
		if (currentTime - lastTime >= 1.0){ // If last prinf() was more than 1sec ago
			// printf and reset
			printf("%f ms/frame\n", 1000.0 / double(nbFrames));
			nbFrames = 0;
			lastTime += 1.0;
		}

		// DRAGGING: move current (picked) vertex with cursor
		if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT))
			moveVertex();

		// DRAWING SCENE
		createObjects();	// re-evaluate curves in case vertices have been moved
		drawScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}
