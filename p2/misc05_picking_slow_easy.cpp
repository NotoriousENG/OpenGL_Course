// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <stack>   
#include <sstream>
#include <string>
#define _USE_MATH_DEFINES
#include <math.h>
// Include GLEW
#include <GL/glew.h>
// Include GLFW
#include <GLFW/glfw3.h>
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
#include <iostream>
#define TINYOBJLOADER_IMPLEMENTATION
#include "tiny_obj_loader.h"
#include <glm/gtx/matrix_decompose.hpp>


const int window_width = 1024, window_height = 768;

enum class SceneObjects
{
	None,
	Camera,
	Pen,
	Base,
	Top,
	Arm1,
	Arm2
};
/// <summary>
/// position[4], color[4], normal[3]
/// </summary>
typedef struct Vertex {
	float Position[4];
	float Color[4];
	float Normal[3];
	void SetPosition(float* coords) {
		Position[0] = coords[0];
		Position[1] = coords[1];
		Position[2] = coords[2];
		Position[3] = 1.0;
	}
	void SetColor(float* color) {
		Color[0] = color[0];
		Color[1] = color[1];
		Color[2] = color[2];
		Color[3] = color[3];
	}
	void SetNormal(float* coords) {
		Normal[0] = coords[0];
		Normal[1] = coords[1];
		Normal[2] = coords[2];
	}
};

// function prototypes
int initWindow(void);
void initOpenGL(void);
void loadObject(char*, glm::vec4, Vertex*&, GLushort*&, int);
void createVAOs(Vertex Vertices[], unsigned short Indices[], int ObjectId);
void createObjects(void);
void pickObject(void);
void renderScene(void);
void cleanup(void);
static void keyCallback(GLFWwindow*, int, int, int, int);
static void mouseCallback(GLFWwindow*, int, int, int);

// GLOBAL VARIABLES
GLFWwindow* window;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;

GLuint gPickedIndex = -1;
std::string gMessage;

GLuint programID;
GLuint pickingProgramID;

const GLuint NumObjects = 10;	// ATTN: THIS NEEDS TO CHANGE AS YOU ADD NEW OBJECTS
GLuint VertexArrayId[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
GLuint VertexBufferId[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
GLuint IndexBufferId[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

size_t NumIndices[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
size_t VertexBufferSize[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };
size_t IndexBufferSize[NumObjects] = { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 };

GLuint MatrixID;
GLuint ModelMatrixID;
GLuint ViewMatrixID;
GLuint ProjMatrixID;
GLuint PickingMatrixID;
GLuint pickingColorID;
GLuint LightID1;
GLuint LightID2;
GLuint CamPosID;

GLint gX = 0.0;
GLint gZ = 0.0;

// animation control
bool animation = false;
// GLfloat phi = 0.0;

float colorWhite[] = { 1, 1, 1, 1 };
float defaultNormal[] = { 0, 0 , 1 };

glm::vec3 camPos = vec3(-10, 10, 10);
glm::vec3 camSpherical = vec3(-10, 10, 10);
glm::vec3 origin = glm::vec3(0, 0, 0);

glm::vec3 worldUp = glm::vec3(0, 1, 0);
glm::vec3 worldFront = glm::vec3(0, 0, -1);

SceneObjects selectedObj = SceneObjects::None;

float cameraSpeed = 5;

float deltaTime = 0.0f;
float lastFrame = 0.0f;

const float PI = 3.1417f;

float theta = 0;
float phi = 0;

// for loading in object
std::string inputfile = "HapticDevice.obj";
char* secondfile = "ball.obj";
tinyobj::attrib_t attrib;

std::vector<tinyobj::shape_t> shapes;
std::vector<tinyobj::material_t> materials;

std::string warn;
std::string err;

typedef std::vector<Vertex> objVertices;
typedef std::vector<uint16_t> objIndices;

std::vector<objVertices> objVerticesList;
std::vector<objIndices> objIndicesList;

glm::mat4x4 BaseModelMatrix = glm::mat4(1.0);
glm::mat4x4 TopModelMatrix = glm::mat4(1.0);
glm::mat4x4 Arm1ModelMatrix = glm::mat4(1.0);
glm::mat4x4 JointModelMatrix = glm::mat4(1.0);
glm::mat4x4 Arm2ModelMatrix = glm::mat4(1.0);
glm::mat4x4 PenModelMatrix = glm::mat4(1.0);
glm::mat4x4 ButtonModelMatrix = glm::mat4(1.0);

//std::vector< glm::mat4x4> Hierarchy
//{
//	BaseModelMatrix, TopModelMatrix, Arm1ModelMatrix, JointModelMatrix,
//	Arm2ModelMatrix, PenModelMatrix, ButtonModelMatrix
//};
//

std::vector<vec3> Positions =
{
	vec3(0,0,0),vec3(0,.1,0),vec3(0,.15,0),vec3(0,.25,0),vec3(0,.05,0),vec3(0,.14,0),vec3(0.015,.1,0), vec3(0,0,0)
	/*vec3(0,.094574, 0),
	vec3(0,1.0911, 0),
	vec3(0,1.5303, 0),
	vec3(0,3.3124, 0),
	vec3(0,3.6398, 0),
	vec3(0,45249, 0),
	vec3(0.002362, 4.5249,  0.42222)*/
};

std::vector<vec3> EulerRotations =
{
	vec3(0,0,0),vec3(0,0,0),vec3(36,0,0),vec3(0,0,0),vec3(2,0,0),vec3(0,0,0),vec3(0,0,0),vec3(0,0,0)
	/*vec3(0,0,0),
	vec3(0,0,0),
	vec3(36,-0,0),
	vec3(129,-0,0),
	vec3(0,0,0),
	vec3(-110,-0,-0),
	vec3(0,0,0)*/
};

std::vector<vec3> h_pos;


glm::vec2 inputVector;

float moveSpeed = 50;
float rotationSpeed = 60;

bool canTwistPen = false;

glm::vec3 temp_scale;
glm::quat temp_rotation;
glm::vec3 temp_translation;
glm::vec3 temp_skew;
glm::vec4 temp_perspective;

float pickingColor[NumObjects];

Vertex* monkeyVerts;
GLushort* monkeyInds;

bool isJumping;
vec3 JumpPos = vec3(0, 0, 0);
vec3 JumpVelocity = vec3(0, 0, 0);
vec3 Gravity = vec3(0, -4, 0);
float initialJumpVelocity = 4;

glm::vec3 sphericalToCartesian(vec3 sphericalCoords) {
	auto rho = sphericalCoords.x;
	auto theta = sphericalCoords.y;
	auto phi = sphericalCoords.z;
	return glm::vec3(rho * sin(phi) * sin(theta), rho * cos(phi), rho * sin(phi) * cos(theta));
}

glm::vec3 cartesianToSpherical(vec3 cartesianCoords) {
	auto x = cartesianCoords.x;
	auto y = cartesianCoords.y;
	auto z = cartesianCoords.z;
	return glm::vec3(sqrt(x * x + z * z + y * y), atan(sqrt(x * x + z * z) / y), atan(z / x));
}

void swapYZ(std::vector<vec3> v)
{
	for (int i = 0; i < v.size(); i++)
	{
		auto temp = v[i].y;
		v[i].y = v[i].z;
		v[i].z = temp;
	}
}

void swapPosRotYZ()
{
	swapYZ(Positions);
	swapYZ(EulerRotations);
}

void updatePositions(int rootObj, vec3 speed)
{
	for (int s = 0; s < shapes.size(); s++)
	{
		h_pos[s] += speed;
	}
}

void setCamera()
{
	camPos = sphericalToCartesian(camSpherical);
	gViewMatrix = glm::lookAt(camPos,	// eye
		origin,	// center
		worldUp);	// up
	// std::cout << "CAMERA: " << camPos.x << "," << camPos.y << "," << camPos.z << std::endl;
}

void moveCamera()
{
	if (inputVector.length() == 0)
	{
		camSpherical.z = clamp(camSpherical.z, 0.001f,(float) M_PI - 0.001f);
		return;
	}

	double v = (double)cameraSpeed * deltaTime;
	if (inputVector.y < 0)
	{
		// Longitude
		camSpherical.z = (float)min(camSpherical.z + v, M_PI - 0.001);
	}
	if (inputVector.y > 0)
	{
		camSpherical.z = (float)max(camSpherical.z - v, 0.001);
	}
	if (inputVector.x > 0)
	{
		// Latitude
		camSpherical.y = (float)fmod(camSpherical.y + v, M_PI * 2);
	}
	if (inputVector.x < 0)
	{
		camSpherical.y = (float)fmod(camSpherical.y - v, M_PI * 2);
	}
}
void moveBase()
{
	vec3 translation = vec3(inputVector.x, 0, inputVector.y) * 0.01f/* moveSpeed * deltaTime*/;
	Positions[0] += translation;
}

void rotate_yAxis(int i, float speedMultiplier)
{
	if (inputVector.x != 0)
	{
		vec3 rotationAxis = vec3(0, inputVector.x, 0) * 0.01f;
		EulerRotations[i] += rotationAxis;
	}
}

void rotate_Up(int i, float speedMultiplier)
{
	if (inputVector.x != 0)
	{
		vec3 rotationAxis = vec3(inputVector.x, 0, 0) * 0.01f;
		EulerRotations[i] += rotationAxis;
	}
}
void rotate_Pen(int i)
{
	vec3 rotationAxis = vec3(inputVector.x, 0, inputVector.y) * 0.01f;
	EulerRotations[i] += rotationAxis;
}

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
/// Set colors evenly spaced by hues given RGBA[], 
/// index and total number of vertices
/// Algorithm: https://www.rapidtables.com/convert/color/hsv-to-rgb.html
/// </summary>
void setColor(float RGBA[], int index, int lastIndex)
{
	// Allow hues from red to magenta
	// Saturation and Value are set to 100%
	float h = ((float)((float)index / (lastIndex - 1)) * 300.0f);
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

void tinyLoadObj()
{
	bool ret = tinyobj::LoadObj(&attrib, &shapes, &materials, &warn, &err, inputfile.c_str());

	if (!warn.empty()) {
		std::cout << warn << std::endl;
	}

	if (!err.empty()) {
		std::cerr << err << std::endl;
	}

	if (!ret) {
		exit(1);
	}

	// Loop over shapes
	for (size_t s = 0; s < shapes.size(); s++) {
		objVertices verts;
		objIndices inds;
		// Loop over faces(polygon)
		size_t index_offset = 0;
		for (size_t f = 0; f < shapes[s].mesh.num_face_vertices.size(); f++) {
			int fv = shapes[s].mesh.num_face_vertices[f];


			// Loop over vertices in the face.
			for (size_t v = 0; v < fv; v++) {
				Vertex vertex = {};
				// access to vertex
				tinyobj::index_t idx = shapes[s].mesh.indices[index_offset + v];
				tinyobj::real_t vx = attrib.vertices[3 * idx.vertex_index + 0];
				tinyobj::real_t vy = attrib.vertices[3 * idx.vertex_index + 1];
				tinyobj::real_t vz = attrib.vertices[3 * idx.vertex_index + 2];
				tinyobj::real_t nx = attrib.normals[3 * idx.normal_index + 0];
				tinyobj::real_t ny = attrib.normals[3 * idx.normal_index + 1];
				tinyobj::real_t nz = attrib.normals[3 * idx.normal_index + 2];
				tinyobj::real_t tx = attrib.texcoords[2 * idx.texcoord_index + 0];
				tinyobj::real_t ty = attrib.texcoords[2 * idx.texcoord_index + 1];
				// Optional: vertex colors
				tinyobj::real_t red = attrib.colors[3 * idx.vertex_index + 0];
				tinyobj::real_t green = attrib.colors[3 * idx.vertex_index + 1];
				tinyobj::real_t blue = attrib.colors[3 * idx.vertex_index + 2];
				float coords[] = { vx, vy, vz };
				vertex.SetPosition(coords);

				// std::cout << vx << "," << vy << "," << vz << std::endl;

				float norms[] = { nx, ny, nz };
				vertex.SetNormal(norms);

				float cols[] = { red, green, blue, 1 };
				setColor(cols, s, shapes.size());
				vertex.SetColor(cols);

				verts.push_back(vertex);
				inds.push_back(inds.size());
			}
			index_offset += fv;

			// per-face material
			shapes[s].mesh.material_ids[f];
		}

		// for each shape

		NumIndices[2 + s] = inds.size();
		VertexBufferSize[2 + s] = sizeof(verts[0]) * verts.size();
		IndexBufferSize[2 + s] = sizeof(GLushort) * inds.size();

		objVerticesList.push_back(verts);
		objIndicesList.push_back(inds);
	}
}


void loadObject(char* file, glm::vec4 color, Vertex*& out_Vertices, GLushort*& out_Indices, int ObjectId)
{
	// Read our .obj file
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	bool res = loadOBJ(file, vertices, normals);

	std::vector<GLushort> indices;
	std::vector<glm::vec3> indexed_vertices;
	std::vector<glm::vec2> indexed_uvs;
	std::vector<glm::vec3> indexed_normals;
	indexVBO(vertices, normals, indices, indexed_vertices, indexed_normals);

	const size_t vertCount = indexed_vertices.size();
	const size_t idxCount = indices.size();

	// populate output arrays
	out_Vertices = new Vertex[vertCount];
	for (int i = 0; i < vertCount; i++) {
		out_Vertices[i].SetPosition(&indexed_vertices[i].x);
		out_Vertices[i].SetNormal(&indexed_normals[i].x);
		out_Vertices[i].SetColor(&color[0]);
	}
	out_Indices = new GLushort[idxCount];
	for (int i = 0; i < idxCount; i++) {
		out_Indices[i] = indices[i];
	}

	// set global variables!!
	NumIndices[ObjectId] = idxCount;
	VertexBufferSize[ObjectId] = sizeof(out_Vertices[0]) * vertCount;
	IndexBufferSize[ObjectId] = sizeof(GLushort) * idxCount;
}

void setGridVerts(Vertex& v, float pos[4])
{
	v.SetColor(colorWhite);
	v.SetNormal(defaultNormal);
	v.SetPosition(pos);
}

void createGrid(Vertex vertices[], int length)
{
	int vertexIndex = 0;

	// lines perpendidular to z axis
	for (int i = -length; i <= length; i++)
	{
		float start_X[] = { -length , 0 , i , 0 };
		setGridVerts(vertices[vertexIndex], start_X);
		vertexIndex++;

		float end_X[] = { length , 0 , i , 0 };
		setGridVerts(vertices[vertexIndex], end_X);
		vertexIndex++;

		float start_Z[] = { i , 0 , -length , 0 };
		setGridVerts(vertices[vertexIndex], start_Z);
		vertexIndex++;

		float end_Z[] = { i , 0 , length , 0 };
		setGridVerts(vertices[vertexIndex], end_Z);
		vertexIndex++;
	}
}

void createObjects(void)
{
	//-- COORDINATE AXES --//
	Vertex CoordVerts[] =
	{
		{ { 0.0, 0.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 5.0, 0.0, 0.0, 1.0 }, { 1.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 5.0, 0.0, 1.0 }, { 0.0, 1.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 0.0, 1.0 }, { 0.0, 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
		{ { 0.0, 0.0, 5.0, 1.0 }, { 0.0, 0.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } },
	};

	VertexBufferSize[0] = sizeof(CoordVerts);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(CoordVerts, NULL, 0);

	//-- GRID --//

	Vertex GridVerts[5 * 2 * 4 + 4];

	createGrid(GridVerts, 5);

	VertexBufferSize[1] = sizeof(GridVerts);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	createVAOs(GridVerts, NULL, 1);

	//-- .OBJs --//

	// load in our object
	tinyLoadObj();


	// createVAOs()

	// ATTN: load your models here
	//Vertex* Verts;
	//GLushort* Idcs;
	//loadObject("models/base.obj", glm::vec4(1.0, 0.0, 0.0, 1.0), Verts, Idcs, ObjectID);
	for (int s = 0; s < shapes.size(); s++)
	{
		int i = 0;
		createVAOs(&objVerticesList[s][0], &objIndicesList[s][0], 2 + s);
	}
	// createVAOs(&objVertices[0], &objIndices[0], 2);

	loadObject(secondfile, vec4(0,0,0,1), monkeyVerts, monkeyInds, 9);
	createVAOs(monkeyVerts, monkeyInds, 9);
}

void drawMesh()
{
	glm::mat4x4 ModelMatrix = glm::mat4(1.0);
	for (int s = 0; s < shapes.size(); s++)
	{
		for (int i = 0; i < shapes.size(); i++)
		{
			ModelMatrix = glm::translate(ModelMatrix, Positions[s]);
			ModelMatrix = glm::rotate(ModelMatrix, EulerRotations[s].x, { 1, 0, 0 });
			ModelMatrix = glm::rotate(ModelMatrix, EulerRotations[s].y, { 0, 1, 0 });
			ModelMatrix = glm::rotate(ModelMatrix, EulerRotations[s].z, { 0, 0, 1 });
		}

		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

		// draw Obj
		glBindVertexArray(VertexArrayId[2 + s]);
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[2 + s]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[2 + s]);

		glEnableVertexAttribArray(0); // position
		glEnableVertexAttribArray(1); // color
		glEnableVertexAttribArray(2); // normal

		// glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &Hierarchy[s][0][0]);

		glDrawElements(GL_TRIANGLES, IndexBufferSize[2 + s], GL_UNSIGNED_SHORT, (void*)0);
	}
	
}

void drawMonkey()
{
	if (isJumping)
	{
		glm::mat4x4 ModelMatrix = glm::mat4(1.0);
		ModelMatrix = glm::translate(ModelMatrix, JumpPos);

		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

		// draw Obj
		glBindVertexArray(VertexArrayId[9]);
		glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[9]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[9]);

		glEnableVertexAttribArray(0); // position
		glEnableVertexAttribArray(1); // color
		glEnableVertexAttribArray(2); // normal

		// glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &Hierarchy[s][0][0]);

		glDrawElements(GL_TRIANGLES, IndexBufferSize[9], GL_UNSIGNED_SHORT, (void*)0);
	}
}

void doJump() 
{
	isJumping = true;

	glm::mat4x4 backMat = glm::mat4(1.0);
	for (int i = 0; i < shapes.size() - 1; i++) {
		backMat = glm::translate(backMat, Positions[i]);
		backMat = glm::rotate(backMat, EulerRotations[i].x, { 1, 0, 0 });
		backMat = glm::rotate(backMat, EulerRotations[i].y, { 0, 1, 0 });
		backMat = glm::rotate(backMat, EulerRotations[i].z, { 0, 0, 1 });
	}
	glm::vec4 PenBack = backMat[3];

	glm::mat4x4 tipMat = glm::mat4(1.0);
	tipMat = glm::translate(backMat, { 0.0, -1.5, 0});
	glm::vec4 PenTip = tipMat[3];

	glm::vec4 tangent4 = PenBack - PenTip;
	glm::vec3 tangent3 = { tangent4.x, tangent4.y, tangent4.z };
	glm::vec3 tangent = glm::normalize(tangent3);

	JumpPos = { PenBack.x / PenBack.w, PenBack.y / PenBack.w, PenBack.z / PenBack.w };
	JumpVelocity = tangent * initialJumpVelocity;

}

void getJumpMovement()
{
	if (isJumping) {
		JumpPos += JumpVelocity * deltaTime;
		JumpVelocity += Gravity * deltaTime;
		if (JumpPos.y < 0) {
			isJumping = false;
			Positions[0] = vec3(JumpPos.x, 0, JumpPos.z);
		}
	}
}

void renderScene(void)
{
	//ATTN: DRAW YOUR SCENE HERE. MODIFY/ADAPT WHERE NECESSARY!


	// Dark blue background
	glClearColor(0.0f, 0.0f, 0.2f, 0.0f);
	// Re-clear the screen for real rendering
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	glUseProgram(programID);
	{
		glm::vec3 lightPos1 = glm::vec3(2, 4, 2);
		glUniform3f(LightID1, lightPos1.x, lightPos1.y, lightPos1.z);
		glm::vec3 lightPos2 = glm::vec3(1, 2, -2);
		glUniform3f(LightID2, lightPos2.x, lightPos2.y, lightPos2.z);

		glm::mat4x4 ModelMatrix = glm::mat4(1.0);
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
		glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
		vec3 camPos = sphericalToCartesian(camSpherical);

		glUniform3f(CamPosID, camPos.x, camPos.y, camPos.z);

		glBindVertexArray(VertexArrayId[0]);	// draw CoordAxes
		glDrawArrays(GL_LINES, 0, 6);

		glBindVertexArray(VertexArrayId[1]);	// draw CoordAxes

		glDrawArrays(GL_LINES, 0, 5 * 2 * 4 + 4);

		drawMesh();

		if (isJumping)
		{
			drawMonkey();
		}
		

		glBindVertexArray(0);

	}
	glUseProgram(0);
	// Draw GUI
	TwDraw();

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
}

void pickObject(void)
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

		// ATTN: DRAW YOUR PICKING SCENE HERE. REMEMBER TO SEND IN A DIFFERENT PICKING COLOR FOR EACH OBJECT BEFOREHAND
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

	// Convert the color back to an integer ID
	gPickedIndex = int(data[0]);

	if (gPickedIndex == 255) { // Full white, must be the background !
		gMessage = "background";
	}
	else {
		std::ostringstream oss;
		oss << "point " << gPickedIndex;
		gMessage = oss.str();
	}

	// Uncomment these lines to see the picking shader in effect
	//glfwSwapBuffers(window);
	//continue; // skips the normal rendering
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
	glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

	// Open a window and create its OpenGL context
	window = glfwCreateWindow(window_width, window_height, "O'Connell, Michael(9990-6691)", NULL, NULL);
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
	TwBar* GUI = TwNewBar("Picking");
	TwSetParam(GUI, NULL, "refresh", TW_PARAM_CSTRING, 1, "0.1");
	TwAddVarRW(GUI, "Last picked object", TW_TYPE_STDSTRING, &gMessage, NULL);

	// Set up inputs
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetKeyCallback(window, keyCallback);
	glfwSetMouseButtonCallback(window, mouseCallback);

	return 0;
}

void initOpenGL(void)
{

	// Enable depth test
	glEnable(GL_DEPTH_TEST);
	// Accept fragment if it closer to the camera than the former one
	glDepthFunc(GL_LESS);
	// Cull triangles which normal is not towards the camera
	glEnable(GL_CULL_FACE);

	// Projection matrix : 45� Field of View, 4:3 ratio, display range : 0.1 unit <-> 100 units
	gProjectionMatrix = glm::perspective(45.0f, 4.0f / 3.0f, 0.1f, 100.0f);
	// Or, for an ortho camera :
	//gProjectionMatrix = glm::ortho(-4.0f, 4.0f, -3.0f, 3.0f, 0.0f, 100.0f); // In world coordinates


	moveCamera();
	setCamera();

	//// Camera matrix
	//gViewMatrix = glm::lookAt(camPos,	// eye
	//	origin,	// center
	//	worldUp);	// up

	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");
	pickingProgramID = LoadShaders("Picking.vertexshader", "Picking.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ProjMatrixID = glGetUniformLocation(programID, "P");
	CamPosID = glGetUniformLocation(programID, "camPos");

	PickingMatrixID = glGetUniformLocation(pickingProgramID, "MVP");
	// Get a handle for our "pickingColorID" uniform
	pickingColorID = glGetUniformLocation(pickingProgramID, "PickingColor");
	// Get a handle for our "LightPosition" uniform
	LightID1 = glGetUniformLocation(programID, "LightPosition1_worldspace");
	LightID2 = glGetUniformLocation(programID, "LightPosition2_worldspace");

	createObjects();
}

void createVAOs(Vertex Vertices[], unsigned short Indices[], int ObjectId) {

	GLenum ErrorCheckValue = glGetError();
	const size_t VertexSize = sizeof(Vertices[0]);
	const size_t RgbOffset = sizeof(Vertices[0].Position);
	const size_t Normaloffset = sizeof(Vertices[0].Color) + RgbOffset;

	// Create Vertex Array Object
	glGenVertexArrays(1, &VertexArrayId[ObjectId]);	//
	glBindVertexArray(VertexArrayId[ObjectId]);		//

	// Create Buffer for vertex data
	glGenBuffers(1, &VertexBufferId[ObjectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[ObjectId]);
	glBufferData(GL_ARRAY_BUFFER, VertexBufferSize[ObjectId], Vertices, GL_STATIC_DRAW);

	// Create Buffer for indices
	if (Indices != NULL) {
		glGenBuffers(1, &IndexBufferId[ObjectId]);
		glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[ObjectId]);
		glBufferData(GL_ELEMENT_ARRAY_BUFFER, IndexBufferSize[ObjectId], Indices, GL_STATIC_DRAW);
	}

	// Assign vertex attributes
	glVertexAttribPointer(0, 4, GL_FLOAT, GL_FALSE, VertexSize, 0);
	glVertexAttribPointer(1, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)RgbOffset);
	glVertexAttribPointer(2, 4, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)Normaloffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color
	glEnableVertexAttribArray(2);	// normal

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

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// ATTN: MODIFY AS APPROPRIATE
	if (action == GLFW_PRESS) {
		switch (key)
		{
		case GLFW_KEY_P:
			selectedObj = SceneObjects::Pen;
			canTwistPen = false;
			break;
		case GLFW_KEY_B:
			selectedObj = SceneObjects::Base;
			//std::cout << "BASE" << std::endl;
			break;
		case GLFW_KEY_T:
			selectedObj = SceneObjects::Top;
			break;
		case GLFW_KEY_1:
			selectedObj = SceneObjects::Arm1;
			break;
		case GLFW_KEY_2:
			selectedObj = SceneObjects::Arm2;
			break;
		case GLFW_KEY_C:
			selectedObj = SceneObjects::Camera;
			break;
		case GLFW_KEY_J:
			if (!isJumping)
			{
				doJump();
			}
		case GLFW_KEY_LEFT_SHIFT:
			if (selectedObj == SceneObjects::Pen)
			{
				canTwistPen = !canTwistPen;
			}
			break;
		default:
			break;
		}


	}
	else if (action == GLFW_REPEAT)
	{
		switch (key)
		{
		case GLFW_KEY_DOWN:
			inputVector.y = -1;
			inputVector.x = 0;
			theta -= cameraSpeed * deltaTime;
			break;
		case GLFW_KEY_UP:
			inputVector.y = 1;
			inputVector.x = 0;
			theta += cameraSpeed * deltaTime;
			break;
		case GLFW_KEY_LEFT:
			inputVector.x = -1;
			inputVector.y = 0;
			phi -= cameraSpeed * deltaTime;
			break;
		case GLFW_KEY_RIGHT:
			inputVector.x = 1;
			inputVector.y = 0;
			phi += cameraSpeed * deltaTime;
			break;
		default:
			inputVector = glm::vec2(0, 0);
			break;
		}
		if (selectedObj == SceneObjects::Camera)
		{
			setCamera();
			moveCamera();
		}
		else if (selectedObj == SceneObjects::Base)
		{
			moveBase();
		}
		else if (selectedObj == SceneObjects::Top)
		{
			rotate_yAxis(1, .5f);
		}
		else if (selectedObj == SceneObjects::Arm1)
		{
			rotate_Up(2, 0.1f);
		}
		else if (selectedObj == SceneObjects::Arm2)
		{
			rotate_Up(4, 0.1f);
		}
		else if (selectedObj == SceneObjects::Pen)
		{
			if (!canTwistPen)
			{
				rotate_Pen(5);
			}
			else
			{
				//std::cout << "TWISTY TIME" << std::endl;
				rotate_yAxis(5, 0.01);
			}
		}

	}
	else if (action != GLFW_REPEAT)
	{
		inputVector = glm::vec2(0, 0);
	}
}

void set_deltaTime()
{
	float currentFrame = glfwGetTime();
	deltaTime = currentFrame - lastFrame;
	lastFrame = currentFrame;
}

static void mouseCallback(GLFWwindow* window, int button, int action, int mods)
{
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS) {
		pickObject();
	}
}

int main(void)
{
	swapPosRotYZ();

	// initialize window
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	cartesianToSpherical(vec3(10, 10, 10));

	// initialize OpenGL pipeline
	initOpenGL();

	/*for (auto& v : objVertices)
	{
		std::cout << v.Position[0] << "," << v.Position[1] << "," << v.Position[2] << "," << v.Position[3] << std::endl;
	}*/

	// For speed computation
	double lastTime = glfwGetTime();
	int nbFrames = 0;
	do {
		//// Measure speed
		//double currentTime = glfwGetTime();
		//nbFrames++;
		//if (currentTime - lastTime >= 1.0){ // If last prinf() was more than 1sec ago
		//	// printf and reset
		//	printf("%f ms/frame\n", 1000.0 / double(nbFrames));
		//	nbFrames = 0;
		//	lastTime += 1.0;
		//}

		setCamera();

		// std::cout << inputVector.x << "," << inputVector.y << std::endl;

		set_deltaTime();

		// getLookAt();

		if (animation) {
			phi += 0.01;
			if (phi > 360)
				phi -= 360;
		}

		getJumpMovement();

		// DRAWING POINTS
		renderScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}