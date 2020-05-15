// Include standard headers
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <array>
#include <stack>   
#include <sstream>
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

#include <common/shader.hpp>
#include <common/controls.hpp>
#include <common/objloader.hpp>
#include <common/vboindexer.hpp>
#include <common/texture.hpp>
#include <iostream>

// TODO: MAKE TEXTURE
// TODO: SUBDIV AAaaAAaa

const int window_width = 1024, window_height = 768;

typedef struct Vertex {
	float Position[4];
	float Color[4];
	float Normal[3];
	float UV[2];
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
	void SetUV(float* uv) {
		UV[0] = uv[0];
		UV[1] = uv[1];
	}
};


glm::vec3 sphericalToCartesian(float rho, float theta, float phi) {
	return glm::vec3(
		rho * sin(phi) * sin(theta),
		rho * cos(phi),
		rho * sin(phi) * cos(theta)
	);
}

glm::vec3 cartesianToSpherical(float x, float y, float z) {
	return glm::vec3(
		sqrt(x * x + z * z + y * y),
		atan(sqrt(x * x + z * z) / y),
		atan(z / x)
	);
}

// function prototypes
int initWindow(void);
void initOpenGL(void);
void loadObject(char*, glm::vec4, Vertex*&, GLushort*&, int);
void createVAOs(Vertex[], unsigned short[], int);
void createObjects(void);
void renderScene(void);
void cleanup(void);
static void keyCallback(GLFWwindow*, int, int, int, int);

// GLOBAL VARIABLES
GLFWwindow* window;

glm::mat4 gProjectionMatrix;
glm::mat4 gViewMatrix;

GLuint programID;

enum class Drawable {
	Axis = 0,
	Grid = 1,
	Head = 2,
	Head2 = 3
};

const GLuint NumObjects = 4;	// ATTN: THIS NEEDS TO CHANGE AS YOU ADD NEW OBJECTS

GLuint VertexArrayId[NumObjects] = { 0 };
GLuint VertexBufferId[NumObjects] = { 0 };
GLuint IndexBufferId[NumObjects] = { 0 };

size_t NumIndices[NumObjects] = { 0 };
size_t VertexBufferSize[NumObjects] = { 0 };
size_t IndexBufferSize[NumObjects] = { 0 };

Vertex* verts[NumObjects];
GLushort* idcs[NumObjects];

GLuint MatrixID;
GLuint ModelMatrixID;
GLuint ViewMatrixID;
GLuint ProjMatrixID;
GLuint CameraPosID;
GLuint TextureID;
GLboolean UseTextureID;

// Grid
const int GRID_WIDTH = 5;
const int NUM_GRID_VERTICES = (GRID_WIDTH * 8) + 4;

// Camera Pos
float CameraRadius;
float CameraLatitude;
float CameraLongitude;
double CameraMoveSpeed = 3;

// Relative Object Positions
glm::vec3 Positions[NumObjects] = {
	{ 0.0, 0.0, 0.0 }, // Axis
	{ 0.0, 0.0, 0.0 }, // Grid
	{ 0.0, 0.0, 0.0 } // Head
};

// Object Rotations (radians)
glm::vec3 Rotations[NumObjects] = {
	{ 0.0, 0.0, 0.0 }, // Axis
	{ 0.0, 0.0, 0.0 }, // Grid
	{ 0.0, 0.0, 0.0 } // Head
};

// Object Scale
glm::vec3 Scalings[NumObjects] = {
	{ 1.0, 1.0, 1.0 }, // Axis
	{ 1.0, 1.0, 1.0 }, // Grid
	{ 1, 1, 1 } // Head
};

// Head Model
bool ShowHead;
GLuint Texture;

//float mix1d(float a, float b, float t)
//{
//	return a * (1.0f - t) + b * t;
//}
//
//float* mix3d(float* a, float* b, float t)
//{
//	// degree 1
//	auto ax = mix1d(a[0], b[0], t);
//	auto ay = mix1d(a[1], b[1], t);
//	auto az = mix1d(a[2], b[3], t);
//	float mixed[3] = { ax, ay, az };
//}
//
//float* Bezier(float* A, float* B, float* C, float t)
//{
//	//degree 2
//	float* AB = mix3d(A, B, t);
//	float* BC = mix3d(B, C, t);
//	return mix3d(AB, BC, t);
//}

void setInitialVals() {
	ShowHead = true;
	glm::vec3 initialCameraPos = cartesianToSpherical(10, 10, 10);
	CameraRadius = initialCameraPos.x;
	CameraLatitude = initialCameraPos.y;
	CameraLongitude = initialCameraPos.z;

}



void loadObject(char* file, glm::vec4 color, Vertex*& out_Vertices, GLushort*& out_Indices, int ObjectId)
{
	// Read our .obj file
	std::vector<glm::vec3> vertices;
	std::vector<glm::vec2> uvs;
	std::vector<glm::vec3> normals;
	bool res = loadOBJ(file, vertices, uvs, normals);

	std::vector<GLushort> indices;
	std::vector<glm::vec3> indexed_vertices;
	std::vector<glm::vec2> indexed_uvs;
	std::vector<glm::vec3> indexed_normals;
	indexVBO(vertices, normals, uvs, indices, indexed_vertices, indexed_normals, indexed_uvs);

	const size_t vertCount = indexed_vertices.size();
	const size_t idxCount = indices.size();

	// populate output arrays
	out_Vertices = new Vertex[vertCount];
	for (int i = 0; i < vertCount; i++) {
		out_Vertices[i].SetPosition(&indexed_vertices[i].x);
		out_Vertices[i].SetColor(&color[0]);
		out_Vertices[i].SetNormal(&indexed_normals[i].x);
		out_Vertices[i].SetUV(&indexed_uvs[i].x);
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

	unsigned short coordIndices[6];
	for (int i = 0; i < 6; i++) {
		coordIndices[i] = i;
	}

	VertexBufferSize[(int)Drawable::Axis] = sizeof(CoordVerts);	// ATTN: this needs to be done for each hand-made object with the ObjectID (subscript)
	IndexBufferSize[(int)Drawable::Axis] = sizeof(coordIndices);
	NumIndices[(int)Drawable::Axis] = 6;
	verts[(int)Drawable::Axis] = CoordVerts;
	idcs[(int)Drawable::Axis] = coordIndices;

	createVAOs(CoordVerts, coordIndices, 0);

	//-- GRID --//
	Vertex GridVerts[NUM_GRID_VERTICES];
	int vertIndex = 0;
	for (int i = -GRID_WIDTH; i <= GRID_WIDTH; i++) {
		Vertex xVertexBack = { { -GRID_WIDTH, 0.0, i, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
		Vertex xVertexFront = { { GRID_WIDTH, 0.0, i, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
		Vertex zVertexBack = { { i, 0.0, -GRID_WIDTH, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
		Vertex zVertexFront = { { i, 0.0, GRID_WIDTH, 1.0 }, { 1.0, 1.0, 1.0, 1.0 }, { 0.0, 0.0, 1.0 } };
		GridVerts[vertIndex++] = xVertexBack;
		GridVerts[vertIndex++] = xVertexFront;
		GridVerts[vertIndex++] = zVertexBack;
		GridVerts[vertIndex++] = zVertexFront;
	}

	unsigned short gridIndices[NUM_GRID_VERTICES];
	for (int i = 0; i < NUM_GRID_VERTICES; i++) {
		gridIndices[i] = i;
	}

	VertexBufferSize[(int)Drawable::Grid] = sizeof(GridVerts);
	IndexBufferSize[(int)Drawable::Grid] = sizeof(gridIndices);
	NumIndices[(int)Drawable::Grid] = NUM_GRID_VERTICES;
	verts[(int)Drawable::Grid] = GridVerts;
	idcs[(int)Drawable::Grid] = gridIndices;

	createVAOs(GridVerts, gridIndices, 1);

	//-- .OBJs --//
	// ATTN: load your models here
	loadObject("finalface.obj", { 0,0,0,0 }, verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
	createVAOs(verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
}

void drawObject(int objectId, GLenum mode) {
	glBindVertexArray(VertexArrayId[objectId]);
	glBindBuffer(GL_ARRAY_BUFFER, VertexBufferId[objectId]);
	glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, IndexBufferId[objectId]);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color
	glEnableVertexAttribArray(2);	// normal
	glEnableVertexAttribArray(3);	// uv

	glDrawElements(mode, IndexBufferSize[objectId], GL_UNSIGNED_SHORT, (void*)0);
}

void drawLines(Drawable drawable) {
	glUniform1i(UseTextureID, false);
	glm::mat4x4 ModelMatrix = glm::mat4(1.0);
	glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);
	drawObject((int)drawable, GL_LINES);
}

void drawMesh(Drawable drawable, bool texture = false) {
	glUniform1i(UseTextureID, texture);

	glm::mat4x4 ModelMatrix = glm::mat4(1.0);

	for (int i = 0; i <= (int)drawable; i++) {
		ModelMatrix = glm::scale(ModelMatrix, Scalings[i]);
		ModelMatrix = glm::translate(ModelMatrix, Positions[i]);
		ModelMatrix = glm::rotate(ModelMatrix, Rotations[i].x, { 1, 0, 0 });
		ModelMatrix = glm::rotate(ModelMatrix, Rotations[i].y, { 0, 1, 0 });
		ModelMatrix = glm::rotate(ModelMatrix, Rotations[i].z, { 0, 0, 1 });
	}

	glUniformMatrix4fv(ModelMatrixID, 1, GL_FALSE, &ModelMatrix[0][0]);

	drawObject((int)drawable, GL_TRIANGLES);
}

void drawTexturedMesh(Drawable drawable) {
	// Bind our texture in Texture Unit 0
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, Texture);
	// Set our "myTextureSampler" sampler to use Texture Unit 0
	glUniform1i(TextureID, 0);
	drawMesh(drawable, true);
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
		glUniformMatrix4fv(ViewMatrixID, 1, GL_FALSE, &gViewMatrix[0][0]);
		glUniformMatrix4fv(ProjMatrixID, 1, GL_FALSE, &gProjectionMatrix[0][0]);
		glm::vec3 cameraPos = sphericalToCartesian(CameraRadius, CameraLatitude, CameraLongitude);
		glUniform3f(CameraPosID, cameraPos.x, cameraPos.y, cameraPos.z);

		for (int i = 0; i < NumObjects; i++) {
			Drawable drawable = static_cast<Drawable>(i);
			if (drawable == Drawable::Axis || drawable == Drawable::Grid) { // hardcode these to draw as lines
				drawLines(drawable);
			}
			else if (drawable == Drawable::Head) { // hardcode head
				if (ShowHead) {
					loadObject("finalface.obj", { 0,0,0,0 }, verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
					createVAOs(verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
					drawTexturedMesh(drawable);
				}
				else
				{
					loadObject("finalface2.obj", { 0,0,0,0 }, verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
					createVAOs(verts[(int)Drawable::Head], idcs[(int)Drawable::Head], (int)Drawable::Head);
					drawTexturedMesh(drawable);
				}
				// drawTexturedMesh(drawable);
			}
			else {
				drawTexturedMesh(drawable);
			}
		}

	}
	glUseProgram(0);

	// Swap buffers
	glfwSwapBuffers(window);
	glfwPollEvents();
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
	//glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

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

	// Set up inputs
	glfwSetCursorPos(window, window_width / 2, window_height / 2);
	glfwSetKeyCallback(window, keyCallback);

	return 0;
}

void setCameraPos(void) {
	glm::vec3 cameraPos = sphericalToCartesian(CameraRadius, CameraLatitude, CameraLongitude);
	gViewMatrix = glm::lookAt(glm::vec3(cameraPos),	// eye
		glm::vec3(0.0, 0.0, 0.0),	// center
		glm::vec3(0.0, 1.0, 0.0));	// up
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

	// Camera matrix
	setCameraPos();

	// Create and compile our GLSL program from the shaders
	programID = LoadShaders("StandardShading.vertexshader", "StandardShading.fragmentshader");

	// Get a handle for our "MVP" uniform
	MatrixID = glGetUniformLocation(programID, "MVP");
	ModelMatrixID = glGetUniformLocation(programID, "M");
	ViewMatrixID = glGetUniformLocation(programID, "V");
	ProjMatrixID = glGetUniformLocation(programID, "P");
	CameraPosID = glGetUniformLocation(programID, "cameraPos");

	// Load the texture using any two methods
	Texture = loadBMP_custom("facer.bmp");

	// Get a handle for our "myTextureSampler" uniform
	TextureID = glGetUniformLocation(programID, "myTextureSampler");

	// Get a handle for our "useTexture" uniform
	UseTextureID = glGetUniformLocation(programID, "useTexture");

	createObjects();
}

void createVAOs(Vertex Vertices[], unsigned short Indices[], int ObjectId) {

	GLenum ErrorCheckValue = glGetError();
	const size_t VertexSize = sizeof(Vertices[0]);
	const size_t RgbOffset = sizeof(Vertices[0].Position);
	const size_t NormalOffset = sizeof(Vertices[0].Color) + RgbOffset;
	const size_t UvOffset = sizeof(Vertices[0].Normal) + NormalOffset;

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
	glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)NormalOffset);
	glVertexAttribPointer(3, 2, GL_FLOAT, GL_FALSE, VertexSize, (GLvoid*)UvOffset);

	glEnableVertexAttribArray(0);	// position
	glEnableVertexAttribArray(1);	// color
	glEnableVertexAttribArray(2);	// normal
	glEnableVertexAttribArray(3);	// uv

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

	// Close OpenGL window and terminate GLFW
	glfwTerminate();
}

static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
	// ATTN: MODIFY AS APPROPRIATE
	if (action == GLFW_PRESS) {
		switch (key)
		{
		case GLFW_KEY_R:
		{
			setInitialVals();
			break;
		}
		case GLFW_KEY_F:
		{
			ShowHead = !ShowHead;
			break;
		}
		}
	}
}

void moveCameraPos(double deltaT) {
	double dist = deltaT * CameraMoveSpeed;
	if (glfwGetKey(window, GLFW_KEY_DOWN)) {
		CameraLongitude = (float)min(CameraLongitude + dist, M_PI - 0.0001); // minus a small value to prevent edge case at exactly pi
	}
	if (glfwGetKey(window, GLFW_KEY_UP)) {
		CameraLongitude = (float)max(CameraLongitude - dist, 0 + 0.0001);  // plus a small value to prevent edge case at exactly 0
	}
	if (glfwGetKey(window, GLFW_KEY_RIGHT)) {
		CameraLatitude = (float)fmod(CameraLatitude + dist, M_PI * 2);
	}
	if (glfwGetKey(window, GLFW_KEY_LEFT)) {
		CameraLatitude = (float)fmod(CameraLatitude - dist, M_PI * 2);
	}
}


int main(void)
{
	// initialize window
	int errorCode = initWindow();
	if (errorCode != 0)
		return errorCode;

	// initialize OpenGL pipeline
	initOpenGL();

	setInitialVals();

	// For speed computation
	double lastTime = glfwGetTime();
	do {
		// Measure speed
		double currentTime = glfwGetTime();
		double deltaT = currentTime - lastTime;
		lastTime = currentTime;

		// Move Camera
		moveCameraPos(deltaT);

		// Set Camera
		setCameraPos();

		// DRAWING POINTS
		renderScene();


	} // Check if the ESC key was pressed or the window was closed
	while (glfwGetKey(window, GLFW_KEY_ESCAPE) != GLFW_PRESS &&
		glfwWindowShouldClose(window) == 0);

	cleanup();

	return 0;
}