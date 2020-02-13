// Author : Przemyslaw Zaworski
// Compile: nvcc -o demo.exe demo.cu -lopengl32 -arch=sm_30  user32.lib gdi32.lib
// Usage pattern: demo filename mapsize steps cminr cming cminb cmaxr cmaxg cmaxb intensity threshold crangemin crangemax screenwidth screenheight mode voxelfile offx offy offz
// Example usage: demo obraz3LW_A.vtk 128 128 0.0 0.0 1.0 0.1 0.0 0.0 1.0 1.0 0.0 1.0 1280.0 720.0 0 image.bin 0.25 0.25 0.0

#include <windows.h>
#include <GL/gl.h>
#include <stddef.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <shellapi.h>

// Screen settings
#define FieldOfView 60.0f
#define NearClip 0.01f
#define FarClip 1000.0f
#define VerticalSync 0

// Load selected Win32 API & OpenGL functions
typedef GLuint(WINAPI *PFNGLCREATEPROGRAMPROC) ();
typedef GLuint(WINAPI *PFNGLCREATESHADERPROC) (GLenum t);
typedef void(WINAPI *PFNGLSHADERSOURCEPROC) (GLuint s, GLsizei c, const char*const*string, const GLint* i);
typedef void(WINAPI *PFNGLCOMPILESHADERPROC) (GLuint s);
typedef void(WINAPI *PFNGLATTACHSHADERPROC) (GLuint p, GLuint s);
typedef void(WINAPI *PFNGLLINKPROGRAMPROC) (GLuint p);
typedef void(WINAPI *PFNGLUSEPROGRAMPROC) (GLuint p);
typedef void(WINAPI *PFNGLGENBUFFERSPROC) (GLsizei n, GLuint *b);
typedef void(WINAPI *PFNGLBINDBUFFERPROC) (GLenum t, GLuint b);
typedef void(WINAPI *PFNGLBUFFERDATAPROC) (GLenum t, ptrdiff_t s, const GLvoid *d, GLenum u);
typedef void(WINAPI *PFNGLBINDVERTEXARRAYPROC) (GLuint a);
typedef void(WINAPI *PFNGLENABLEVERTEXATTRIBARRAYPROC) (GLuint i);
typedef void(WINAPI *PFNGLVERTEXATTRIBPOINTERPROC) (GLuint i, GLint s, GLenum t, GLboolean n, GLsizei k, const void *p);
typedef void(WINAPI *PFNGLDISABLEVERTEXATTRIBARRAYPROC) (GLuint i);
typedef int(WINAPI *PFNWGLSWAPINTERVALEXTPROC) (int i);
typedef int(WINAPI *PFNGLGETUNIFORMLOCATIONPROC) (GLuint p, const char *n);
typedef void(WINAPI *PFNGLGENVERTEXARRAYSPROC) (GLsizei n, GLuint *a);
typedef void(WINAPI *PFNGLUNIFORMMATRIX4FVPROC) (GLint l, GLsizei c, GLboolean t, const GLfloat *v);
typedef void(WINAPI *PFNGLUNIFORM1IPROC) (GLint l, GLint v);
typedef void(WINAPI *PFNGLACTIVETEXTUREPROC) (GLenum t);
typedef void(WINAPI *PFNGLUNIFORM3FPROC) (GLint location, float v0, float v1, float v2);
typedef void(WINAPI *PFNGLGETSHADERIVPROC) (GLuint s, GLenum v, GLint *p);
typedef void(WINAPI *PFNGLGETSHADERINFOLOGPROC) (GLuint s, GLsizei b, GLsizei *l, char *i);
typedef void(WINAPI *PFNGLUNIFORM1FPROC) (GLint l, GLfloat v0);
typedef void(WINAPI *PFNGLTEXIMAGE3DPROC) (GLenum a, GLint l, GLint i, GLsizei w, GLsizei h, GLsizei d, GLint b, GLenum f, GLenum t, const void *p);

PFNGLCREATEPROGRAMPROC glCreateProgram;
PFNGLCREATESHADERPROC glCreateShader;
PFNGLSHADERSOURCEPROC glShaderSource;
PFNGLCOMPILESHADERPROC glCompileShader;
PFNGLATTACHSHADERPROC glAttachShader;
PFNGLLINKPROGRAMPROC glLinkProgram;
PFNGLUSEPROGRAMPROC glUseProgram;
PFNGLGENBUFFERSPROC glGenBuffers;
PFNGLBINDBUFFERPROC glBindBuffer;
PFNGLBUFFERDATAPROC glBufferData;
PFNGLBINDVERTEXARRAYPROC glBindVertexArray;
PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray;
PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer;
PFNGLDISABLEVERTEXATTRIBARRAYPROC glDisableVertexAttribArray;
PFNWGLSWAPINTERVALEXTPROC wglSwapIntervalEXT;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation;
PFNGLGENVERTEXARRAYSPROC glGenVertexArrays;
PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv;
PFNGLUNIFORM1IPROC glUniform1i;
PFNGLACTIVETEXTUREPROC glActiveTexture;
PFNGLUNIFORM3FPROC glUniform3f;
PFNGLGETSHADERIVPROC glGetShaderiv;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog;
PFNGLUNIFORM1FPROC glUniform1f;
PFNGLTEXIMAGE3DPROC glTexImage3D;

void glInit()
{
	glCreateProgram = (PFNGLCREATEPROGRAMPROC)wglGetProcAddress("glCreateProgram");
	glCreateShader = (PFNGLCREATESHADERPROC)wglGetProcAddress("glCreateShader");
	glShaderSource = (PFNGLSHADERSOURCEPROC)wglGetProcAddress("glShaderSource");
	glCompileShader = (PFNGLCOMPILESHADERPROC)wglGetProcAddress("glCompileShader");
	glAttachShader = (PFNGLATTACHSHADERPROC)wglGetProcAddress("glAttachShader");
	glLinkProgram = (PFNGLLINKPROGRAMPROC)wglGetProcAddress("glLinkProgram");
	glUseProgram = (PFNGLUSEPROGRAMPROC)wglGetProcAddress("glUseProgram");
	glGenBuffers = (PFNGLGENBUFFERSPROC)wglGetProcAddress("glGenBuffers");
	glBindBuffer = (PFNGLBINDBUFFERPROC)wglGetProcAddress("glBindBuffer");
	glBufferData = (PFNGLBUFFERDATAPROC)wglGetProcAddress("glBufferData");
	glBindVertexArray = (PFNGLBINDVERTEXARRAYPROC)wglGetProcAddress("glBindVertexArray");
	glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC)wglGetProcAddress("glEnableVertexAttribArray");
	glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC)wglGetProcAddress("glVertexAttribPointer");
	glDisableVertexAttribArray = (PFNGLDISABLEVERTEXATTRIBARRAYPROC)wglGetProcAddress("glDisableVertexAttribArray");
	wglSwapIntervalEXT = (PFNWGLSWAPINTERVALEXTPROC)wglGetProcAddress("wglSwapIntervalEXT");
	glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC)wglGetProcAddress("glGetUniformLocation");
	glGenVertexArrays = (PFNGLGENVERTEXARRAYSPROC)wglGetProcAddress("glGenVertexArrays");
	glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC)wglGetProcAddress("glUniformMatrix4fv");
	glUniform1i = (PFNGLUNIFORM1IPROC)wglGetProcAddress("glUniform1i");
	glActiveTexture = (PFNGLACTIVETEXTUREPROC)wglGetProcAddress("glActiveTexture");
	glUniform3f = (PFNGLUNIFORM3FPROC)wglGetProcAddress("glUniform3f");
	glGetShaderiv = (PFNGLGETSHADERIVPROC)wglGetProcAddress("glGetShaderiv");
	glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC)wglGetProcAddress("glGetShaderInfoLog");
	glUniform1f = (PFNGLUNIFORM1FPROC)wglGetProcAddress("glUniform1f");
	glTexImage3D = (PFNGLTEXIMAGE3DPROC) wglGetProcAddress("glTexImage3D");
}

// Global variables
unsigned int VertexBuffer, VertexArrayID;
unsigned char* device;
float* input;
float* output;
float* inputColors;
float* outputColors;	
float offsetX = 0.0f, offsetZ = 0.0f;
float CameraRotYX[4][4], CameraRotYXZ[4][4]; 
float CameraTR[4][4], CameraMatrix[4][4], ViewMatrix[4][4];
float ProjectionViewMatrix[4][4], MVP[4][4];
int cells = 0;
int points = 0;
float sminx = 0.01f;
float smaxx = 0.99f;
float sminy = 0.01f;
float smaxy = 0.99f;
float sminz = 0.01f;
float smaxz = 0.99f;
bool DebugImage = false;

struct DataSet 
{
	float* vertices;
	float* scalars;
};

// 3D coordinates of cube in object space
static const GLfloat vertices[] = 
{
	-0.5f,-0.5f,-0.5f,
	-0.5f,-0.5f, 0.5f,
	-0.5f, 0.5f, 0.5f, 
	 0.5f, 0.5f,-0.5f, 
	-0.5f,-0.5f,-0.5f,
	-0.5f, 0.5f,-0.5f, 
	 0.5f,-0.5f, 0.5f,
	-0.5f,-0.5f,-0.5f,
	 0.5f,-0.5f,-0.5f,
	 0.5f, 0.5f,-0.5f,
	 0.5f,-0.5f,-0.5f,
	-0.5f,-0.5f,-0.5f,
	-0.5f,-0.5f,-0.5f,
	-0.5f, 0.5f, 0.5f,
	-0.5f, 0.5f,-0.5f,
	 0.5f,-0.5f, 0.5f,
	-0.5f,-0.5f, 0.5f,
	-0.5f,-0.5f,-0.5f,
	-0.5f, 0.5f, 0.5f,
	-0.5f,-0.5f, 0.5f,
	 0.5f,-0.5f, 0.5f,
	 0.5f, 0.5f, 0.5f,
	 0.5f,-0.5f,-0.5f,
	 0.5f, 0.5f,-0.5f,
	 0.5f,-0.5f,-0.5f,
	 0.5f, 0.5f, 0.5f,
	 0.5f,-0.5f, 0.5f,
	 0.5f, 0.5f, 0.5f,
	 0.5f, 0.5f,-0.5f,
	-0.5f, 0.5f,-0.5f,
	 0.5f, 0.5f, 0.5f,
	-0.5f, 0.5f,-0.5f,
	-0.5f, 0.5f, 0.5f,
	 0.5f, 0.5f, 0.5f,
	-0.5f, 0.5f, 0.5f,
	 0.5f,-0.5f, 0.5f
};

// Remap value x in range(a,b) to range (c,d)
float remap (float x, float a, float b, float c, float d)  
{
	return (x-a)/(b-a)*(d-c) + c; 
}

// Convert degrees to radians
float deg2rad(float x) 
{
	return (x * 3.14159265358979323846f / 180.0f);
}

// 4x4 matrices multiplication
void Mul(float mat1[][4], float mat2[][4], float res[][4])
{
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			res[i][j] = 0;
			for (int k = 0; k < 4; k++) 
			{
				res[i][j] += mat1[i][k]*mat2[k][j];
			}
		}
	}
}

// Matrix 4x4 inversion
void Inverse( float param[4][4], float k[4][4])
{
	float invOut[16];
	float m[16] = 
	{
		param[0][0],param[0][1],param[0][2],param[0][3],
		param[1][0],param[1][1],param[1][2],param[1][3],
		param[2][0],param[2][1],param[2][2],param[2][3],
		param[3][0],param[3][1],param[3][2],param[3][3]
	};
	float inv[16], det;
	int i;
	inv[0]  =  m[5]*m[10]*m[15]-m[5]*m[11]*m[14]-m[9]*m[6]*m[15]+m[9]*m[7]*m[14]+m[13]*m[6]*m[11]-m[13]*m[7]*m[10];
	inv[4]  = -m[4]*m[10]*m[15]+m[4]*m[11]*m[14]+m[8]*m[6]*m[15]-m[8]*m[7]*m[14]-m[12]*m[6]*m[11]+m[12]*m[7]*m[10];
	inv[8]  =  m[4] *m[9]*m[15]-m[4]*m[11]*m[13]-m[8]*m[5]*m[15]+m[8]*m[7]*m[13]+m[12]*m[5]*m[11]-m[12]*m[7] *m[9];
	inv[12] = -m[4] *m[9]*m[14]+m[4]*m[10]*m[13]+m[8]*m[5]*m[14]-m[8]*m[6]*m[13]-m[12]*m[5]*m[10]+m[12]*m[6] *m[9];
	inv[1]  = -m[1]*m[10]*m[15]+m[1]*m[11]*m[14]+m[9]*m[2]*m[15]-m[9]*m[3]*m[14]-m[13]*m[2]*m[11]+m[13]*m[3]*m[10];
	inv[5]  =  m[0]*m[10]*m[15]-m[0]*m[11]*m[14]-m[8]*m[2]*m[15]+m[8]*m[3]*m[14]+m[12]*m[2]*m[11]-m[12]*m[3]*m[10];
	inv[9]  = -m[0] *m[9]*m[15]+m[0]*m[11]*m[13]+m[8]*m[1]*m[15]-m[8]*m[3]*m[13]-m[12]*m[1]*m[11]+m[12]*m[3] *m[9];
	inv[13] =  m[0] *m[9]*m[14]-m[0]*m[10]*m[13]-m[8]*m[1]*m[14]+m[8]*m[2]*m[13]+m[12]*m[1]*m[10]-m[12]*m[2] *m[9];
	inv[2]  =  m[1] *m[6]*m[15]-m[1] *m[7]*m[14]-m[5]*m[2]*m[15]+m[5]*m[3]*m[14]+m[13]*m[2] *m[7]-m[13]*m[3] *m[6];
	inv[6]  = -m[0] *m[6]*m[15]+m[0] *m[7]*m[14]+m[4]*m[2]*m[15]-m[4]*m[3]*m[14]-m[12]*m[2] *m[7]+m[12]*m[3] *m[6];
	inv[10] =  m[0] *m[5]*m[15]-m[0] *m[7]*m[13]-m[4]*m[1]*m[15]+m[4]*m[3]*m[13]+m[12]*m[1] *m[7]-m[12]*m[3] *m[5];
	inv[14] = -m[0] *m[5]*m[14]+m[0] *m[6]*m[13]+m[4]*m[1]*m[14]-m[4]*m[2]*m[13]-m[12]*m[1] *m[6]+m[12]*m[2] *m[5];
	inv[3]  = -m[1] *m[6]*m[11]+m[1] *m[7]*m[10]+m[5]*m[2]*m[11]-m[5]*m[3]*m[10] -m[9]*m[2] *m[7] +m[9]*m[3] *m[6];
	inv[7]  =  m[0] *m[6]*m[11]-m[0] *m[7]*m[10]-m[4]*m[2]*m[11]+m[4]*m[3]*m[10] +m[8]*m[2] *m[7] -m[8]*m[3] *m[6];
	inv[11] = -m[0] *m[5]*m[11]+m[0] *m[7]*m[9] +m[4]*m[1]*m[11]-m[4]*m[3] *m[9] -m[8]*m[1] *m[7] +m[8]*m[3] *m[5];
	inv[15] =  m[0] *m[5]*m[10]-m[0] *m[6]*m[9] -m[4]*m[1]*m[10]+m[4]*m[2] *m[9] +m[8]*m[1] *m[6] -m[8]*m[2] *m[5];
	det = m[0] * inv[0] + m[1] * inv[4] + m[2] * inv[8] + m[3] * inv[12];
	det = 1.0 / det;
	for (i = 0; i < 16; i++) invOut[i] = inv[i] * det;	
	k[0][0] = invOut[0];  k[0][1] = invOut[1];  k[0][2] = invOut[2];  k[0][3] = invOut[3];
	k[1][0] = invOut[4];  k[1][1] = invOut[5];  k[1][2] = invOut[6];  k[1][3] = invOut[7];
	k[2][0] = invOut[8];  k[2][1] = invOut[9];  k[2][2] = invOut[10]; k[2][3] = invOut[11];
	k[3][0] = invOut[12]; k[3][1] = invOut[13]; k[3][2] = invOut[14]; k[3][3] = invOut[15];  
}

float clamp(float x, float a, float b)
{
	return fmaxf(a, fminf(b, x));
}

// GLSL vertex shader code
static const char* VertexShader = \
	"#version 430 core\n"
	"layout (location=0) in vec3 vertexPosition;"	
	"out vec3 world;"
	"uniform mat4 MVP;"
	"void main()"
	"{"	
		"gl_Position = MVP * vec4(vertexPosition,1.0);"
		"world = vertexPosition;"
	"}";

// GLSL fragment(pixel) shader code	
static const char* FragmentShader = \
	"#version 430 core\n"
	"out vec4 color;"
	"in vec3 world;"
	"uniform int steps;"
	"uniform sampler3D _MainTex;"	
	"uniform vec3 _WorldSpaceCameraPos;"
	"uniform float _Intensity;"
	"uniform float _Threshold;"
	"uniform vec3 _SliceMin ;"
	"uniform vec3 _SliceMax ;"
	"uniform vec3 _ColorMin;"
	"uniform vec3 _ColorMax;"	
	"mat4 _AxisRotationMatrix = mat4(1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,0.0,1.0);"
	"float SampleVolume(vec3 uv, vec3 p)"
	"{"
		"vec3 axis = (   _AxisRotationMatrix * vec4(p, 0) ).xyz + 0.5;"
		"float min = step(_SliceMin.x, axis.x) * step(_SliceMin.y, axis.y) * step(_SliceMin.z, axis.z);"
		"float max = step(axis.x, _SliceMax.x) * step(axis.y, _SliceMax.y) * step(axis.z, _SliceMax.z);"
		"return texture(_MainTex, uv).r * _Intensity * min * max;"
	"}"
	"void main()"
	"{"	
		"vec3 ro = world;"
		"vec3 rd = normalize(world - _WorldSpaceCameraPos);	"		 
		"vec3 AABBmin = vec3(-0.5, -0.5, -0.5);"
		"vec3 AABBmax = vec3(0.5, 0.5, 0.5);"
		"vec3 tbot = (1.0 / rd) * (AABBmin - ro);"
		"vec3 ttop = (1.0 / rd) * (AABBmax - ro);"
		"vec3 tmin = min(ttop, tbot);"
		"vec3 tmax = max(ttop, tbot);"
		"vec2 a = max(tmin.xx, tmin.yz);"
		"float tnear = max(0.0,max(a.x, a.y));"
		"vec2 b = min(tmax.xx, tmax.yz);"
		"float tfar = min(b.x, b.y);"
		"vec3 end = ro + rd * tfar;"
		"vec3 d = normalize(end - ro) * (abs(tfar - tnear) / float(steps));"
		"vec4 t = vec4(0, 0, 0, 0);"
		"for (int i = 0; i < steps; i++)"
		"{"
			"float v = SampleVolume(ro+0.5, ro);"
			"vec4 s = vec4(v, v, v, v);"
			"s.a *= 0.5;"
			"s.rgb *= s.a;"
			"t = (1.0 - t.a) * s + t;"
			"ro += d;"
			"if (t.a > _Threshold) break;"
		"}"
		"color =  clamp(t, 0.0,1.0).rgba;"
		"if (t.r>0.005) "
		"color = vec4(mix(_ColorMin,_ColorMax,t.r),t.a);"
		"else "
		"color = vec4(0.0,0.0,0.0,t.a);"
	"}";

// Compute determinant of 4x4 matrix on the GPU
__device__ float Determinant (float m[4][4]) 
{
	return
	m[0][3] * m[1][2] * m[2][1] * m[3][0] - m[0][2] * m[1][3] * m[2][1] * m[3][0] -
	m[0][3] * m[1][1] * m[2][2] * m[3][0] + m[0][1] * m[1][3] * m[2][2] * m[3][0] +
	m[0][2] * m[1][1] * m[2][3] * m[3][0] - m[0][1] * m[1][2] * m[2][3] * m[3][0] -
	m[0][3] * m[1][2] * m[2][0] * m[3][1] + m[0][2] * m[1][3] * m[2][0] * m[3][1] +
	m[0][3] * m[1][0] * m[2][2] * m[3][1] - m[0][0] * m[1][3] * m[2][2] * m[3][1] -
	m[0][2] * m[1][0] * m[2][3] * m[3][1] + m[0][0] * m[1][2] * m[2][3] * m[3][1] +
	m[0][3] * m[1][1] * m[2][0] * m[3][2] - m[0][1] * m[1][3] * m[2][0] * m[3][2] -
	m[0][3] * m[1][0] * m[2][1] * m[3][2] + m[0][0] * m[1][3] * m[2][1] * m[3][2] +
	m[0][1] * m[1][0] * m[2][3] * m[3][2] - m[0][0] * m[1][1] * m[2][3] * m[3][2] -
	m[0][2] * m[1][1] * m[2][0] * m[3][3] + m[0][1] * m[1][2] * m[2][0] * m[3][3] +
	m[0][2] * m[1][0] * m[2][1] * m[3][3] - m[0][0] * m[1][2] * m[2][1] * m[3][3] -
	m[0][1] * m[1][0] * m[2][2] * m[3][3] + m[0][0] * m[1][1] * m[2][2] * m[3][3];
}

// Compute sign on the GPU
__device__ float Sign (float x)
{
	return ( (0.0f < x) - (x < 0.0f) );
}

// Compute on the GPU whether given point is inside tetrahedron
__device__ bool InsideTetrahedronLegacy (float v1[3], float v2[3], float v3[3], float v4[3], float p[3])
{
	float D0[4][4] = 
	{
		{v1[0],v1[1],v1[2],1.0f},{v2[0],v2[1],v2[2],1.0f},{v3[0],v3[1],v3[2],1.0f},{v4[0],v4[1],v4[2],1.0f}
	};
	float D1[4][4] = 
	{
		{p[0],p[1],p[2],1.0f},{v2[0],v2[1],v2[2],1.0f},{v3[0],v3[1],v3[2],1.0f},{v4[0],v4[1],v4[2],1.0f}
	};
	float D2[4][4] = 
	{
		{v1[0],v1[1],v1[2],1.0f},{p[0],p[1],p[2],1.0f},{v3[0],v3[1],v3[2],1.0f},{v4[0],v4[1],v4[2],1.0f}
	};
	float D3[4][4] = 
	{
		{v1[0],v1[1],v1[2],1.0f},{v2[0],v2[1],v2[2],1.0f},{p[0],p[1],p[2],1.0f},{v4[0],v4[1],v4[2],1.0f}
	};
	float D4[4][4] = 
	{
		{v1[0],v1[1],v1[2],1.0f},{v2[0],v2[1],v2[2],1.0f},{v3[0],v3[1],v3[2],1.0f},{p[0],p[1],p[2],1.0f}
	};
	float a = Determinant(D0);
	float b = Determinant(D1);
	float c = Determinant(D2);
	float d = Determinant(D3);
	float e = Determinant(D4);
	return ( (Sign(a)==Sign(b)) &&  (Sign(a)==Sign(c)) && (Sign(a)==Sign(d)) && (Sign(a)==Sign(e)) );
}

// Compute on the GPU whether given point is inside tetrahedron - faster version	
__device__ bool InsideTetrahedron (float a[3], float b[3], float c[3], float d[3], float p[3])
{
	float vap[3] = {p[0] - a[0], p[1] - a[1], p[2] - a[2]};
	float vbp[3] = {p[0] - b[0], p[1] - b[1], p[2] - b[2]};
	float vab[3] = {b[0] - a[0], b[1] - a[1], b[2] - a[2]};
	float vac[3] = {c[0] - a[0], c[1] - a[1], c[2] - a[2]};
	float vad[3] = {d[0] - a[0], d[1] - a[1], d[2] - a[2]};
	float vbc[3] = {c[0] - b[0], c[1] - b[1], c[2] - b[2]};
	float vbd[3] = {d[0] - b[0], d[1] - b[1], d[2] - b[2]};
	float x[3] = {vbd[1]*vbc[2]-vbd[2]*vbc[1], vbd[2]*vbc[0]-vbd[0]*vbc[2], vbd[0]*vbc[1]-vbd[1]*vbc[0]};
	float y[3] = {vac[1]*vad[2]-vac[2]*vad[1], vac[2]*vad[0]-vac[0]*vad[2], vac[0]*vad[1]-vac[1]*vad[0]};
	float z[3] = {vad[1]*vab[2]-vad[2]*vab[1], vad[2]*vab[0]-vad[0]*vab[2], vad[0]*vab[1]-vad[1]*vab[0]};
	float w[3] = {vab[1]*vac[2]-vab[2]*vac[1], vab[2]*vac[0]-vab[0]*vac[2], vab[0]*vac[1]-vab[1]*vac[0]};		
	float va6 = vbp[0] * x[0] + vbp[1] * x[1] + vbp[2] * x[2];     
	float vb6 = vap[0] * y[0] + vap[1] * y[1] + vap[2] * y[2];   
	float vc6 = vap[0] * z[0] + vap[1] * z[1] + vap[2] * z[2];
	float vd6 = vap[0] * w[0] + vap[1] * w[1] + vap[2] * w[2];
	float q[3] = {vac[1]*vad[2]-vac[2]*vad[1], vac[2]*vad[0]-vac[0]*vad[2], vac[0]*vad[1]-vac[1]*vad[0]};
	float v6 = 1.0f / (vab[0] * q[0] + vab[1] * q[1] + vab[2] * q[2]);
	float k[4] =  {va6*v6, vb6*v6, vc6*v6, vd6*v6};
	return ((k[0] >= 0.0) && (k[0] <= 1.0) && (k[1] >= 0.0) && (k[1] <= 1.0) && (k[2] >= 0.0) && (k[2] <= 1.0) && (k[3] >= 0.0) && (k[3] <= 1.0)) ? true : false;
}

// CUDA kernel, to generate voxel map
__global__ void GenerateVoxelMap (unsigned char* color, float* v, float* s, int k, int msize) 
{
	unsigned int x = blockIdx.x * blockDim.x + threadIdx.x;   
	unsigned int y = blockIdx.y * blockDim.y + threadIdx.y;
	unsigned int z = blockIdx.z * blockDim.z + threadIdx.z;
	unsigned int i = x + y * msize + z * msize * msize ;
	if (color[i]>0) return;    
	float3 resolution = make_float3(msize,msize,msize);
	float3 coordinates = make_float3((float)x, (float)y,float(z));
	float3 uv = make_float3 (coordinates.x / resolution.x, coordinates.y / resolution.y, coordinates.z / resolution.z );
	float a[3] = {v[k*12+0], v[k*12+1], v[k*12+2]};
	float b[3] = {v[k*12+3], v[k*12+4], v[k*12+5]};
	float c[3] = {v[k*12+6], v[k*12+7], v[k*12+8]};
	float d[3] = {v[k*12+9], v[k*12+10], v[k*12+11]};
	float p[3] = {uv.x,uv.y,uv.z};
	if (InsideTetrahedron(a, b, c, d, p))
	{
		float f = (s[k*4+0] + s[k*4+1] + s[k*4+2] + s[k*4+3]) / 4.0f;
		color[i] = (unsigned char)(f*255);
	}
	else 
	{
		color[i] = 0;
	}
}

// Declare matrices for building ModelViewProjection matrix.
float ModelMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f
};

float CameraTranslationMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,-5.0f,
	0.0f,0.0f,0.0f,1.0f
};

float CameraRotationYMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f
};

float CameraRotationXMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f
};

float CameraRotationZMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f
};

float CameraScaleMatrix[4][4] = 
{
	1.0f,0.0f,0.0f,0.0f,
	0.0f,1.0f,0.0f,0.0f,
	0.0f,0.0f,-1.0f,0.0f,
	0.0f,0.0f,0.0f,1.0f
};

float ProjectionMatrix[4][4] = 
{
	0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,0.0f,0.0f,
	0.0f,0.0f,-1.0f,0.0f
};

// Shader debugging
void Debug(int sh)
{
	GLint isCompiled = 0;
	glGetShaderiv(sh,0x8B82,&isCompiled);
	if(isCompiled == GL_FALSE)
	{
		GLint length = 0;
		glGetShaderiv(sh,0x8B84,&length);
		GLsizei q = 0;
		char* log = (char*)malloc(sizeof(char)*length);
		glGetShaderInfoLog(sh,length,&q,log);
		if (length>1)
		{
			FILE *file = fopen ("debug.log","a");
			fprintf (file,"%s\n%s\n",(char*)glGetString(0x8B8C),log);
			fclose (file);
			ExitProcess(0);
		}
		free(log);
	}
}

// Compile shaders	
int MakeShader(const char* VS, const char* FS)
{
	int p = glCreateProgram();
	int s1 = glCreateShader(0x8B31);
	int s2 = glCreateShader(0x8B30);
	glShaderSource(s1,1,&VS,0);
	glShaderSource(s2,1,&FS,0);	
	glCompileShader(s1);
	glCompileShader(s2);
	glAttachShader(p,s1);
	glAttachShader(p,s2);
	glLinkProgram(p);
	Debug(s2);
	return p;
}

// Load data from VTK file.
struct DataSet StreamData (const char* path, float ColorRangeMin, float ColorRangeMax, float ox, float oy, float oz)
{
	struct DataSet grid; 
	FILE* fp = fopen(path, "r");
	char line[200];
	int v=0;
	int m=0;
	int o=0;
	int i=0;
	char *token;
	int state = 0;
	float MinX = 100.0f;
	float MaxX = -100.0f;
	float MinY = 100.0f;
	float MaxY = -100.0f;
	float MinZ = 100.0f;
	float MaxZ = -100.0f;
	float CMin = 100.0f;
	float CMax = -100.0f;
	
	while (1) //read number of cells
	{
		if (fgets(line,150, fp) == NULL) break;
		if(strcmp(line, "\n") == 0) 
		{
			state++;
			continue;
		}
		i++;
		if (i<5) continue;
		if (state==0)
		{
			continue;
		}
		if (state==1)
		{
			token = strtok(line, " ");
			if(strcmp(token, "CELLS") == 0) 
			{
				token = strtok(NULL," ");
				int hh = atoll(token);
				cells = hh;
				points = 4 * hh;
				break;
			}		
		}	
		if (state==2)
		{
			continue;
		}
		if (state==3)
		{
			continue;
		}
	} 

	printf( "Number of cells: %d\n", cells );
	fseek(fp, 0, 0)	;
	v=m=o=i=state=0;
	int p=0;
	
	input = (float*)malloc(4*3*points);
	output = (float*)malloc(4*3*points);
	inputColors = (float*)malloc(4*points);
	outputColors = (float*)malloc(4*points);
	
	while (1) //read vertices and calculate bounds
	{
		if (fgets(line,150, fp) == NULL) break;
		if(strcmp(line, "\n") == 0) 
		{
			state++;
			continue;
		}
		i++;
		if (i<5) continue;
		if (state==0)
		{
			token = strtok(line, " ");
			if(strcmp(token, "POINTS") == 0) continue;
			while( token != NULL ) 
			{
				float h = atof(token);
				input[v] = h;
				v++;
				token = strtok(NULL," ");
			}
		}
		
		if (state==1)
		{
			token = strtok(line, " ");
			if(strcmp(token, "CELLS") == 0) 
			{	
				token = strtok(NULL," ");
				int hh = atoll(token);
				cells = hh;
				points = 4 * hh;
				continue;
			}
			p=0;
			while( token != NULL ) 
			{
				if (p==0) 
				{
					p++;
					token = strtok(NULL," ");
					continue;
				}
				int h = atoll(token);
				output[3*m] = input[3*h];
				output[3*m+1] = input[3*h+1];
				output[3*m+2] = input[3*h+2];
				if (input[3*m]<MinX) MinX = input[3*m];
				if (input[3*m]>MaxX) MaxX = input[3*m];
				if (input[3*m+1]<MinY) MinY = input[3*m+1];
				if (input[3*m+1]>MaxY) MaxY = input[3*m+1];
				if (input[3*m+2]<MinZ) MinZ = input[3*m+2];
				if (input[3*m+2]>MaxZ) MaxZ = input[3*m+2];
				m++;
				token = strtok(NULL," ");
			}		
		}		
		if (state==2)
		{
			continue;
		}
		if (state==3)
		{
			token = strtok(line, " ");
			if(strcmp(token, "POINT_DATA") == 0) continue;
			if(strcmp(token, "SCALARS") == 0) continue;
			if(strcmp(token, "LOOKUP_TABLE") == 0) continue;
			while( token != NULL ) 
			{
				float h = atof(token);
				inputColors[o] = h;
				if (h<CMin) CMin = h;
				if (h>CMax) CMax = h;
				o++;
				token = strtok(NULL," ");
			}
		}
	} 

	printf( "Scalars range: %f %f [remapped to range %f %f]\n", CMin, CMax, ColorRangeMin, ColorRangeMax );
	fseek(fp, 0, 0)	;	
	v=m=o=i=state=p=0;

	float deltaA = abs(MaxX - MinX);
	float deltaB = abs(MaxY - MinY);
	float deltaC = abs(MaxZ - MinZ);
	float MaxDelta = max(deltaA,max(deltaB,deltaC));
	float ScaleX = deltaA/MaxDelta;
	float ScaleY = deltaB/MaxDelta;
	float ScaleZ = deltaC/MaxDelta;
	
	while (1) //set output vertices and colors
	{
		if (fgets(line,150, fp) == NULL) break;
		if(strcmp(line, "\n") == 0) 
		{
			state++;
			continue;
		}
		i++;
		if (i<5) continue;
		if (state==0)
		{
			continue;
		}
		
		if (state==1)
		{
			token = strtok(line, " ");
			if(strcmp(token, "CELLS") == 0) 
			{	
				token = strtok(NULL," ");
				int hh = atoll(token);
				cells = hh;
				points = 4 * hh;
				continue;
			}
			p=0;
			while( token != NULL ) 
			{
				if (p==0) 
				{
					p++;
					token = strtok(NULL," ");
					continue;
				}
				int h = atoll(token);
				output[3*m] = remap(input[3*h],MinX,MaxX,0.0f,ScaleX) + ox;
				output[3*m+1] = remap(input[3*h+1],MinY,MaxY,0.0f,ScaleY) + oy;
				output[3*m+2] = remap(input[3*h+2],MinZ,MaxZ,0.0f,ScaleZ) + oz;
				outputColors[m] = remap(inputColors[h],CMin,CMax,ColorRangeMin,ColorRangeMax);
				m++;
				token = strtok(NULL," ");
			}
		}
		
		if (state==2)
		{
			continue;
		}

		if (state==3)
		{
			continue;
		}
	} 
	
	printf( "Bounds X range: %f %f [remapped to range 0.0 %f]\n", MinX,MaxX,ScaleX );
	printf( "Bounds Y range: %f %f [remapped to range 0.0 %f]\n", MinY,MaxY,ScaleY );
	printf( "Bounds Z range: %f %f [remapped to range 0.0 %f]\n", MinZ,MaxZ,ScaleZ );
	grid.vertices = output;
	grid.scalars = outputColors;
	return grid;
}

// Execute CUDA kernel
unsigned char* LoadCUDATexture(const char* filename, int msize, float CRMin, float CRMax, float ox, float oy, float oz)
{
	unsigned char* host = (unsigned char*) malloc(msize*msize*msize*sizeof(int)); 
	cudaMalloc(&device, msize*msize*msize*sizeof(int));	
	struct DataSet k = StreamData(filename, CRMin, CRMax, ox, oy, oz);
	float* coords;
	float* colors;
	cudaMalloc(&coords, 3*points*sizeof(float));
	cudaMemcpy(coords, k.vertices, 3*points*sizeof(float), cudaMemcpyHostToDevice);
	cudaMalloc(&colors, points*sizeof(float));	
	cudaMemcpy(colors, k.scalars, points*sizeof(float), cudaMemcpyHostToDevice);
	dim3 block(8,8,8); 
	dim3 grid(msize / block.x, msize / block.y, msize / block.z);
	printf( "%s\n", "Voxelizing scene..." );
	for (int j=0; j<cells; j++)
	{
		printf("\rPlease wait... %d percents...", (int)((float)j/(float)cells*100.0f));
		GenerateVoxelMap <<< grid, block >>>(device,coords,colors,j,msize);
	}
	cudaMemcpy(host, device, msize*msize*msize*sizeof(int), cudaMemcpyDeviceToHost);
	cudaFree(coords);
	cudaFree(colors);
	return host;
}

unsigned char* LoadVoxelMapFromFile (const char* filename, int msize)
{
	FILE* binFile = fopen(filename,"rb");
	unsigned char* source = (unsigned char*) malloc(msize*msize*msize*sizeof(int));
	fread(source, sizeof(unsigned char), msize*msize*msize, binFile);
	fclose(binFile);
	return source;
}

void SaveVoxelMapToFile (const char* filename, unsigned char* data, int msize)
{
	FILE* binFile = fopen(filename,"wb");
	fwrite(data,sizeof(unsigned char),msize*msize*msize,binFile);
	fclose(binFile);		
}

// Set 3D volumetric texture in OpenGL environment
void SetTexture(int unit, int id, int shader, const char *name, unsigned char* source,int msize)
{
	glActiveTexture(unit);
	glBindTexture(0x806F, id);	
	glTexImage3D(0x806F,0,0x8229,msize,msize,msize,0,GL_RED,GL_UNSIGNED_BYTE,source);
	glTexParameteri(0x806F,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(0x806F,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	int loc = glGetUniformLocation(shader, name);
	glUniform1i(loc, id);
}

void MouseLook(HWND hwnd, float w, float h)
{	
	if (GetForegroundWindow()!=hwnd) return;
	POINT point;
	int mx = (int)w >> 1;
	int my = (int)h >> 1;
	GetCursorPos(&point);
	if( (point.x == mx) && (point.y == my) ) return;
	SetCursorPos(mx, my);	
	float deltaZ = (float)((mx - point.x)) ;
	float deltaX = (float)((my - point.y)) ;
	if (deltaX>0.0f) offsetX-=0.5f; 
	if (deltaX<0.0f) offsetX+=0.5f; 
	if (deltaZ>0.0f) offsetZ-=0.5f; 
	if (deltaZ<0.0f) offsetZ+=0.5f; 
	CameraRotationXMatrix[1][1] = cos(deg2rad(offsetX));
	CameraRotationXMatrix[1][2] = (-1.0f)*sin(deg2rad(offsetX));
	CameraRotationXMatrix[2][1] = sin(deg2rad(offsetX));
	CameraRotationXMatrix[2][2] = cos(deg2rad(offsetX));				
	CameraRotationYMatrix[0][0] = cos(deg2rad(offsetZ));
	CameraRotationYMatrix[0][2] = sin(deg2rad(offsetZ));
	CameraRotationYMatrix[2][0] = (-1.0f)*sin(deg2rad(offsetZ));
	CameraRotationYMatrix[2][2] = cos(deg2rad(offsetZ));
}
		
void KeyboardMovement(HWND hwnd)
{
	if (GetForegroundWindow()!=hwnd) return;
	float forward[3] = {ViewMatrix[2][0],ViewMatrix[2][1],ViewMatrix[2][2]};
	float strafe[3] = {ViewMatrix[0][0],ViewMatrix[1][0],ViewMatrix[2][0]};
	float dz = 0.0f;
	float dx = 0.0f;
	if (GetAsyncKeyState(0x57)) dz =  2.0f;
	if (GetAsyncKeyState(0x53)) dz = -2.0f ;
	if (GetAsyncKeyState(0x44)) dx =  2.0f;
	if (GetAsyncKeyState(0x41)) dx = -2.0f ;
	if (GetAsyncKeyState(0x45)) CameraTranslationMatrix[1][3] += 0.001f ;
	if (GetAsyncKeyState(0x51)) CameraTranslationMatrix[1][3] -= 0.001f ; 
	float eyeVector[3] = {CameraTranslationMatrix[0][3],CameraTranslationMatrix[1][3] ,CameraTranslationMatrix[2][3]};
	eyeVector[0] += (-dz * forward[0] + dx * strafe[0]) * 0.001f;
	eyeVector[1] += (-dz * forward[1] + dx * strafe[1]) * 0.001f;
	eyeVector[2] += (-dz * forward[2] + dx * strafe[2]) * 0.001f;
	CameraTranslationMatrix[0][3] = eyeVector[0];
	CameraTranslationMatrix[1][3] = eyeVector[1];
	CameraTranslationMatrix[2][3] = eyeVector[2];
	DebugImage = false;
	if (GetAsyncKeyState(0x58) & 0x8000) DebugImage = true;
}

// Slice volumetric texture
void Slicer (int smin, int smax)
{
	if (GetAsyncKeyState(0x55)) {sminx -= 0.001f; clamp(sminx,0.0f,1.0f);} //u
	if (GetAsyncKeyState(0x49)) {sminx += 0.001f; clamp(sminx,0.0f,1.0f);} //i
	if (GetAsyncKeyState(0x4F)) {smaxx -= 0.001f; clamp(smaxx,0.0f,1.0f);} //o
	if (GetAsyncKeyState(0x50)) {smaxx += 0.001f; clamp(smaxx,0.0f,1.0f);} //p
	if (GetAsyncKeyState(0x48)) {sminy -= 0.001f; clamp(sminy,0.0f,1.0f);} //h
	if (GetAsyncKeyState(0x4A)) {sminy += 0.001f; clamp(sminy,0.0f,1.0f);} //j
	if (GetAsyncKeyState(0x4B)) {smaxy -= 0.001f; clamp(smaxy,0.0f,1.0f);} //k
	if (GetAsyncKeyState(0x4C)) {smaxy += 0.001f; clamp(smaxy,0.0f,1.0f);} //l
	if (GetAsyncKeyState(0x56)) {sminz -= 0.001f; clamp(sminz,0.0f,1.0f);} //v
	if (GetAsyncKeyState(0x42)) {sminz += 0.001f; clamp(sminz,0.0f,1.0f);} //b
	if (GetAsyncKeyState(0x4E)) {smaxz -= 0.001f; clamp(smaxz,0.0f,1.0f);} //n
	if (GetAsyncKeyState(0x4D)) {smaxz += 0.001f; clamp(smaxz,0.0f,1.0f);} //m
	glUniform3f(smin, sminx, sminy, sminz);
	glUniform3f(smax, smaxx, smaxy, smaxz);	
}

// Diagnostic
void PrintDeviceName()
{
	cudaDeviceProp device;
	cudaGetDeviceProperties(&device, 0);
	printf( "%s\n", device.name );
	printf( "CUDA compute capability: %d\n", device.major );
}

void Release ()
{
	free (input);
	free (output);
	free (inputColors);
	free (outputColors);
	cudaFree(device);
}

static LRESULT CALLBACK WindowProc( HWND hWnd, UINT uMsg, WPARAM wParam, LPARAM lParam )
{
	if( uMsg==WM_SYSCOMMAND && (wParam==SC_SCREENSAVE || wParam==SC_MONITORPOWER) )
		return 0;
	if( uMsg==WM_CLOSE || uMsg==WM_DESTROY || (uMsg==WM_KEYDOWN && wParam==VK_ESCAPE) )
	{
		PostQuitMessage(0);
		return 0;
	}
	if( uMsg==WM_SIZE )
	{
		glViewport( 0, 0, lParam&65535, lParam>>16 );
	}
	if( uMsg==WM_CHAR || uMsg==WM_KEYDOWN)
	{
		if( wParam==VK_ESCAPE )
		{
			PostQuitMessage(0);
			return 0;
		}
	}
	return(DefWindowProc(hWnd,uMsg,wParam,lParam));
}

// Main function
int WINAPI WinMain( HINSTANCE hInstance, HINSTANCE hPrevInstance, LPSTR lpCmdLine, int nCmdShow )
{
	AllocConsole();
	freopen("CONOUT$", "w+", stdout);
	int argc;
	char** argv;
	LPWSTR* lpArgv = CommandLineToArgvW( GetCommandLineW(), &argc );
	argv = (char**)malloc( argc*sizeof(char*) );
	int size = 0;
	for( int i=0; i < argc; ++i )
	{
		size = wcslen( lpArgv[i] ) + 1;
		argv[i] = (char*)malloc( size );
		wcstombs( argv[i], lpArgv[i], size );
	}
	if (argc < 21)
	{
		printf( "Not enough arguments. Exit..." );
		ExitProcess(0);
	}
	LocalFree(lpArgv);
	PrintDeviceName();
	float ScreenWidth = atof(argv[14]);
	float ScreenHeight = atof(argv[15]);
	unsigned char* Source;
	if (atoi(argv[16])==0) Source = LoadCUDATexture(argv[1],atoi(argv[2]),atof(argv[12]),atof(argv[13]),atof(argv[18]),atof(argv[19]),atof(argv[20]));
	if (atoi(argv[16])==1) Source = LoadVoxelMapFromFile(argv[17],atoi(argv[2])),atoi(argv[2]);
	if (atoi(argv[16])==2) 
	{	
		SaveVoxelMapToFile(argv[17], LoadCUDATexture(argv[1],atoi(argv[2]),atof(argv[12]),atof(argv[13]),atof(argv[18]),atof(argv[19]),atof(argv[20])), atoi(argv[2]));
		printf( "Save To file. Exit..." );
		ExitProcess(0);
	}	
	MSG msg;
	int exit = 0;
	ShowCursor(0);	
	PIXELFORMATDESCRIPTOR pfd = { 0,0,PFD_DOUBLEBUFFER };
	WNDCLASS wc;
	ZeroMemory( &wc, sizeof(WNDCLASS) );
	wc.style = CS_OWNDC|CS_HREDRAW|CS_VREDRAW;
	wc.lpfnWndProc = WindowProc;
	wc.hInstance = 0;
	wc.lpszClassName = "demo";
	wc.hbrBackground =(HBRUSH)CreateSolidBrush(0x00102030);
	RegisterClass(&wc);
	HWND hwnd = CreateWindowEx(0, wc.lpszClassName, "Demo", WS_VISIBLE|WS_OVERLAPPEDWINDOW, 0, 0, ScreenWidth, ScreenHeight, 0, 0, 0, 0);
	HDC hdc = GetDC(hwnd);
	SetPixelFormat(hdc,ChoosePixelFormat(hdc,&pfd),&pfd);
	wglMakeCurrent(hdc,wglCreateContext(hdc));
	glInit();
	wglSwapIntervalEXT (VerticalSync);
	glGenVertexArrays (1, &VertexArrayID);
	glBindVertexArray (VertexArrayID);	
	glGenBuffers(1, &VertexBuffer);
	glBindBuffer(0x8892, VertexBuffer);
	glBufferData(0x8892, sizeof(vertices), vertices, 0x88E4);	
	int PS = MakeShader(VertexShader,FragmentShader);
	SetTexture(0x84C0, 0, PS, "_MainTex", Source, atoi(argv[2]));
	glUseProgram(PS);	
	int MatrixID = glGetUniformLocation(PS,"MVP"); 
	int WorldSpaceID = glGetUniformLocation(PS,"_WorldSpaceCameraPos");
	int SliceMin = glGetUniformLocation(PS,"_SliceMin");
	int SliceMax = glGetUniformLocation(PS,"_SliceMax");
	int StepsID = glGetUniformLocation(PS,"steps");	
	int ColorMin = glGetUniformLocation(PS,"_ColorMin");
	glUniform3f(ColorMin, atof(argv[4]),atof(argv[5]),atof(argv[6]));	
	int ColorMax = glGetUniformLocation(PS,"_ColorMax");
	glUniform3f(ColorMax, atof(argv[7]),atof(argv[8]),atof(argv[9]));
	int Intensity = glGetUniformLocation(PS,"_Intensity");
	glUniform1f(Intensity, atof(argv[10]));
	int Threshold = glGetUniformLocation(PS,"_Threshold");
	glUniform1f(Threshold, atof(argv[11]));	
	ProjectionMatrix[0][0] = ((1.0f/tan(deg2rad(FieldOfView/2.0f)))/(ScreenWidth/ScreenHeight));
	ProjectionMatrix[1][1] = (1.0f/tan(deg2rad(FieldOfView/2.0f)));
	ProjectionMatrix[2][2] = (-1.0f)* (FarClip+NearClip)/(FarClip-NearClip);
	ProjectionMatrix[2][3] = (-1.0f)*(2.0f*FarClip*NearClip)/(FarClip-NearClip)	;
	while( !exit )
	{
		if (DebugImage)
		{
			glDisable (GL_BLEND);
		}
		else
		{
			glEnable (GL_BLEND); 
			glBlendFunc (GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
		}
		while(PeekMessage(&msg, 0, 0, 0, PM_REMOVE) )
		{
			if( msg.message==WM_QUIT ) exit = 1;
			TranslateMessage( &msg );
			DispatchMessage( &msg );
		}
		MouseLook(hwnd,ScreenWidth,ScreenHeight);
		Slicer(SliceMin,SliceMax);
		Mul(CameraRotationYMatrix,CameraRotationXMatrix,CameraRotYX);
		Mul(CameraRotYX,CameraRotationZMatrix,CameraRotYXZ);	
		Mul(CameraTranslationMatrix,CameraRotYXZ,CameraTR);
		Mul(CameraTR,CameraScaleMatrix,CameraMatrix);
		Inverse(CameraMatrix,ViewMatrix);
		Mul(ProjectionMatrix,ViewMatrix,ProjectionViewMatrix);
		Mul(ProjectionViewMatrix,ModelMatrix,MVP);	
		float MVPT[4][4] = 
		{
			MVP[0][0], MVP[1][0], MVP[2][0], MVP[3][0],
			MVP[0][1], MVP[1][1], MVP[2][1], MVP[3][1],
			MVP[0][2], MVP[1][2], MVP[2][2], MVP[3][2],
			MVP[0][3], MVP[1][3], MVP[2][3], MVP[3][3]
		};
		glEnable(GL_DEPTH_TEST);
		glDepthFunc(GL_LESS);		
		glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
		KeyboardMovement(hwnd);
		if (DebugImage) 
			glClearColor(0.0f, 0.0f, 1.0f, 1.0f);
		else
			glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
		glUniform1i(StepsID, atoi(argv[3]) );
		glUniformMatrix4fv(MatrixID, 1, GL_FALSE, &MVPT[0][0]);
		glUniform3f(WorldSpaceID, CameraTranslationMatrix[0][3], CameraTranslationMatrix[1][3], CameraTranslationMatrix[2][3]);
		glEnableVertexAttribArray(0);
		glBindBuffer(0x8892, VertexBuffer);
		glVertexAttribPointer(0,3, GL_FLOAT, GL_FALSE, 0,(void*)0 );
		glDrawArrays(GL_TRIANGLES, 0, 12*3);
		glDisableVertexAttribArray(0);
		wglSwapLayerBuffers(hdc, WGL_SWAP_MAIN_PLANE);
	}
	Release ();
	free(Source);
	free(argv);
	return 0;
}