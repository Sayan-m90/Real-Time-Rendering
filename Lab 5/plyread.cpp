////////////////////////////////////////////////////////
//
// ply model reading sample program
//
// Han-Wei Shen
//
////////////////////////////////////////////////////////

#include <stdio.h>
#include <cstdlib>
#include <stack> 
#include <vector>
#include <math.h> 

#ifdef __APPLE__
#include <GLUT/glut.h> 
#include <OpenGL/gl.h> 
#else
#include <GL/glew.h> 
#include <GL/glut.h> 
#include <GL/gl.h> 
#endif
#include <time.h>
#include<Windows.h>
#include<glm/glm.hpp>
#include<glm/gtx/transform.hpp>

#include "ply.h"

typedef struct Vertex {
  float x, y, z;
  float nx, ny, nz;
} Vertex;

typedef struct
{
	float location[4];
	float normal[4];
	float color[4];
	float texCoord[2];
}myVertex;


myVertex sqverts[4]; 
GLuint sqindices[6]; 
int sq_nindices = 6; 

///////////////////////////////////////////////////////////////////////////
#define SetMaterialColor(d, r, g, b, a)  glUniform4f(d, r, g, b, a); 
//////////////////////////////////////////////////////////////////////////

GLuint programObject;
GLuint SetupGLSL(char*);
//myVertex sqverts[4];
//GLuint sqindices[6];
//int sq_nindices = 6;

typedef struct Face {
  unsigned int count;
  unsigned int *vertices;
  float nx, ny, nz;
} Face;

char* string_list[] = {
  "x", "y", "z", "nx", "ny", "nz", "vertex_indices"
};

#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

float cx, cy, cz; 
float x_min, x_max, y_min, y_max, z_min, z_max; 

Vertex** pl_vertices = 0;
Face** pl_faces = 0;
myVertex *pl_verts;
GLuint *pl_indices;
unsigned int pl_vertexcount;
unsigned int pl_facecount;
int pl_vertexnormals = 0;
int pl_facenormals = 0;
int pl_nindices;


Vertex** pl_vertices2 = 0;
Face** pl_faces2 = 0;
myVertex *pl_verts2;
GLuint *pl_indices2;
unsigned int pl_vertexcount2;
unsigned int pl_facecount2;
int pl_vertexnormals2 = 0;
int pl_facenormals2 = 0;
int pl_nindices2;


Vertex** pl_vertices3 = 0;
Face** pl_faces3 = 0;
myVertex *pl_verts3;
GLuint *pl_indices3;
unsigned int pl_vertexcount3;
unsigned int pl_facecount3;
int pl_vertexnormals3 = 0;
int pl_facenormals3 = 0;
int pl_nindices3;

/////////////////////////////////////////////////////// VBO

char fname1[] = "negx.ppm";
char fname2[] = "negy.ppm";
char fname3[] = "negz.ppm";
char fname4[] = "posx.ppm";
char fname5[] = "posy.ppm";
char fname6[] = "posz.ppm";

GLuint pl_vboHandle[1];   // a VBO that contains interleaved positions and cube_colors 
GLuint pl_indexVBO; 
GLuint pl_vboHandle2[1];   // a VBO that contains interleaved positions and cube_colors 
GLuint pl_indexVBO2; 
GLuint pl_vboHandle3[1];   // a VBO that contains interleaved positions and cube_colors 
GLuint pl_indexVBO3; 


GLuint cubemap_texture;
/////////////////////////////////////////////////////////////

int xform_mode_cyl = 0;
int xform_mode_cube = 0;
int press_x_cyl, press_y_cyl; 
int release_x_cyl, release_y_cyl; 
int press_x_cube, press_y_cube; 
int release_x_cube, release_y_cube; 
float x_angle_rot = 0.0;		//rotate whole scene by mouse
float x_angle_body = 0.0;		// rotate body by keyboard
float z_angle = 0.0;
float x_angle = 0.0;
float scale_size = 1;
bool  x_angle_body_flag = true;
float y_angle_rot = 0.0; 
float scale_size_cyl = 1; 
bool samescale = true;

//////////////////////////////////////////////////////////////////////////////////
/////////////////// CYlinder		/////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

typedef struct 
{
  float cyl_location[4]; 
  float normal[4]; 
  float cyl_color[4]; 
} Vertexcyl; 

int nindices; 
Vertexcyl *cyverts;   // cylinder vertices
GLuint *cindices; 

GLuint cyl_vboHandle[1];   // a VBO that contains interleaved positions and cyl_colors 
GLuint cyl_indexVBO; 


bool rotate_x = false;
bool rotate_y = false;
void mymouse(int, int, int, int);
void mymotion(int, int);
void mykey(unsigned char, int, int);
float ply_min[3],ply_max[3];
float ply_min2[3],ply_max2[3];
float ply_min3[3],ply_max3[3];
glm::mat4 modelM = glm::mat4(1.0f);

///////////////////////////////////////////////////////

float vertices[] = {
	-0.5, -0.5, -0.5, 1.0, //location				1
	0, 0, -1.0, 0,			//normal
	0, 0.6, 0, 1.0,			//color
	0.5,0.5, -0.5,1.0,		//location
	0,0,-1,0,				//normal				2
	0,0.6,0,1,				//color
	
	0.5, -0.5, -0.5, 1.0, //location			3
	0, 0, -1.0, 0,			//normal			
	0, 0.6, 0, 1.0,			//color
	-0.5, -0.5, -0.5,1.0,		//location			4
	0,0,-1,0,				//normal
	0,0.6, 0, 1,				//color
	
	-0.5, 0.5, -0.5, 1.0, //location				5
	0, 0, -1.0, 0,			//normal
	0, 0.6, 0, 1.0,			//color
	0.5,0.5, -0.5,1.0,		//location				6
	0,0,-1,0,				//normal
	0,0.6,0,1,				//color
	
	-0.5, -0.5, 0.5, 1.0, //location				7
	0, 0, 1.0, 0,			//normal
	0, 0.6, 0, 1.0,			//color
	0.5,-0.5, 0.5,1.0,		//location							8
	0,0,1,0,				//normal
	0,0.6,0,1,				//color
	
	0.5, 0.5, 0.5, 1.0, //location					9
	0, 0, 1.0, 0,			//normal
	0, 0, 0.6, 1.0,			//color
	-0.5,-0.5, 0.5,1.0,		//location						10
	0,0,1,0,				//normal
	0,0.6,0,1,				//color
	
	0.5, 0.5, 0.5, 1.0, //location					11
	0, 0, 1.0, 0,			//normal					
	0, 0, 0.6, 1.0,			//color
	-0.5,0.5, 0.5,1.0,		//location				12
	0,0,1,0,				//normal
	0,0,0.6,1,				//color

	
	0.5, -0.5, -0.5, 1.0, //location					13
	1.0, 0, 0, 0,			//normal
	0, 0, 0.6, 1.0,			//color
	0.5,0.5, -0.5,1.0,		//location					14
	1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	0.5, 0.5, 0.5, 1.0, //location						15
	1, 0, 0, 0,			//normal							
	0, 0, 0.6, 1.0,			//color
	0.5,-0.5, -0.5,1.0,		//location					16
	1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	0.5, 00.5, 00.5, 1.0, //location					17
	1, 0, 0.0, 0,			//normal
	0, 0, 0.6, 1.0,			//color							
	0.5,-0.5, 0.5,1.0,		//location						18
	1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	-0.5, -0.5, -0.5, 1.0, //location					19
	-1.0, 0, 0, 0,			//normal							
	0, 0, 0.6, 1.0,			//color
	-0.5, 0.5, 0.5,1.0,		//location					20
	-1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	-0.5, 0.5, -0.5, 1.0, //location				21
	-1, 0, 0, 0,			//normal					
	0, 0, 0.6, 1.0,			//color
	-0.5, -0.5, -0.5,1.0,		//location			22
	-1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	-0.5, -0.5, 00.5, 1.0, //location				23
	-1.0,  0.0, 0, 0,			//normal
	0, 0, 0.6, 1.0,			//color					
	-0.5,0.5, 0.5,1.0,		//location			24
	-1,0,0,0,				//normal
	0,0,0.6,1,				//color
	
	-0.5, -0.5, -0.5, 1.0, //location			25
	0, -1.0, 0, 0,			//normal
	0, 0, 0.6, 1.0,			//color
	0.5,-0.5, -0.5,1.0,		//location				26
	0,-1,0,0,				//normal
	0,0,0.6,1,				//color
	
	0.5, -0.5, 0.5, 1.0, //location					27
	0, -1.0, 0, 0,			//normal
	0, 0, 0.6, 1.0,			//color						
	-0.5,-0.5, -0.5,1.0,		//location			28
	0,-1,0,0,				//normal
	0,0,0.6,1,				//color
	
	0.5, -0.5, 0.5, 1.0, //location				29
	0, -1.0, 0, 0,			//normal
	0, 0, 0.6, 1.0,			//color
	-0.5, -0.5, 0.5,1.0,		//location			30
	0, -1, 0, 0,				//normal
	0,0, 0.6,1,				//color
	
	-0.5, 0.5, -0.5, 1.0, //location			31
	0, 1.0, 0, 0,			//normal				
	0, 0, 0.6, 1.0,			//color
	0.5,0.5, 0.5,1.0,		//location		32
	0,1,0,0,				//normal
	0,0,0.6,1,				//color
	
	0.5, 0.5, -0.5, 1.0, //location			33
	0, 1.0, 0, 0,			//normal
	0, 0,0.6, 1.0,			//color				
	-0.5,0.5, -0.5,1.0,		//location		34
	0,1,0,0,				//normal
	0,0,0.6,1,				//color
	
	-0.5, 0.5, 0.5, 1.0, //location			35
	0, 1.0, 0, 0,			//normal
	0, 0, 0.6,1.0,			//color
	0.5,0.5, 0.5,1.0,		//location			36
	0,1,0,0,				//normal
	0,0, 0.6,1				//color

};


typedef struct 
{
  float sph_location[4]; 
  float sph_normal[4]; 
  float sph_color[4]; 
} Vertexsphere; 

int sphindices; 
Vertexsphere *sphverts;   // sphere vertices
GLuint *sp_indices; 

GLuint sphere_vboHandle[1];   // a VBO that contains interleaved positions and sph_colors 
GLuint sphere_indexVBO; 

///////////////////////////////////////////////////////////////
// Define Light Properties -  Ia, Id, Is, and light position 
////////////////////////////////////////////////////////////////

GLfloat light_ambient[4] = {0.8,0.8,0.8,1};  //Ia 
GLfloat light_diffuse[4] = {0.4,0.8,0.8,1};  //Id
GLfloat light_specular[4] = {1,1,1,1};  //Is
GLfloat light_pos [4] = {-0.0995,-5.089,-3.00994,1};

/////////////////////////////////////////////////////////////////
// Define Default Material Properties -  Ka, Kd, Ks, Shininess 
////////////////////////////////////////////////////////////////

GLfloat mat_ambient[4] = {0.2,0.2,0.2,1};  //Ka 
GLfloat mat_diffuse[4] = {0.8,0.8,0.8,1};  //Kd
GLfloat mat_specular[4] = {1,1,1,1};  //Ks
GLfloat mat_shine[1] = {100}; 

////////////////////////////////////////////////
typedef struct 
{
  float location[4]; 
  float normal[4]; 
 float color[4]; 
  float texCoord[2]; 
} Vertexcube; 
//
Vertexcube cverts[24];       // cube vertices 


GLubyte tindices[36];   // cube face indices: 6 faces, 6 indices each face 

GLuint square_vboHandle[1];   // a VBO that contains interleaved positions and square_colors 
GLuint square_indexVBO; 

GLuint cube_vboHandle[1];
GLuint cube_indexVBO;
////////////////////////////////////////////////////////////////////////////////




float M_PI = 3.14127;
glm::mat4 model = glm::mat4(1.0f); 
glm::mat4 model1 = glm::mat4(1.0f); 
glm::mat4 model2 = glm::mat4(1.0f); 
glm::mat4 model3 = glm::mat4(1.0f); 
/////////////////////////////////
// glut mouse control 
// 

std::stack<glm::mat4> st;

float r_upper_arm = 0.0f;
float r_lower_arm = 0.0f;
float l_upper_arm = 0.0f;
float l_lower_arm = 0.0f;
float r_lower_arm_x = -2.1;
float r_lower_arm_y =  1.0;
float l_lower_arm_x = 2.1;
float l_lower_arm_y =  1.0;
float hand_move = 0.0;
float lower_hand_move = 0.0;
bool hand_move_flag = false;
bool hand_move_forward = true;
int times,timed;

bool leg_move_flag = false;
bool leg_move_forward = true;
float leg_move = true;

bool r_arm_up = true;
bool l_arm_up = true;


float head = 0.0f;
bool head_flag = false;
//float headl = 0.0f;


bool draw_line = false; 
int which_sph_color = 0; // 0: T  1: N  2: B

////////////////////////////////////////////////////////////////////////////////////
//
//  cross product 
//
void cross(float a0, float a1, float a2, 
		   float b0, float b1, float b2, 
		   float*a_cross_b)
{
  a_cross_b[0] = a1*b2-b1*a2; 
  a_cross_b[1] = a2*b0-b2*a0; 
  a_cross_b[2] = a0*b1-a1*b0; 
}

////////////////////////////////////////////////////////////////////////////////////
///////////////////	Read Image		///////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////
GLubyte checkImage[256][256][4];
GLubyte readImage[600][600][4]; 
GLubyte texImage[512][512][4]; 
GLubyte gradientImage[512][512][4]; 

bool use_gradient = false; 
void read_Image(char* fname) 
{
  //  FILE* in = fopen("normal.ppm", "r"); 
  //  FILE* in = fopen("rough.ppm", "r"); 
  //  FILE* in = fopen("brick.ppm", "r"); 
  //  FILE* in = fopen("earth.ppm", "r"); 

    FILE* in = fopen(fname, "r"); 

  int height, width, ccv; 
  char header[100]; 
  fscanf(in, "%s %d %d %d", header, &width, &height, &ccv); 

  printf("%s %d %d %d\n", header, width, height, ccv);
  int r, g, b; 

  for (int i=height-1; i>=0; i--)
     for (int j=0; j<width; j++)
{
      fscanf(in, "%d %d %d", &r, &g, &b); 
      readImage[i][j][0] = (GLubyte)r; 
      readImage[i][j][1] = (GLubyte)g; 
      readImage[i][j][2] = (GLubyte)b; 
      readImage[i][j][3] = 255; 
    }

  for (int i=0; i<512; i++)
    for ( int j=0; j<512; j++) {
      if (i<height && j <width) {
	texImage[i][j][0] = readImage[i][j][0]; texImage[i][j][1] = readImage[i][j][1];
	texImage[i][j][2] = readImage[i][j][2];	texImage[i][j][3] = 255; 
      }
      else {
      	texImage[i][j][0] = 0; 	texImage[i][j][1] = 0; 	texImage[i][j][2] = 0; 
	texImage[i][j][3] = 255; 
      }
    }
  //////////////////////////////////
  // Use central difference to calculate the gradients 
  for (int i=0; i<512; i++)
    for ( int j=0; j<512; j++) {
      gradientImage[i][j][0] = (unsigned char)((texImage[(i+1)%512][j][0]-
						    texImage[(i-1)%512][j][0])/2.0+256); 
      gradientImage[i][j][1] = (unsigned char)((texImage[i][(j+1)%512][0]-
				     texImage[i][(j-1)%512][0])/2.0+256); 
      gradientImage[i][j][2] = 20;  
      gradientImage[i][j][3] = 255; 
      }
  fclose(in); 
}

void InitSphere(int nslices, int nstacks, float r, float g, float b) 
{

  float istep = 2 * 3.14159/(float)(nslices); 
  float jstep = 2 * 3.14159 /(float)(nstacks); 

  int nvertices = nslices * nstacks; 
  sphverts = new Vertexsphere[nvertices]; 

  for (int j =0; j<nstacks; j++)
    for (int i=0; i<nslices; i++) {
      int idx = j*nslices + i; // mesh[j][i] 
	  jstep =  j *2* 3.14159/(float)(nstacks);
	  istep =  i *2* 3.14159/(float)(nslices);
      float s = istep; 
      float t = jstep; 
      float x = 4*cos(s)*sin(t);//3*cos(s)+cos(t)*cos(s); 
      float y = 4*sin(t)*sin(s);//3*cos(s)+cos(t)*sin(s); 
	  float z = 4*cos(t);//sin(t); 
	 
      sphverts[idx].sph_location[0] = x; 
      sphverts[idx].sph_location[1] = y; 
      sphverts[idx].sph_location[2] = z; 

      float tangentx = 4*cos(s)*cos(t);//-3*sin(s)-cos(t)*sin(s); 
      float tangenty =  4*cos(s)*sin(t);//3*cos(s)+cos(t)*cos(s); 
      float tangentz =  4*cos(s); 
	  float dx_dt = -4*sin(t)*sin(s);//-sin(t)*cos(s); 
	  float dy_dt = 4*cos(t)*sin(s);//-sin(t)*sin(s); 
      float dz_dt = 0; 

      float sph_normal[3]; 
      float bisph_normal[3]; 

      cross(tangentx, tangenty, tangentz, dx_dt, dy_dt, dz_dt, 
	    sph_normal); 
      cross(sph_normal[0], sph_normal[1], sph_normal[2], tangentx, tangenty, 
	    tangentz, bisph_normal); 

      sphverts[idx].sph_normal[0] = (sph_normal[0]); 
      sphverts[idx].sph_normal[1] = (sph_normal[1]); 
      sphverts[idx].sph_normal[2] = (sph_normal[2]); 

      sphverts[idx].sph_location[3] = 1.0;  sphverts[idx].sph_normal[3] = 0.0; 

      if (which_sph_color == 0) {
      float mag = sqrt(tangentx*tangentx+tangenty*tangenty+tangentz*tangentz); 
      sphverts[idx].sph_color[0] = tangentx/mag; 
      sphverts[idx].sph_color[1] = tangenty/mag; 
      sphverts[idx].sph_color[2] = tangentz/mag; 
      sphverts[idx].sph_color[3] = 1.0; 
      }
      if (which_sph_color == 1) {
      float mag = sqrt(sph_normal[0]*sph_normal[0]+sph_normal[1]*sph_normal[1]+sph_normal[2]*sph_normal[2]); 
      sphverts[idx].sph_color[0] = sph_normal[0]/mag; 
      sphverts[idx].sph_color[1] = sph_normal[1]/mag; 
      sphverts[idx].sph_color[2] = sph_normal[2]/mag; 
      sphverts[idx].sph_color[3] = 1.0; 
      }
      if (which_sph_color == 2) {
      float mag = sqrt(bisph_normal[0]*bisph_normal[0]+bisph_normal[1]*bisph_normal[1]+bisph_normal[2]*bisph_normal[2]); 
      sphverts[idx].sph_color[0] = bisph_normal[0]/mag; 
      sphverts[idx].sph_color[1] = bisph_normal[1]/mag; 
      sphverts[idx].sph_color[2] = bisph_normal[2]/mag; 
      sphverts[idx].sph_color[3] = 1.0; 
      }
    }

  // now create the index array 
  sphindices = (nstacks)*(nslices)*6; 
  sp_indices = new GLuint[sphindices]; 
  int n = 0; 
  for (int j =0; j<nstacks; j++)
    for (int i=0; i<nslices; i++) {
      
      int i1 = i; 
      int i2 = (i+1) % nslices; 
      int j1 = j; 
      int j2 = (j+1) % nstacks; 

      int idx1 = j1* nslices + i1; 
      int idx2 = j1* nslices + i2; 
      int idx3 = j2* nslices + i2; 
      int idx4 = j2* nslices + i1; 

      sp_indices[n++] = idx1; 
      sp_indices[n++] = idx2; 
      sp_indices[n++] = idx4; 

      sp_indices[n++] = idx2; 
      sp_indices[n++] = idx3; 
      sp_indices[n++] = idx4;  
    }
}

////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////
//
// create VBO objects and send the triangle vertices/sph_colors to the graphics card
// 
void InitVBO_sphere(int nslices, int nstacks) 
{
  int nvertices = nslices * nstacks; 
  sphindices = (nstacks)*6*(nslices); 

  glGenBuffers(1, &sphere_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, sphere_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexsphere)*nvertices, sphverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &sphere_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*sphindices, sp_indices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

  // by now, we moved the position and sph_color data over to the graphics card. There will be no redundant data copy at drawing time 
} 


void InitCylinder(int nslices, int nstacks, float r, float g, float b) 
{
  int nvertices = nslices * nstacks; 
  cyverts = new Vertexcyl[nvertices]; 

 
  float Dangle = 2*M_PI/(float)(nslices-1); 
  
  for (int j =0; j<nstacks; j++)
  {
    for (int i=0; i<nslices; i++) {
      int idx = j*nslices + i; // mesh[j][i] 
      float angle = Dangle * i; 
	  //printf("angle:%f",angle);
      cyverts[idx].cyl_location[0] = cyverts[idx].normal[0] = cos(angle); 
      cyverts[idx].cyl_location[1] = cyverts[idx].normal[1] = sin(angle); 
      cyverts[idx].cyl_location[2] = j*1.0/(float)(nstacks-1); 
      cyverts[idx].normal[2] = 0.0; 
      cyverts[idx].cyl_location[3] = 1.0;  cyverts[idx].normal[3] = 0.0; 
      
	  cyverts[idx].cyl_color[0] = r;
	  cyverts[idx].cyl_color[1] = g;
	  cyverts[idx].cyl_color[2] = b; 
      cyverts[idx].cyl_color[3] = 1.0; 
    }
  if(r>0.1)		r-=0.02;
  if(g<0.9)		g+=0.02;
  if(b>=0.1)	b-=0.02;
  }
   //now create the index array 

  nindices = (nstacks-1)*2*(nslices+1); 
  cindices = new GLuint[nindices]; 
  int n = 0; 
  for (int j =0; j<nstacks-1; j++)
    for (int i=0; i<=nslices; i++) {
      int mi = i % nslices;  
      int idx = j*nslices + mi; // mesh[j][mi] 
      int idx2 = (j+1) * nslices + mi; 
      cindices[n++] = idx; 
      cindices[n++] = idx2; 
    }
}

void InitVBO_cyl(int nslices, int nstacks) 
{
  InitCylinder(nslices, nstacks, 1.0, 1.0, 0.0);
	int nvertices = nslices * nstacks; 
  nindices = (nstacks-1)*2*(nslices+1); 

  glGenBuffers(1, &cyl_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cyl_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexcyl)*nvertices, cyverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &cyl_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cyl_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*nindices, cindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

   //by now, we moved the position and cyl_color data over to the graphics card. There will be no redundant data copy at drawing time 
}

//////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////
//
// create VBO objects and send the triangle vertices/square_colors to the graphics card
// 
void Initsquare()
{ float size = 1.0;
	
  sqverts[0].location[0] = -0.5; sqverts[0].location[1] = -0.5;  sqverts[0].location[2] = 0;  sqverts[0].location[3] = 1; 
  sqverts[0].normal[0] = 0; sqverts[0].normal[1] = 0;  sqverts[0].normal[2] = 1;  sqverts[0].normal[3] = 0; 
  sqverts[0].texCoord[0] = 0; sqverts[0].texCoord[1] = 0; 
  sqverts[0].color[0] = 0.5; sqverts[0].color[1] = 0.5;  sqverts[0].color[2] = 1;  sqverts[0].color[3] = 1; 

  sqverts[1].location[0] = 0.5; sqverts[1].location[1] = -0.5;  sqverts[1].location[2] = 0;  sqverts[1].location[3] = 1; 
  sqverts[1].normal[0] = 0; sqverts[1].normal[1] = 0;  sqverts[1].normal[2] = 1;  sqverts[1].normal[3] = 0; 
  sqverts[1].texCoord[0] = 1.0; sqverts[1].texCoord[1] = 0; 
  sqverts[1].color[0] = 0.5; sqverts[1].color[1] = 0.5;  sqverts[1].color[2] = 1;  sqverts[1].color[3] = 1; 

  sqverts[2].location[0] = 0.5; sqverts[2].location[1] = 0.5;  sqverts[2].location[2] = 0;  sqverts[2].location[3] = 1; 
  sqverts[2].normal[0] = 0; sqverts[2].normal[1] = 0;  sqverts[2].normal[2] = 1;  sqverts[2].normal[3] = 0; 
  sqverts[2].texCoord[0] = 1.0; sqverts[2].texCoord[1] = 1.0; 
  sqverts[2].color[0] = 0.5; sqverts[2].color[1] = 0.5;  sqverts[2].color[2] = 1;  sqverts[2].color[3] = 1; 

  sqverts[3].location[0] = -0.5; sqverts[3].location[1] = 0.5;  sqverts[3].location[2] = 0;  sqverts[3].location[3] = 1; 
  sqverts[3].normal[0] = 0; sqverts[3].normal[1] = 0;  sqverts[3].normal[2] = 1;  sqverts[3].normal[3] = 0; 
  sqverts[3].texCoord[0] = 0.0; sqverts[3].texCoord[1] = 1.0; 
  sqverts[3].color[0] =1; sqverts[3].color[1] = 1;  sqverts[3].color[2] = 1;  sqverts[3].color[3] = 1; 


  sqindices[0] = 0;   sqindices[1] = 1;   sqindices[2] = 2; 
  sqindices[3] = 2;   sqindices[4] = 3;   sqindices[5] = 0; 


}

void InitVBO_square() 
{
		Initsquare();
  glGenBuffers(1, &square_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, square_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(myVertex)*4, sqverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &square_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, square_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*6, sqindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up */


 /* glGenBuffers(1, &square_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexcube)*24, cverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &cube_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*24, tindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 
  */

} 

void InitCube_VBO() 
{
		//Initsquare();
 glGenBuffers(1, &cube_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*432, vertices, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 
  int j=0;
  for(j=0; j<36 ;j++)
	  tindices[j] = j;
  glGenBuffers(1, &cube_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLubyte)*36, tindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up */


} 

 

bool first = true;

///////////////////////////////////////////////////////////////////////
//
void store_ply(PlyFile* input, Vertex ***pl_vertices, Face ***pl_faces,
	       unsigned int* pl_vertexcount, unsigned int* pl_facecount,
	       int* pl_vertexnormals, int* pl_facenormals) {
  int i, j;

  // go through the element types
  for(i = 0; i < input->num_elem_types; i = i + 1) {
    int count;
    
    // setup the element for reading and get the element count
    char* element = setup_element_read_ply(input, i, &count);

    // get vertices
    if(strcmp("vertex", element) == 0) {
      *pl_vertices = (Vertex**)malloc(sizeof(Vertex) * count);
      *pl_vertexcount = count;

      // run through the properties and store them
      for(j = 0; j < input->elems[i]->nprops; j = j + 1) {
	PlyProperty* property = input->elems[i]->props[j];
	PlyProperty setup;

	if(strcmp("x", property->name) == 0 &&
	   property->is_list == PLY_SCALAR) {

	  setup.name = string_list[0];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, x);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	}
	else if(strcmp("y", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[1];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, y);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	}
	else if(strcmp("z", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[2];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, z);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	}
	else if(strcmp("nx", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[3];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, nx);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_vertexnormals = 1;
	}
	else if(strcmp("ny", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[4];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, ny);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_vertexnormals = 1;
	}
	else if(strcmp("nz", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[5];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Vertex, nz);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_vertexnormals = 1;
	}
	// dunno what it is
	else {
	  fprintf(stderr, "unknown property type found in 1%s: %s\n",
		  element, property->name);
	}
      }

      // do this if you want to grab the other data
      // list_pointer = get_other_properties_ply
      //                (input, offsetof(Vertex, struct_pointer));

      // copy the data
      for(j = 0; j < count; j = j + 1) {
	(*pl_vertices)[j] = (Vertex*)malloc(sizeof(Vertex));
	
	get_element_ply(input, (void*)((*pl_vertices)[j]));
      }
    }
    // get faces
    else if(strcmp("face", element) == 0) {
      *pl_faces = (Face**)malloc(sizeof(Face) * count);
      *pl_facecount = count;

      // run through the properties and store them
      for(j = 0; j < input->elems[i]->nprops; j = j + 1) {
	PlyProperty* property = input->elems[i]->props[j];
	PlyProperty setup;

	if(strcmp("vertex_indices", property->name) == 0 &&
	   property->is_list == PLY_LIST) {

	  setup.name = string_list[6];
	  setup.internal_type = Uint32;
	  setup.offset = offsetof(Face, vertices);
	  setup.count_internal = Uint32;
	  setup.count_offset = offsetof(Face, count);

	  setup_property_ply(input, &setup);
	}
	else if(strcmp("nx", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[3];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Face, nx);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_facenormals = 1;
	}
	else if(strcmp("ny", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[4];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Face, ny);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_facenormals = 1;
	}
	else if(strcmp("nz", property->name) == 0 &&
		property->is_list == PLY_SCALAR) {

	  setup.name = string_list[5];
	  setup.internal_type = Float32;
	  setup.offset = offsetof(Face, nz);
	  setup.count_internal = 0;
	  setup.count_offset = 0;

	  setup_property_ply(input, &setup);
	  *pl_facenormals = 1;
	}
	// dunno what it is
	else {
	  fprintf(stderr, "unknown property type found in %s: %s\n",
		  element, property->name);
	}
      }
	
      // do this if you want to grab the other data
      // list_pointer = get_other_properties_ply
      //                (input, offsetof(Face, struct_pointer));
      
      // copy the data
      for(j = 0; j < count; j = j + 1) {
	(*pl_faces)[j] = (Face*)malloc(sizeof(Face));
	
	get_element_ply(input, (void*)((*pl_faces)[j]));
      }
    }
    // who knows?
    else {
      fprintf(stderr, "unknown element type found: %s\n", element);
    }
  }

  // if you want to grab the other data do this
  // get_other_element_ply(input);
}

//////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////

void InitPLY_VBO(float min[3], float max[3])
{
	int pl_nvertices = pl_vertexcount;
	pl_nindices = pl_facecount * 3;
	pl_verts = new myVertex[pl_nvertices];
	pl_indices = new GLuint[pl_nindices];
	min[0] = min[1] = min[2] = 9999999999999;
	max[0] = max[1] = max[2] = -999999999999;
	for(int i=0;i<pl_vertexcount;i++)
	{
		pl_verts[i].location[0] = pl_vertices[i]->x;
		pl_verts[i].location[1] = pl_vertices[i]->y;
		pl_verts[i].location[2] = pl_vertices[i]->z;
		pl_verts[i].location[3] = 1.0;
		if(pl_verts[i].location[0] < min[0]) min[0] = pl_verts[i].location[0];
		if(pl_verts[i].location[1] < min[1]) min[1] = pl_verts[i].location[1];
		if(pl_verts[i].location[2] < min[2]) min[2] = pl_verts[i].location[2];

		if(pl_verts[i].location[0] > max[0]) max[0] = pl_verts[i].location[0];
		if(pl_verts[i].location[1] > max[1]) max[1] = pl_verts[i].location[1];
		if(pl_verts[i].location[2] > max[2]) max[2] = pl_verts[i].location[2];

		pl_verts[i].normal[0] = pl_vertices[i]->nx;
		pl_verts[i].normal[1] = pl_vertices[i]->ny;
		pl_verts[i].normal[2] = pl_vertices[i]->nz;
		pl_verts[i].normal[3] = 0.0;

		pl_verts[i].color[0] = 1.0;
		pl_verts[i].color[1] = 1.0;
		pl_verts[i].color[2] = 1.0;
		pl_verts[i].color[3] = 1.0;
	}
	for(int i=0;i<pl_facecount;i++)
	{
		//if(faces[i]->count != 3) printf("No\n");
		pl_indices[i*3] = pl_faces[i]->vertices[0];
		pl_indices[i*3+1] = pl_faces[i]->vertices[1];
		pl_indices[i*3+2] = pl_faces[i]->vertices[2];
	}
	glGenBuffers(1,pl_vboHandle);
	glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(myVertex)*pl_nvertices, pl_verts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &pl_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*pl_nindices, pl_indices, GL_STATIC_DRAW);  // load the index data 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
}


void InitPLY_VBO2(float min2[3], float max2[3])
{
	int pl_nvertices2 = pl_vertexcount2;
	pl_nindices2 = pl_facecount2 * 3;
	pl_verts2 = new myVertex[pl_nvertices2];
	pl_indices2 = new GLuint[pl_nindices2];
	min2[0] = min2[1] = min2[2] = 9999999999999;
	max2[0] = max2[1] = max2[2] = -999999999999;
	for(int i=0;i<pl_vertexcount2;i++)
	{
		pl_verts2[i].location[0] = pl_vertices2[i]->x;
		pl_verts2[i].location[1] = pl_vertices2[i]->y;
		pl_verts2[i].location[2] = pl_vertices2[i]->z;
		pl_verts2[i].location[3] = 1.0;
		if(pl_verts2[i].location[0] < min2[0]) min2[0] = pl_verts2[i].location[0];
		if(pl_verts2[i].location[1] < min2[1]) min2[1] = pl_verts2[i].location[1];
		if(pl_verts2[i].location[2] < min2[2]) min2[2] = pl_verts2[i].location[2];

		if(pl_verts2[i].location[0] > max2[0]) max2[0] = pl_verts2[i].location[0];
		if(pl_verts2[i].location[1] > max2[1]) max2[1] = pl_verts2[i].location[1];
		if(pl_verts2[i].location[2] > max2[2]) max2[2] = pl_verts2[i].location[2];

		pl_verts2[i].normal[0] = pl_vertices2[i]->nx;
		pl_verts2[i].normal[1] = pl_vertices2[i]->ny;
		pl_verts2[i].normal[2] = pl_vertices2[i]->nz;
		pl_verts2[i].normal[3] = 0.0;

		pl_verts2[i].color[0] = 1.0;
		pl_verts2[i].color[1] = 1.0;
		pl_verts2[i].color[2] = 1.0;
		pl_verts2[i].color[3] = 1.0;
	}
	for(int i=0;i<pl_facecount2;i++)
	{
		//if(faces[i]->count != 3) printf("No\n");
		pl_indices2[i*3] = pl_faces2[i]->vertices[0];
		pl_indices2[i*3+1] = pl_faces2[i]->vertices[1];
		pl_indices2[i*3+2] = pl_faces2[i]->vertices[2];
	}
	glGenBuffers(1,pl_vboHandle2);
	glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle2[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(myVertex)*pl_nvertices2, pl_verts2, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &pl_indexVBO2); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO2); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*pl_nindices2, pl_indices2, GL_STATIC_DRAW);  // load the index data 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
}

void InitPLY_VBO3(float min3[3], float max3[3])
{
	int pl_nvertices3 = pl_vertexcount3;
	pl_nindices3 = pl_facecount3 * 3;
	pl_verts3 = new myVertex[pl_nvertices3];
	pl_indices3 = new GLuint[pl_nindices3];
	min3[0] = min3[1] = min3[2] = 9999999999999;
	max3[0] = max3[1] = max3[2] = -999999999999;
	for(int i=0;i<pl_vertexcount3;i++)
	{
		pl_verts3[i].location[0] = pl_vertices3[i]->x;
		pl_verts3[i].location[1] = pl_vertices3[i]->y;
		pl_verts3[i].location[2] = pl_vertices3[i]->z;
		pl_verts3[i].location[3] = 1.0;
		if(pl_verts3[i].location[0] < min3[0]) min3[0] = pl_verts3[i].location[0];
		if(pl_verts3[i].location[1] < min3[1]) min3[1] = pl_verts3[i].location[1];
		if(pl_verts3[i].location[2] < min3[2]) min3[2] = pl_verts3[i].location[2];

		if(pl_verts3[i].location[0] > max3[0]) max3[0] = pl_verts3[i].location[0];
		if(pl_verts3[i].location[1] > max3[1]) max3[1] = pl_verts3[i].location[1];
		if(pl_verts3[i].location[2] > max3[2]) max3[2] = pl_verts3[i].location[2];

		pl_verts3[i].normal[0] = pl_vertices3[i]->nx;
		pl_verts3[i].normal[1] = pl_vertices3[i]->ny;
		pl_verts3[i].normal[2] = pl_vertices3[i]->nz;
		pl_verts3[i].normal[3] = 0.0;

		pl_verts3[i].color[0] = 1.0;
		pl_verts3[i].color[1] = 1.0;
		pl_verts3[i].color[2] = 1.0;
		pl_verts3[i].color[3] = 1.0;
	}
	for(int i=0;i<pl_facecount3;i++)
	{
		//if(faces[i]->count != 3) printf("No\n");
		pl_indices3[i*3] = pl_faces3[i]->vertices[0];
		pl_indices3[i*3+1] = pl_faces3[i]->vertices[1];
		pl_indices3[i*3+2] = pl_faces3[i]->vertices[2];
	}
	glGenBuffers(1,pl_vboHandle3);
	glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle3[0]);
	glBufferData(GL_ARRAY_BUFFER, sizeof(myVertex)*pl_nvertices3, pl_verts3, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &pl_indexVBO3); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO3); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*pl_nindices3, pl_indices3, GL_STATIC_DRAW);  // load the index data 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);
}

void draw_sphere(glm::mat4 local2clip, glm::mat4 local2eye, float* world2eye, float color[3], GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4) 
{
    
  glBindBuffer(GL_ARRAY_BUFFER, sphere_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_indexVBO);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Vertexsphere),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(Vertexsphere),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Vertexsphere),(char*) NULL+16); // normal

  glm::mat4 normal_matrix = glm::inverse(local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) &local2clip[0][0]);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) &local2eye[0][0]);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix 
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  
  glDrawElements(GL_TRIANGLE_STRIP, sphindices, GL_UNSIGNED_INT, (char*) NULL+0); 

}


void draw_square(float* local2clip, float* local2eye, float* world2eye, float color[3],  GLuint c0, GLuint c1, GLuint c2, GLuint c3, GLuint m1, GLuint m2, GLuint m3, GLuint m4)
{

  glBindBuffer(GL_ARRAY_BUFFER, square_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, square_indexVBO);
  glEnableVertexAttribArray(c3);
  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+16); // normal 
  glVertexAttribPointer(c3,2,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+48); // texture

  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix

  glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, (char*)NULL+0); 

}

void draw_cube(float* local2clip, float* local2eye, float* world2eye, float color[3],  GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4, GLuint m5)
{

  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO);
//  glEnableVertexAttribArray(c3);
  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(float)*12,(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(float)*12,(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(float)*12,(char*) NULL+16); // normal 
  //glVertexAttribPointer(c3,2,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+48); // texture

  glm::mat4 eye2world = glm::inverse((glm::mat4)*world2eye);
  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m5, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*)NULL+0); 

}
void draw_cylinder(glm::mat4 local2clip, glm::mat4 local2eye, float* world2eye, float color[3], GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4) 
{
    
  glBindBuffer(GL_ARRAY_BUFFER, cyl_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cyl_indexVBO);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Vertexcyl),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(Vertexcyl),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Vertexcyl),(char*) NULL+16); // normal

  glm::mat4 normal_matrix = glm::inverse(local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) &local2clip[0][0]);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) &local2eye[0][0]);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix 
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 

  glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0); 

}
////////////////////////////////////////////////////////////////////////////////////
void Init_texture(int tex) {

  GLuint renderTex; 
  GLuint gradientTex; 

  /*char fname1[] = "negx.ppm";
char fname2[] = "negy.ppm";
char fname3[] = "negz.ppm";
char fname4[] = "posx.ppm";
char fname5[] = "posy.ppm";
char fname6[] = "posz.ppm";
*/

  read_Image(fname1); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE0); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	//////////////////	SECOND TEXTURE
   read_Image(fname2); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE3); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
  

	//////////////////	THIRD TEXTURE
   read_Image(fname3); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE4); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	//////////////////	FOURTH TEXTURE
   read_Image(fname4); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE5); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	////////////////ea//	FIFTH TEXTURE
   read_Image(fname5); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE6); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 

	//////////////////	SIXTH TEXTURE
   read_Image(fname6); 

  glGenTextures(1, &renderTex); 
  glActiveTexture(GL_TEXTURE7); 
  glBindTexture(GL_TEXTURE_2D, renderTex); 

  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, 512, 
	            512,0, GL_RGBA, GL_UNSIGNED_BYTE, 
				texImage);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR); 
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR); 
}

///////////////////////////////////////////////////////////////////////////////
void init_cubemap()
{

	glEnable(GL_TEXTURE_CUBE_MAP);
	glGenTextures(1,&cubemap_texture);
	printf("cubemap_texture %d \n",&cubemap_texture);
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_CUBE_MAP,cubemap_texture);
	glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_S,GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_T,GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_WRAP_R,GL_REPEAT);
	glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_CUBE_MAP,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT,GL_NICEST);
	glTexEnvf(GL_TEXTURE_ENV,GL_TEXTURE_ENV_MODE,GL_REPLACE);
	for(int i=0;i<6;i++)
	{
		glTexImage2D(GL_TEXTURE_CUBE_MAP_POSITIVE_X+i,0,GL_RGBA,512,512,0,GL_RGBA,GL_UNSIGNED_BYTE,texImage);
	}
}
void draw_ply(float* local2clip, float *local2eye, float *world2eye, float color[3], GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4, GLuint m5)
{
	glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+16); // normal

  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glm::mat4 eye2world = glm::inverse((glm::mat4)*world2eye);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix 
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  glUniformMatrix4fv(m5, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  
  glDrawElements(GL_TRIANGLES, pl_nindices, GL_UNSIGNED_INT, (char*) NULL+0); 
}

void draw_ply2(float* local2clip, float *local2eye, float *world2eye, float color[3], GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4, GLuint m5)
{
  glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle2[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO2);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+16); // normal

  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glm::mat4 eye2world = glm::inverse((glm::mat4)*world2eye);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix 
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  glUniformMatrix4fv(m5, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  
  glDrawElements(GL_TRIANGLES, pl_nindices2, GL_UNSIGNED_INT, (char*) NULL+0); 
}


void draw_ply3(float* local2clip, float *local2eye, float *world2eye, float color[3], GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4, GLuint m5)
{
  glBindBuffer(GL_ARRAY_BUFFER, pl_vboHandle3[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, pl_indexVBO3);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+0);  // position 
  glVertexAttribPointer(c1,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+32); // color
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(myVertex),(char*) NULL+16); // normal

  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glm::mat4 eye2world = glm::inverse((glm::mat4)*world2eye);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix 
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  glUniformMatrix4fv(m5, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix 
  
  glDrawElements(GL_TRIANGLES, pl_nindices3, GL_UNSIGNED_INT, (char*) NULL+0); 
}

///////////////////////////////////////////////////////////////

void display()
{
	glClearColor(0,0,0,1); 
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
	float color[3];
	  timed = time(NULL);
 if(timed-times >= 5)
	 {
		float temp;
		//temp = 0.2;
		temp = light_diffuse[0];
		 light_diffuse[0] = light_diffuse[1];  //Id
		 light_diffuse[1] = temp;
		 times=time(NULL);
	 }
	glPolygonMode(GL_FRONT_AND_BACK,GL_FILL);
   glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects 
   glUseProgram(programObject);
   glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 
  /////////////////////////////////////////////
  ////////////////////////    Shader 
  /////////////////////////////////////////////



  GLuint c0 = glGetAttribLocation(programObject, "position");
  GLuint c1 = glGetAttribLocation(programObject, "color");
  GLuint c2 = glGetAttribLocation(programObject, "normal");
  GLuint c3 = glGetAttribLocation(programObject, "texCoord");

  GLuint m1 = glGetUniformLocation(programObject, "local2clip");
  GLuint m2 = glGetUniformLocation(programObject, "local2eye");
  GLuint m3 = glGetUniformLocation(programObject, "normal_matrix");
  GLuint m4 = glGetUniformLocation(programObject, "world2eye");  
  GLuint m5 = glGetUniformLocation(programObject, "eye2world");  

  GLuint Ia = glGetUniformLocation(programObject, "light_ambient");
  GLuint Id = glGetUniformLocation(programObject, "light_diffuse");
  GLuint Is = glGetUniformLocation(programObject, "light_specular");
  GLuint Lpos = glGetUniformLocation(programObject, "light_pos");

  GLuint Ka = glGetUniformLocation(programObject, "mat_ambient");
  GLuint Kd = glGetUniformLocation(programObject, "mat_diffuse");
  GLuint Ks = glGetUniformLocation(programObject, "mat_specular");
  GLuint Shine = glGetUniformLocation(programObject, "mat_shine"); 

  GLuint b1= glGetUniformLocation(programObject,"use_texture");
  int use_texture = 0, tex_loc;
  glUniform1i(b1,use_texture);
  glUniform4f(Ia, light_ambient[0], light_ambient[1], light_ambient[2], light_ambient[3]);
  glUniform4f(Id, light_diffuse[0], light_diffuse[1], light_diffuse[2], light_diffuse[3]);
  glUniform4f(Is, light_specular[0], light_specular[1], light_specular[2], light_specular[3]);
  glUniform4f(Lpos, light_pos[0], light_pos[1], light_pos[2], light_pos[3]);

  glUniform4f(Ka, mat_ambient[0], mat_ambient[1], mat_ambient[2], mat_ambient[3]);
  glUniform4f(Kd, mat_diffuse[0], mat_diffuse[1], mat_diffuse[2], mat_diffuse[3]);
  glUniform4f(Ks, mat_specular[0], mat_specular[1], mat_specular[2], mat_specular[3]);
  glUniform1f(Shine, mat_shine[0]); 
  
  glEnableVertexAttribArray(c0);
  glEnableVertexAttribArray(c1);
  glEnableVertexAttribArray(c2);
  //glEnableVertexAttribArray(c3);
    /////////////////////////////////////////////
  ///////////////////////////   End of Shader
  ////////////////////////////////////////////



  ///////////////////////////////////////////////////////////////
  ///////////////// 3D veiw and projection	/////////////////////
  /////////////////////////////////////////////////////////////
  float angle = 90.0*3.141/180;
  glm::mat4 mvp;
  glm::mat4 mv;
  float toradian = 100*3.141/180;
  glm::mat4 projection = glm::perspective(toradian,1.0f,.1f,100.0f); 
 

  glm::mat4 view = glm::lookAt(glm::vec3(5, 5, 0.2), 
			       glm::vec3(0.0, 0.0, 0.0), 
			       glm::vec3(0.0, 0.0, 1.0)); 

  //glm::mat4 model = glm::mat4(1.0f);

  if(rotate_x == true) 
  {
	  modelM = glm::rotate(modelM, x_angle_rot, glm::vec3(0.0f, 0.0f, 1.0f));
	  x_angle_rot = 0;
	  rotate_x = false;

  }

  if(rotate_y == true) 
	{modelM = glm::rotate(modelM, y_angle_rot, glm::vec3(-1.0f, 0.0f, 0.0f)); 

	 rotate_y = false;
	 y_angle_rot = 0.0f;
    }
   
   if(samescale==false) 
	  { modelM = glm::scale(modelM, glm::vec3(scale_size_cyl, scale_size_cyl, scale_size_cyl)); 

		samescale = true;
		}

	
////////////////////////draw a light
  
  glm::mat4 lightM = modelM; 
  lightM = glm::translate(modelM, glm::vec3(light_pos[0], light_pos[1], light_pos[2])); 
  //lightM = glm::scale(lightM, glm::vec3(0.05, 0.05, 0.05));

 mvp = projection*view*lightM;
 mv = view*lightM;

  
  SetMaterialColor(Kd, 1, 1, 1, 1); 
  draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3, m4);
  

  ////////////////////////////////////////////////
  /////////////	draw ply model	/////////////////
  //////////////////////////////////////////////
  
 
  float ply_center[3], ply_size[3];
  float ply_center2[3], ply_size2[3];
  float ply_center3[3], ply_size3[3];
  ply_center[0] = (ply_min[0]+ply_max[0])/2.0;
  ply_center[1] = (ply_min[1]+ply_max[1])/2.0;
  ply_center[2] = (ply_min[0]+ply_max[2])/2.0;
  ply_size[0] = (ply_max[0] - ply_min[0]);
  ply_size[1] = (ply_max[1] - ply_min[1]);
  ply_size[2] = (ply_max[2] - ply_min[2]);


  ply_center2[0] = (ply_min2[0]+ply_max2[0])/2.0;
  ply_center2[1] = (ply_min2[1]+ply_max2[1])/2.0;
  ply_center2[2] = (ply_min2[0]+ply_max2[2])/2.0;
  ply_size2[0] = (ply_max2[0] - ply_min2[0]);
  ply_size2[1] = (ply_max2[1] - ply_min2[1]);
  ply_size2[2] = (ply_max2[2] - ply_min2[2]);


  ply_center3[0] = (ply_min3[0]+ply_max3[0])/2.0;
  ply_center3[1] = (ply_min3[1]+ply_max3[1])/2.0;
  ply_center3[2] = (ply_min3[0]+ply_max3[2])/2.0;
  ply_size3[0] = (ply_max3[0] - ply_min3[0]);
  ply_size3[1] = (ply_max3[1] - ply_min3[1]);
  ply_size3[2] = (ply_max3[2] - ply_min3[2]);

  float torad;
  torad= 90.0*3.141/180;
  glm::mat4 plyM;
  ///////////////////////////		BUDDHA		/////////////////////////////////
  float rad = (200.0*3.141/180);
  use_texture = 0;
 	glUniform1i(b1,use_texture);
    tex_loc = glGetUniformLocation(programObject,"cubeMap");
    glUniform1i(tex_loc,1);
  plyM = glm::rotate(modelM,torad,glm::vec3(1.0,0.0,0.0f));
  plyM = glm::translate(plyM,glm::vec3(-1*ply_center[0], -1*ply_center[1], -1*ply_center[2]));
  plyM = glm::translate(plyM,glm::vec3(0.9,-8.5,-3.5f));
  plyM = glm::scale(plyM,glm::vec3(2,2,2));
  plyM = glm::translate(plyM,glm::vec3(-0.75,-0.80,4.01f));
  plyM = glm::scale(plyM,glm::vec3(30,30,30));
  plyM = glm::rotate(plyM,rad,glm::vec3(0.0,1.0,0.0));
  //plyM = glm::rotate(plyM,rad,glm::vec3(0,1,0));
  mvp = projection*view*plyM;
  mv = view*plyM;
  SetMaterialColor(Kd,0.96,0.77,0.17,1);
  draw_ply(&mvp[0][0],&mv[0][0],&view[0][0],color,c0,c1,c2,m1,m2,m3,m4,m5);





  ///////////////COPTER
  use_texture = 2;
  glUniform1i(b1,use_texture);
  tex_loc = glGetUniformLocation(programObject,"cubeMap");
  glUniform1i(tex_loc,1);
  //plyM = glm::rotate(modelM,torad,glm::vec3(1.0,0.0,0.0f));
  plyM = glm::translate(modelM,glm::vec3(-1*ply_center2[0], -1*ply_center2[1], -1*ply_center2[2]));
  plyM = glm::translate(plyM,glm::vec3(4.0,10.7,4.0f));
  plyM = glm::scale(plyM,glm::vec3(5/(ply_size2[0]),5/(ply_size2[1]),5/(ply_size2[2])));
  
  float ang = 90.0;
  SetMaterialColor(Kd,0.6,0.3,0.3,1);
  plyM = plyM ;//* modelM;
  mvp = projection *view*plyM;
  mv = view*plyM;
  draw_ply2(&mvp[0][0],&mv[0][0],&view[0][0],color,c0,c1,c2,m1,m2,m3,m4,m5);
  

    ///////////////DRAGON
  use_texture = 2;
  glUniform1i(b1,use_texture);
  tex_loc = glGetUniformLocation(programObject,"cubeMap");
  glUniform1i(tex_loc,1);
  float angle2= 110*3.14/180;
  plyM = glm::translate(modelM,glm::vec3(-1*ply_center3[0], -1*ply_center3[1], -1*ply_center3[2]));
  plyM = glm::translate(plyM,glm::vec3(9,-4.0,-14.5f));
  plyM = glm::rotate(plyM,angle2,glm::vec3(1.0,0.0,0.0));
  angle2= 10*3.14/180;
 plyM = glm::rotate(plyM,angle2,glm::vec3(0.0,1.0,0.0));
  plyM = glm::scale(plyM,glm::vec3(10/(ply_size3[0]),3/(ply_size3[1]),3/(ply_size3[2]))); 
  SetMaterialColor(Kd,0.8,0.2,0.2,1);
  plyM = plyM ;
  mvp = projection *view*plyM;
  mv = view*plyM;
  draw_ply3(&mvp[0][0],&mv[0][0],&view[0][0],color,c0,c1,c2,m1,m2,m3,m4,m5);

      ///////////////DRAGON 2
  use_texture = 2;
  glUniform1i(b1,use_texture);
  tex_loc = glGetUniformLocation(programObject,"cubeMap");
  glUniform1i(tex_loc,1);
  angle2= 110*3.14/180;
  plyM = glm::translate(modelM,glm::vec3(-1*ply_center3[0], -1*ply_center3[1], -1*ply_center3[2]));
  plyM = glm::translate(plyM,glm::vec3(-9,-4.0,-14.5f));
  plyM = glm::rotate(plyM,angle2,glm::vec3(1.0,0.0,0.0));
  angle2 = 170*3.14/180;
  plyM = glm::rotate(plyM,angle2,glm::vec3(0.0,1.0,0.0));
  plyM = glm::scale(plyM,glm::vec3(10/(ply_size3[0]),3/(ply_size3[1]),3/(ply_size3[2]))); 
  SetMaterialColor(Kd,0.8,0.2,0.2,1);
  plyM = plyM ;
  mvp = projection *view*plyM;
  mv = view*plyM;
  draw_ply3(&mvp[0][0],&mv[0][0],&view[0][0],color,c0,c1,c2,m1,m2,m3,m4,m5);
  /* 
  */

  /////////////////////////////////////////////////////////
  //////////////////////////	Cube	//////////////////
  ///////////////////////////////////////////////////////
    st.push(modelM);
	modelM = st.top();

	
  	use_texture = 2;
 	glUniform1i(b1,use_texture);
    tex_loc = glGetUniformLocation(programObject,"cubeMap");
    glUniform1i(tex_loc,1);


	modelM = glm::translate(modelM,glm::vec3(-0.1,-6.9,-10.5));
	st.push(modelM);
	modelM = glm::scale(modelM,glm::vec3(3.7,3.7,8.5));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 	//stand
	modelM = st.top();
	modelM = glm::scale(modelM,glm::vec3(3.4,5.0,8.5));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 	//stand2
	modelM = st.top();
	modelM = glm::scale(modelM,glm::vec3(3.1,6.3,8.5));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 	//stand3
		modelM = st.top();
	modelM = glm::scale(modelM,glm::vec3(2.8,7.6,8.5));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 	//stand4
			modelM = st.top();
	modelM = glm::scale(modelM,glm::vec3(2.5,8.9,8.5));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 	//stand5
	st.pop();
	
	
	modelM = st.top();		//pops out stand data



	modelM = glm::translate(modelM,glm::vec3(-0.1,-6.9,-6.5));	//pedestral 1
	st.push(modelM);
	modelM = glm::scale(modelM,glm::vec3(8.18,8.18,2.7));
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 
	
	modelM = st.top();

  	modelM = glm::scale(modelM,glm::vec3(6.18,10.18,2.7));		//pedestral 2
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 


	modelM = st.top();
  	modelM = glm::scale(modelM,glm::vec3(4.18,12.18,2.7));		//pedestral 2
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
	mvp = projection *view*modelM;
	mv = view*modelM;
	draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4, m5); 
	 st.pop();
  ///////////////////////////////////////////////////////////////////////
  ////////////////////////	Wall				////////////////////////
  /////////////////////////////////////////////////////////////////////
   	modelM = st.top();


    use_texture = 3;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex3");
    glUniform1i(tex_loc,3);
	
	
	SetMaterialColor(Kd, 1.0, 1.0, 1.0, 1); 
  modelM = glm::translate(modelM,glm::vec3(-0.0, -0.0, -15.0));		//negative y
   modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 

 
 
  ////////////////////////////////////////////////////////////////////////////////////////////
      use_texture = 6;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex6");
    glUniform1i(tex_loc,6);
   modelM = st.top();
  modelM =  glm::translate(modelM,glm::vec3(0.0,0.0,+15.0));	//positive y
  modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 

  ////////////////////////////////////////////////////////////////////////////////////////////
  use_texture = 7;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex7");
    glUniform1i(tex_loc,7);
  modelM = st.top();
  modelM = glm::rotate(modelM, angle, glm::vec3(1.0, 0.0, 0.0));
  modelM =  glm::translate(modelM,glm::vec3(0.0,0.0,+15));	//positive z
 
  modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 


////////////////////////////////////////////////////////////////////////////////////////////
    use_texture = 4;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex4");
    glUniform1i(tex_loc,4);

  modelM = st.top();
  modelM = glm::rotate(modelM, angle, glm::vec3(1.0, 0.0, 0.0));
  modelM =  glm::translate(modelM,glm::vec3(0.0,0.0,-15));	//negative z
  modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 
 
  ////////////////////////////////////////////////////////////////////////////////////////////

        use_texture = 1;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex1");
    glUniform1i(tex_loc,0);

  modelM = st.top();
  modelM = glm::rotate(modelM, angle, glm::vec3(0.0, -1.0, 0.0));
  modelM =  glm::translate(modelM,glm::vec3(0.0,0.0,-15));	//negative x
  modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 

  ////////////////////////////////////////////////////////////////////////////////////////////

        use_texture = 5;
	glUniform1i(b1,use_texture);
	tex_loc = glGetUniformLocation(programObject,"Tex5");
    glUniform1i(tex_loc,5);
    modelM = st.top();
  modelM = glm::rotate(modelM, angle, glm::vec3(0.0, -1.0, 0.0));
  modelM =  glm::translate(modelM,glm::vec3(0.0,0.0,15));	//positive x
  modelM = glm::scale(modelM, glm::vec3(30.2, 30.2, 1.0));
  mvp = projection *view*modelM;
  mv = view*modelM;
  draw_square(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, c3,m1, m2, m3,m4); 
  modelM = st.top();
  

  
  st.pop();
  glDisableClientState(GL_VERTEX_ARRAY);
  glutSwapBuffers();

}
///////////////////////////////////////////////////////////
//    Find the geometry center of the input model 
////////////////////////////////////////////////////////////
void find_center (float& cx, float& cy, float& cz, 
		  float& minx, float& maxx, float&miny, 
		  float &maxy, float &minz, float & maxz)
{
  float x, y, z; 
  float min_x=9999, max_x=-9999, min_y=9999, max_y=-9999; 
  float min_z=9999, max_z=-9999; 

  x = y= z = 0; 
  for(int i = 0; i < pl_vertexcount; i++) {
    x += pl_vertices[i]->x; 
    y += pl_vertices[i]->y; 
    z += pl_vertices[i]->z; 
    if (min_x >pl_vertices[i]->x) min_x = pl_vertices[i]->x; 
    if (max_x <pl_vertices[i]->x) max_x = pl_vertices[i]->x; 
    if (min_y >pl_vertices[i]->y) min_y = pl_vertices[i]->y; 
    if (max_y <pl_vertices[i]->y) max_y = pl_vertices[i]->y; 
    if (min_z >pl_vertices[i]->z) min_z = pl_vertices[i]->z; 
    if (max_z <pl_vertices[i]->z) max_z = pl_vertices[i]->z; 
  }
  cx = x / (float) pl_vertexcount; 
  cy = y / (float) pl_vertexcount; 
  cz = z / (float) pl_vertexcount; 
  minx = min_x; maxx = max_x; 
  miny = min_y; maxy = max_y; 
  minz = min_z; maxz = max_z; 
}

void find_center2 (float& cx, float& cy, float& cz, 
		  float& minx, float& maxx, float&miny, 
		  float &maxy, float &minz, float & maxz)
{
  float x, y, z; 
  float min_x=9999, max_x=-9999, min_y=9999, max_y=-9999; 
  float min_z=9999, max_z=-9999; 

  x = y= z = 0; 
  for(int i = 0; i < pl_vertexcount2; i++) {
    x += pl_vertices2[i]->x; 
    y += pl_vertices2[i]->y; 
    z += pl_vertices2[i]->z; 
    if (min_x >pl_vertices2[i]->x) min_x = pl_vertices2[i]->x; 
    if (max_x <pl_vertices2[i]->x) max_x = pl_vertices2[i]->x; 
    if (min_y >pl_vertices2[i]->y) min_y = pl_vertices2[i]->y; 
    if (max_y <pl_vertices2[i]->y) max_y = pl_vertices2[i]->y; 
    if (min_z >pl_vertices2[i]->z) min_z = pl_vertices2[i]->z; 
    if (max_z <pl_vertices2[i]->z) max_z = pl_vertices2[i]->z; 
  }
  cx = x / (float) pl_vertexcount2; 
  cy = y / (float) pl_vertexcount2; 
  cz = z / (float) pl_vertexcount2; 
  minx = min_x; maxx = max_x; 
  miny = min_y; maxy = max_y; 
  minz = min_z; maxz = max_z; 
}

float angle1=0.0;	float angle2=0.0;	float angle3=0.0;	float angle4=0.0;

void mykey(unsigned char key, int x, int y)
{
        float d_angle = 10; 
	if (key == 'q') exit(1); 
	if (key == 'R') 
	  modelM = glm::rotate(modelM, d_angle, glm::vec3(0.0f, 0.0f, 1.0f)); 
	//  modelM =  rotatez(modelM, d_angle);
	if (key == 'd') 
	  modelM = glm::translate(modelM, glm::vec3(0.1f, 0.0f, 0.0f)); 
	//  modelM =  translate(modelM, .1,0,0);
	if (key == 'a') 
	  modelM = glm::translate(modelM, glm::vec3(-0.1f, 0.0f, 0.0f)); 
	//modelM =  translate(modelM, -.1,0,0);
	if (key == 'w') 
	  modelM = glm::translate(modelM, glm::vec3(0.0f, 0.1f, 0.0f)); 
	//modelM =  translate(modelM, 0,.1,0);
	if (key == 's') 
	  modelM = glm::translate(modelM, glm::vec3(0.0f, -0.1f, 0.0f)); 
	//modelM =  translate(modelM, 0,-.1,0);
	if (key == 'c') {
	  modelM =  glm::mat4(1.0f);
	  angle1 = angle2 = angle3 = angle4 = 0; 
	}
	
	if (key == '+') 
          mat_shine[0] += 1; 	
		if (key == '-') 
          mat_shine[0] -= 1; 	
        if (key == '3') 
		{light_pos[0] += 0.1; printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
        if (key == '4') 
	  {light_pos[0] -= 0.1;  printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
        if (key == '5') 
	  {light_pos[1] += 0.1;  printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
        if (key == '6') 
	  {light_pos[1] -= 0.1;  printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
        if (key == '7') 
	  {light_pos[2] += 0.1;  printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
        if (key == '8') 
	  {light_pos[2] -= 0.1; printf("light_pos[0]: %f, light_pos[1]: %f, light_pos[2]: %f\n",light_pos[0],light_pos[1],light_pos[2]);}
		if(key == ' ')
			light_ambient[0]+=0.5;//light_ambient[1]=light_ambient[2]=light_ambient[0]+0.5;
		
	
	glutPostRedisplay(); 
}


void mymouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    press_x_cube = x; press_y_cube = y; press_x_cyl = x; press_y_cyl = y; 
    if (button == GLUT_LEFT_BUTTON)
	{xform_mode_cube = XFORM_ROTATE; xform_mode_cyl = XFORM_ROTATE; }
	 else if (button == GLUT_RIGHT_BUTTON) 
	 {xform_mode_cube = XFORM_SCALE; xform_mode_cyl = XFORM_SCALE; }
  }
  else if (state == GLUT_UP) {
	  { xform_mode_cube = XFORM_NONE;  xform_mode_cyl = XFORM_NONE; }
  }



}



void mymotion(int x, int y)
{
 

 if (xform_mode_cyl==XFORM_ROTATE) {
      x_angle_rot += (x - press_x_cyl)/20.0; 
      if (x_angle_rot > 180) x_angle_rot -= 360; 
      else if (x_angle_rot <-180) x_angle_rot += 360; 
      press_x_cyl = x; 
	   rotate_x = true;
	   rotate_y = true;
      y_angle_rot += (y - press_y_cyl)/40.0; 
      if (y_angle_rot > 180) y_angle_rot -= 360; 
      else if (y_angle_rot <-180) y_angle_rot += 360; 
      press_y_cyl = y; 
     }
	else if (xform_mode_cyl == XFORM_SCALE){
	  samescale = false;
	  float old_size = scale_size_cyl;
      scale_size_cyl *= (1+ (y - press_y_cyl)/60.0); 
      if (scale_size_cyl <0) scale_size_cyl = old_size; 
      press_y_cyl = y; 
    }


    glutPostRedisplay(); 

}

int main(int argc, char** argv) 
{

  PlyFile* input;
  PlyFile* input2;
  PlyFile* input3;
  // get the ply structure and open the file
  input = read_ply(fopen("happy2.ply","r"));
  input2 = read_ply(fopen("hind2.ply","r"));
  input3 = read_ply(fopen("dragon2.ply","r"));
  // read in the data
  store_ply(input, 
	    &pl_vertices, &pl_faces, 
	    &pl_vertexcount, &pl_facecount,
	    &pl_vertexnormals, &pl_facenormals);
  store_ply(input2, 
	    &pl_vertices2, &pl_faces2, 
	    &pl_vertexcount2, &pl_facecount2,
	    &pl_vertexnormals2, &pl_facenormals2);
  store_ply(input3, 
	    &pl_vertices3, &pl_faces3, 
	    &pl_vertexcount3, &pl_facecount3,
	    &pl_vertexnormals3, &pl_facenormals3);
  // close the file*/
  close_ply(input);
  close_ply(input2);
  close_ply(input3);
  glutInit(&argc,argv);
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH);
  glutInitWindowSize(600,600);
  glutCreateWindow("VBO simple ply");
  glutDisplayFunc(display);
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 

  GLenum err = glewInit(); 

  InitPLY_VBO(ply_min,ply_max);
  InitPLY_VBO2(ply_min2,ply_max2);
  InitPLY_VBO3(ply_min3,ply_max3);
  Initsquare();
   InitVBO_square(); 
   InitCube_VBO();
 // Initsquare_VBO();
  //InitSquare();
  //InitSquare_VBO();

   int nslices, nstacks; 
  nslices = 60; 
  nstacks = 20; 
  //model = glm::translate(model,glm::vec3(0,1.0,0));
//  InitCylinder(nslices, nstacks, 1.0, 0.0, 1.0); 

   glewInit(); 
  InitVBO_cyl(nslices,nstacks);
  InitSphere(nslices, nstacks, 1.0, 0.0, 1.0);
  InitVBO_sphere(nslices, nstacks); 


  Init_texture(1);
  init_cubemap();
  programObject = SetupGLSL("plyTextureEnv");
  find_center(cx, cy, cz, x_min, x_max, 
	      y_min, y_max, z_min, z_max); 
  printf("geometry center = %f %f %f \n", cx, cy, cz); 
  printf("geometry bound = x: %f %f y: %f %f z: %f %f\n", 
	 x_min, x_max, y_min, y_max, z_min, z_max); 

  find_center2(cx, cy, cz, x_min, x_max, 
	      y_min, y_max, z_min, z_max); 
  printf("geometry center2 = %f %f %f \n", cx, cy, cz); 
  printf("geometry bound2 = x: %f %f y: %f %f z: %f %f\n", 
	 x_min, x_max, y_min, y_max, z_min, z_max); 
   glutMainLoop();
  //getchar();
}

