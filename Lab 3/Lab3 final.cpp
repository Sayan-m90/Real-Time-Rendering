#include <stdio.h>
#include <cstdlib>
#include <stack> 

#ifdef __APPLE__
#include <GLUT/glut.h> 
#include <OpenGL/gl.h> 
#else
#include <GL/glew.h> 
#include <GL/glut.h> 
#include <GL/gl.h> 
#endif
#include<Windows.h>
#include<glm/glm.hpp>
#include<glm/gtx/transform.hpp>
//#include"cylinder.cpp"


//////////////////////////////////////////////////////////////////////////////////
//
//  Just hard code the geometry 
// 
//////////////////////////////////////////////////////////////////////////////
typedef struct 
{
  float sph_location[4]; 
  float sph_normal[4]; 
  float sph_color[4]; 
} Vertexsphere; 

int sphindices; 
Vertexsphere *sphverts;   // cylinder vertices
GLuint *sp_indices; 

GLuint sphere_vboHandle[1];   // a VBO that contains interleaved positions and sph_colors 
GLuint sphere_indexVBO; 

///////////////////////////////////////////////////////////////

typedef struct 
{
  float cube_location[4]; 
  float cube_color[4]; 
} Vertexcube; 

Vertexcube verts[8];       // cube vertices 

GLubyte tindices[6*6];   // cube face indices: 6 faces, 6 indices each face 

GLuint cube_vboHandle[1];   // a VBO that contains interleaved positions and cube_colors 
GLuint cube_indexVBO; 



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

bool samescale = true;
float M_PI = 3.14127;
glm::mat4 model = glm::mat4(1.0f); 
glm::mat4 model1 = glm::mat4(1.0f); 
glm::mat4 model2 = glm::mat4(1.0f); 
glm::mat4 model3 = glm::mat4(1.0f); 
/////////////////////////////////
// glut mouse control 
// 
int xform_mode_cube = 0; 
#define XFORM_NONE    0 
#define XFORM_ROTATE  1
#define XFORM_SCALE 2 

bool hello_flag = false;
bool rotate_x = false;
bool rotate_y = false;
void mymouse(int, int, int, int);
void mymotion(int, int);
void mykey(unsigned char, int, int);
void display_cyl();
int press_x_cube, press_y_cube; 
int release_x_cube, release_y_cube; 
float z_coord = 0.0f;
bool forwardmove = false;
float x_angle_cube = 0.0; 
float y_angle_cube = 0.0; 
float scale_size_cube = 1; 

std::stack<glm::mat4> mv;

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
bool hand_move_flag = true;
bool hand_move_forward = true;
bool r_arm_up = true;
bool l_arm_up = true;
bool r_leg_front = true;
bool r_lowerleg_front = false;
bool l_lowerleg_front = false;
float legangle1 = 0.0f;
float lowerlegangle1 = 0.0f;
float legangle2 = 0.0f;
float lowerlegangle2 = 0.0f;

float head = 0.0f;
bool head_flag = false;
//float headl = 0.0f;

int xform_mode_cyl = 0;
int press_x_cyl, press_y_cyl; 
int release_x_cyl, release_y_cyl; 
float x_angle_rot = 0.0;		//rotate whole scene by mouse
float x_angle_body = 0.0;		// rotate body by keyboard
bool  x_angle_body_flag = true;
float y_angle_rot = 0.0; 
float scale_size_cyl = 1; 
void InitVBO_cyl(int,int); 
void InitCylinder(int , int, float, float, float) ;

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


void InitSphere(int nslices, int nstacks, float r, float g, float b) 
{

  float istep = 2 * 3.14159/(float)(nslices); 
  float jstep = 2 * 3.14159 /(float)(nstacks); 

  int nvertices = nslices * nstacks; 
  sphverts = new Vertexsphere[nvertices]; 

  for (int j =0; j<nstacks; j++)
    for (int i=0; i<nslices; i++) {
      int idx = j*nslices + i; // mesh[j][i] 

      float s = i * istep; 
      float t = j * jstep; 
      float x = 4*cos(t)*sin(s);//3*cos(s)+cos(t)*cos(s); 
      float y = 4*sin(t)*sin(s);//3*cos(s)+cos(t)*sin(s); 
	  float z = 4*sin(s);//sin(t); 
	 
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

      sphverts[idx].sph_normal[0] = sph_normal[0]; 
      sphverts[idx].sph_normal[1] = sph_normal[1]; 
      sphverts[idx].sph_normal[2] = sph_normal[2]; 

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

  glGenBuffers(1, sphere_vboHandle);   // create an interleaved VBO object
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

//////////////////////////////////////////////////////////////////////////////////
void InitGeometry()
{
  // negative z faces vertices 
  verts[0].cube_location[0] = -0.5; verts[0].cube_location[1] = -0.5;  verts[0].cube_location[2] = -0.5;  verts[0].cube_location[3] = 1; 
  verts[1].cube_location[0] = -0.5; verts[1].cube_location[1] = 0.5;  verts[1].cube_location[2] = -0.5;  verts[1].cube_location[3] = 1; 
  verts[2].cube_location[0] = 0.5; verts[2].cube_location[1] = 0.5;  verts[2].cube_location[2] = -0.5;  verts[2].cube_location[3] = 1; 
  verts[3].cube_location[0] = 0.5; verts[3].cube_location[1] = -0.5;  verts[3].cube_location[2] = -0.5;  verts[3].cube_location[3] = 1; 

  // positive  z faces vertices 
  verts[4].cube_location[0] = -0.5; verts[4].cube_location[1] = -0.5;  verts[4].cube_location[2] = 0.5;  verts[4].cube_location[3] = 1; 
  verts[5].cube_location[0] = -0.5; verts[5].cube_location[1] = 0.5;  verts[5].cube_location[2] = 0.5;  verts[5].cube_location[3] = 1; 
  verts[6].cube_location[0] = 0.5; verts[6].cube_location[1] = 0.5;  verts[6].cube_location[2] = 0.5;  verts[6].cube_location[3] = 1; 
  verts[7].cube_location[0] = 0.5; verts[7].cube_location[1] = -0.5;  verts[7].cube_location[2] = 0.5;  verts[7].cube_location[3] = 1; 

  verts[0].cube_color[0] = 0; verts[0].cube_color[1] = 0;  verts[0].cube_color[2] = 1;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0; verts[1].cube_color[1] = 0;  verts[1].cube_color[2] = 1;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 1; verts[2].cube_color[1] = 1;  verts[2].cube_color[2] = 0;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 1; verts[3].cube_color[1] = 1;  verts[3].cube_color[2] = 0;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 1; verts[4].cube_color[1] = 1;  verts[4].cube_color[2] = 0;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 1; verts[5].cube_color[1] = 1;  verts[5].cube_color[2] = 0;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0; verts[6].cube_color[1] = 0;  verts[6].cube_color[2] = 1;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0; verts[7].cube_color[1] = 0;  verts[7].cube_color[2] = 1;  verts[7].cube_color[3] = 1; 
 // changeCubeColor();
  // create cube vertex indices. 

  // negative z face 
  tindices[0] = 0;   tindices[1] = 1;   tindices[2] = 2; 
  tindices[3] = 2;   tindices[4] = 3;   tindices[5] = 0; 
  // positive z face 
  tindices[6] = 4;   tindices[7] = 5;   tindices[8] = 6; 
  tindices[9] = 6;   tindices[10] = 7;   tindices[11] =4;
  // negative y face 
  tindices[12] = 4;   tindices[13] = 0;   tindices[14] = 3; 
  tindices[15] = 3;   tindices[16] = 7;   tindices[17] = 4;  
  // positive y face 
  tindices[18] = 5;   tindices[19] = 1;   tindices[20] = 2; 
  tindices[21] = 2;   tindices[22] = 6;   tindices[23] = 5;  
  // negative x face 
  tindices[24] = 4;   tindices[25] = 5;   tindices[26] = 1; 
  tindices[27] = 1;   tindices[28] = 0;   tindices[29] = 4;  
  // positive x face 
  tindices[30] = 7;   tindices[31] = 6;   tindices[32] = 2; 
  tindices[33] = 2;   tindices[34] = 3;   tindices[35] = 7;  
}


//////////////////////////////////////////////////////////////////////////////////
//
// create VBO objects and send the triangle vertices/cube_colors to the graphics card
// 
void InitVBO_cube() 
{
  glGenBuffers(1, cube_vboHandle);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexcube)*8, verts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &cube_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLubyte)*36, tindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

} 


void InitVBO_cyl(int nslices, int nstacks) 
{
  int nvertices = nslices * nstacks; 
  nindices = (nstacks-1)*2*(nslices+1); 

  glGenBuffers(1, cyl_vboHandle);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cyl_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexcyl)*nvertices, cyverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &cyl_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cyl_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*nindices, cindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

   //by now, we moved the position and cyl_color data over to the graphics card. There will be no redundant data copy at drawing time 
} 

void initChangeColor(int choice)
{
  InitGeometry();

    if(choice == 0)
  {
  verts[0].cube_color[0] = 0.99; verts[0].cube_color[1] = 0.33;  verts[0].cube_color[2] = 0.1;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.99; verts[1].cube_color[1] = 0.33;  verts[1].cube_color[2] = 0.1;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.99; verts[2].cube_color[1] = 0.78;  verts[2].cube_color[2] = 0.49;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.99; verts[3].cube_color[1] = 0.78;  verts[3].cube_color[2] = 0.49;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.99; verts[4].cube_color[1] = 0.33;  verts[4].cube_color[2] = 0.1;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.99; verts[5].cube_color[1] = 0.33;  verts[5].cube_color[2] = 0.1;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.99; verts[6].cube_color[1] = 0.78;  verts[6].cube_color[2] = 0.49;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.99; verts[7].cube_color[1] = 0.78;  verts[7].cube_color[2] = 0.49;  verts[7].cube_color[3] = 1; 
  }
  if(choice == 1)
  {
  verts[0].cube_color[0] = 0.7; verts[0].cube_color[1] = 0.58;  verts[0].cube_color[2] = 0.35;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.7; verts[1].cube_color[1] = 0.58;  verts[1].cube_color[2] = 0.35;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.85; verts[2].cube_color[1] = 0.74;  verts[2].cube_color[2] = 0.66;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.85; verts[3].cube_color[1] = 0.74;  verts[3].cube_color[2] = 0.66;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.7; verts[4].cube_color[1] = 0.58;  verts[4].cube_color[2] = 0.35;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.7; verts[5].cube_color[1] = 0.58;  verts[5].cube_color[2] = 0.35;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.85; verts[6].cube_color[1] = 0.74;  verts[6].cube_color[2] = 0.66;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.85; verts[7].cube_color[1] = 0.74;  verts[7].cube_color[2] = 0.66;  verts[7].cube_color[3] = 1; 
  }

   if(choice == 8)
  {
  verts[0].cube_color[0] = 0.99; verts[0].cube_color[1] = 0.33;  verts[0].cube_color[2] = 0.1;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.99; verts[1].cube_color[1] = 0.33;  verts[1].cube_color[2] = 0.1;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.99; verts[2].cube_color[1] = 0.98;  verts[2].cube_color[2] = 0.49;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.99; verts[3].cube_color[1] = 0.98;  verts[3].cube_color[2] = 0.49;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.99; verts[4].cube_color[1] = 0.33;  verts[4].cube_color[2] = 0.1;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.99; verts[5].cube_color[1] = 0.33;  verts[5].cube_color[2] = 0.1;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.99; verts[6].cube_color[1] = 0.98;  verts[6].cube_color[2] = 0.49;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.99; verts[7].cube_color[1] = 0.98;  verts[7].cube_color[2] = 0.49;  verts[7].cube_color[3] = 1; 
   }

  if(choice == 2)
  {
  verts[0].cube_color[0] = 0.39; verts[0].cube_color[1] = 0.21;  verts[0].cube_color[2] = 0;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.39; verts[1].cube_color[1] = 0.21;  verts[1].cube_color[2] = 0;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.68; verts[2].cube_color[1] = 0.37;  verts[2].cube_color[2] = 0;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.68; verts[3].cube_color[1] = 0.37;  verts[3].cube_color[2] = 0;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.68; verts[4].cube_color[1] = 0.37;  verts[4].cube_color[2] = 0;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.68; verts[5].cube_color[1] = 0.37;  verts[5].cube_color[2] = 0;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.39; verts[6].cube_color[1] = 0.21;  verts[6].cube_color[2] = 0;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.39; verts[7].cube_color[1] = 0.21;  verts[7].cube_color[2] = 0;  verts[7].cube_color[3] = 1; 
  }
  if(choice == 3)
  {
  verts[0].cube_color[0] = 0.13; verts[0].cube_color[1] = 0.10;  verts[0].cube_color[2] = 0;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.13; verts[1].cube_color[1] = 0.10;  verts[1].cube_color[2] = 0;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.13; verts[2].cube_color[1] = 0.10;  verts[2].cube_color[2] = 0;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.13; verts[3].cube_color[1] = 0.10;  verts[3].cube_color[2] = 0;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.13; verts[4].cube_color[1] = 0.10;  verts[4].cube_color[2] = 0;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.53; verts[5].cube_color[1] = 0.30;  verts[5].cube_color[2] = 0;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.13; verts[6].cube_color[1] = 0.10;  verts[6].cube_color[2] = 0;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.13; verts[7].cube_color[1] = 0.10;  verts[7].cube_color[2] = 0;  verts[7].cube_color[3] = 1; 
  }

  if(choice == 4)
  {
  verts[0].cube_color[0] = 0.95; verts[0].cube_color[1] = 0.72;  verts[0].cube_color[2] = 0.46;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.97; verts[1].cube_color[1] = 0.81;  verts[1].cube_color[2] = 0.64;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.97; verts[2].cube_color[1] = 0.81;  verts[2].cube_color[2] = 0.64;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.95; verts[3].cube_color[1] = 0.72;  verts[3].cube_color[2] = 0.46;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.95; verts[4].cube_color[1] = 0.72;  verts[4].cube_color[2] = 0.46;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.97; verts[5].cube_color[1] = 0.81;  verts[5].cube_color[2] = 0.64;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.97; verts[6].cube_color[1] = 0.81;  verts[6].cube_color[2] = 0.64;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.95; verts[7].cube_color[1] = 0.72;  verts[7].cube_color[2] = 0.46;  verts[7].cube_color[3] = 1; 
  }
  if(choice == 5)
  {
	 
  verts[0].cube_color[0] = 0.97; verts[0].cube_color[1] = 0.81;  verts[0].cube_color[2] = 0.64;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.95; verts[1].cube_color[1] = 0.72;  verts[1].cube_color[2] = 0.46;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.95; verts[2].cube_color[1] = 0.72;  verts[2].cube_color[2] = 0.46;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.97; verts[3].cube_color[1] = 0.81;  verts[3].cube_color[2] = 0.64;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.97; verts[4].cube_color[1] = 0.81;  verts[4].cube_color[2] = 0.64;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.95; verts[5].cube_color[1] = 0.72;  verts[5].cube_color[2] = 0.46;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.95; verts[6].cube_color[1] = 0.72;  verts[6].cube_color[2] = 0.46;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.97; verts[7].cube_color[1] = 0.81;  verts[7].cube_color[2] = 0.64;  verts[7].cube_color[3] = 1; 

  }
   if(choice == 6)
  {
  verts[0].cube_color[0] = 0.95; verts[0].cube_color[1] = 0.72;  verts[0].cube_color[2] = 0.46;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0.95; verts[1].cube_color[1] = 0.72;  verts[1].cube_color[2] = 0.46;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0.95; verts[2].cube_color[1] = 0.72;  verts[2].cube_color[2] = 0.46;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0.95; verts[3].cube_color[1] = 0.72;  verts[3].cube_color[2] = 0.46;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0.95; verts[4].cube_color[1] = 0.72;  verts[4].cube_color[2] = 0.46;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0.95; verts[5].cube_color[1] = 0.72;  verts[5].cube_color[2] = 0.46;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0.95; verts[6].cube_color[1] = 0.72;  verts[6].cube_color[2] = 0.46;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0.95; verts[7].cube_color[1] = 0.72;  verts[7].cube_color[2] = 0.46;  verts[7].cube_color[3] = 1; 
  }

     if(choice == 7)
  {
  verts[0].cube_color[0] = 0; verts[0].cube_color[1] = 0.55;  verts[0].cube_color[2] = 0.30;  verts[0].cube_color[3] = 1;   
  verts[1].cube_color[0] = 0; verts[1].cube_color[1] = 0.55;  verts[1].cube_color[2] = 0.30;  verts[1].cube_color[3] = 1; 
  verts[2].cube_color[0] = 0; verts[2].cube_color[1] = 0.55;  verts[2].cube_color[2] = 0.30;  verts[2].cube_color[3] = 1; 
  verts[3].cube_color[0] = 0; verts[3].cube_color[1] = 0.55;  verts[3].cube_color[2] = 0.30;  verts[3].cube_color[3] = 1; 

  verts[4].cube_color[0] = 0; verts[4].cube_color[1] = 0.72;  verts[4].cube_color[2] = 0.96;  verts[4].cube_color[3] = 1; 
  verts[5].cube_color[0] = 0; verts[5].cube_color[1] = 0.72;  verts[5].cube_color[2] = 0.96;  verts[5].cube_color[3] = 1; 
  verts[6].cube_color[0] = 0; verts[6].cube_color[1] = 0.72;  verts[6].cube_color[2] = 0.96;  verts[6].cube_color[3] = 1; 
  verts[7].cube_color[0] = 0; verts[7].cube_color[1] = 0.72;  verts[7].cube_color[2] = 0.96;  verts[7].cube_color[3] = 1; 
  }

    InitVBO_cube();
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 

  glColorPointer(4, GL_FLOAT, sizeof(Vertexcube), (char*) NULL+ 16); 
  glVertexPointer(4,GL_FLOAT, sizeof(Vertexcube), (char*) NULL+ 0); 
}

bool first = true;
void display_cube() 
{ 

  glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects 
  
  //changeCubeColor();
 
  glm::mat4 projection = glm::perspective(60.0f,1.0f,.1f,100.0f); 
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixf(&projection[0][0]); 

  glm::mat4 view = glm::lookAt(glm::vec3(0.0, 0.0, 5.0), 
			       glm::vec3(0.0, 0.0, 0.0), 
			       glm::vec3(0.0, 1.0, 0.0)); 
 if(first == true)
  {model = glm::translate(model, glm::vec3(0.0, 1.0, 0.0));
   first = false;
  }

 
  if(rotate_x == true) 
  {
	  model = glm::rotate(model, x_angle_rot, glm::vec3(1.0f, 0.0f, 0.0f));
	  rotate_x = false;
	  model2 = glm::rotate(model2, x_angle_rot, glm::vec3(1.0f, 0.0f, 0.0f));
  }

  if(rotate_y == true) 
	{model = glm::rotate(model, y_angle_rot, glm::vec3(0.0f, 1.0f, 0.0f)); 
	 model2 = glm::rotate(model2, y_angle_rot, glm::vec3(0.0f, 1.0f, 0.0f)); 
	 rotate_y = false;
    }
   
   if(samescale==false) 
	  { model = glm::scale(model, glm::vec3(scale_size_cyl, scale_size_cyl, scale_size_cyl)); 
		model2 = glm::scale(model, glm::vec3(scale_size_cyl, scale_size_cyl, scale_size_cyl)); 
		samescale = true;
		}
   if(forwardmove == true) 
   {model = glm::translate(model, glm::vec3(0.0, 0.0, z_coord)); forwardmove = false;}
   if(x_angle_body_flag == true)
   {
	   model = glm::rotate(model, x_angle_body , glm::vec3(0.0f, 1.0f, 0.0f)); 
	   x_angle_body_flag = false;
   }

   mv.push(model);


  if (draw_line) glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
	  else glPolygonMode(GL_FRONT_AND_BACK,GL_FILL); 
	
	  //////////////////////////////////////
	  ////////////////////   Sphere
	  //////////////////////////////////////
  glBindBuffer(GL_ARRAY_BUFFER, sphere_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, sphere_indexVBO); 

  glEnableClientState(GL_VERTEX_ARRAY); // enable the Vertexsphere array on the client side
  glEnableClientState(GL_NORMAL_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY); // enable the sph_color array on the client side

  glVertexPointer(4,GL_FLOAT, sizeof(Vertexsphere),(char*) NULL+0); 
  glNormalPointer(GL_FLOAT, sizeof(Vertexsphere),(char*) NULL+16); 
  glColorPointer(4,GL_FLOAT,  sizeof(Vertexsphere),(char*) NULL+32);

  // Now we are ready to draw, using the triangle indices in the buffer 

  glMatrixMode(GL_PROJECTION); 
   //model = glm::translate(model, glm::vec3(0.0, 0.0, 5.0f));   
  model = glm::translate(model, glm::vec3(0.0, -3.0, -0.1f));   
  model = glm::scale(model, glm::vec3(0.2,0.2,0.2)); 

  glm::mat4 temp = glm::mat4(1.0f); 
  temp = model;
   if(head_flag == true) 
   {temp = glm::rotate(model, head, glm::vec3(0.0,1.0,0.0));	}//	head_flag = false;	}

   glm::mat4 modelview = view*model;
   glm::mat4 modelview1 = view * temp; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview1[0][0]); 
  

  glDrawElements(GL_TRIANGLES, sphindices, GL_UNSIGNED_INT, (char*) NULL+0); 

  glDisableClientState(GL_VERTEX_ARRAY); // enable the Vertexsphere array on the client side
  glDisableClientState(GL_COLOR_ARRAY); // enable the sph_color array on the client side
  
  ////////////////////////////////////////////////////// 
  ////////////////////    Cylinder    /////////////////
  ////////////////////////////////////////////////////
  
  
  glBindBuffer(GL_ARRAY_BUFFER, cyl_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cyl_indexVBO); 

  glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
  glEnableClientState(GL_NORMAL_ARRAY); 
  glEnableClientState(GL_COLOR_ARRAY); // enable the cyl_color array on the client side

  glVertexPointer(4,GL_FLOAT, sizeof(Vertexcyl),(char*) NULL+0); 
  glNormalPointer(GL_FLOAT, sizeof(Vertexcyl),(char*) NULL+16); 
  glColorPointer(4,GL_FLOAT,  sizeof(Vertexcyl),(char*) NULL+32);

  

   model=mv.top();
  
   model = glm::rotate(model, 1.57f, glm::vec3(1.0f, 0.0f, 0.0f)); 
   model = glm::scale(model, glm::vec3(1.0f, 1.0f, 3.0f)); 
   model = glm::translate(model,glm::vec3(0.0f,0.0f,-0.8f));
   
   

   model1 = model;
   model1 = glm::scale(model1, glm::vec3(1.1f, 1.0f, 0.4f));

  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0); //cylinder
  modelview = view * model1; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLE_STRIP, nindices, GL_UNSIGNED_INT, (char*) NULL+0); //cylinder2
  glDisableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
  glDisableClientState(GL_COLOR_ARRAY); // enable the cyl_color array on the client side
  glDisableClientState(GL_NORMAL_ARRAY);
  glDisable(GL_DEPTH_TEST);
  
  //////////////////////////////////////////
  ///////////// Cubes          /////////////
  /////////////////////////////////////////
  
  initChangeColor(0);
  glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects 
  

  glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array on the client side
  glEnableClientState(GL_COLOR_ARRAY); // enable the cube_color array on the client side

 
   view = glm::lookAt(glm::vec3(0.0, 0.0, 5.0), 
			       glm::vec3(0.0, 0.0, 0.0), 
			       glm::vec3(0.0, 1.0, 0.0)); 
  
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 

  glColorPointer(4, GL_FLOAT, sizeof(Vertexcube), (char*) NULL+ 16); 
  glVertexPointer(4,GL_FLOAT, sizeof(Vertexcube), (char*) NULL+ 0); 

  
 
  // now do a push matrix 
  //mv.push(model); // just like glPushMatrix(); 
  glMatrixMode(GL_PROJECTION);
  glLoadMatrixf(&projection[0][0]); 
 
  initChangeColor(1);
  ////////////////////////////////////////
  //////          base
  ///////////////////////////////////////
  //model=mv.top();
  model3 = glm::translate(model2, glm::vec3(0.0,8.0 , 0.0)); 
  model3 = glm::scale(model3, glm::vec3(50.0f, 0.20f, 10.0f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  
  model=mv.top();

  initChangeColor(8);
  ////////////////////////////////////////
  //////          base2
  ///////////////////////////////////////
  //model=mv.top();
  model3 = glm::translate(model2, glm::vec3(-1.7,8.0 , 0.0)); 
  model3 = glm::scale(model3, glm::vec3(10.0f, 0.20f, 10.0f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  
  model=mv.top();

  initChangeColor(0);
  ////////////////////////////////////////
  //////          back
  ///////////////////////////////////////
  //model3 = glm::translate(model2, glm::vec3(0.0,0.0 , -10.0)); 
  model3 = glm::scale(model2, glm::vec3(50.0f, 30.0f, 0.2f)); 

  // construct the modelview and modelview projection matrices 
  modelview = view * model3; 
 glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  
  model=mv.top();
 

  initChangeColor(7);
  ////////////////////////////////////////
  //////         boxes
  ///////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(7.5, 5.2 , 3.5)); 
  model3 = glm::scale(model3, glm::vec3(2.0f, 2.0f, 0.1f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
   
  model=mv.top();
  initChangeColor(2);
  ////////////////////////////////////////
  //////          table
  ///////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-8, 3.2 , 2.5)); 
  model3 = glm::scale(model3, glm::vec3(4.8f, 0.4f, 1.8f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  
  initChangeColor(3);
  //////////////////////////////////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-7,5.0 , 2.3)); 
  model3 = glm::scale(model3, glm::vec3(1.0f, 2.3f, 1.5f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);

  ////////////////////////////////////////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-11, 5.0 , 2.3)); 
  model3 = glm::scale(model3, glm::vec3(1.0f, 2.3f, 1.5f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);


  //////////////////////////////////////////////////////////////
  ///////////////////		chair		///////////////////////
  ////////////////////////////////////////////////////////////
  initChangeColor(2);
  model3 = glm::translate(model2, glm::vec3(8.0, 2.0, 2.9)); 
  model3 = glm::scale(model3, glm::vec3(1.5f, 4.0f, 0.05f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);


  /////////////////////////////////////////////////
  //initChangeColor(2);
  model3 = glm::translate(model2, glm::vec3(8.0, 3.0, 3.26)); 
  model3 = glm::scale(model3, glm::vec3(1.5f, 1.6f, 0.05f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);

    /////////////////////////////////////////////////
  initChangeColor(3);
  model3 = glm::translate(model2, glm::vec3(8.0, 2.0, 3.05)); 
  model3 = glm::scale(model3, glm::vec3(1.8f, 0.30f, 0.5f)); 

  modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);


  model=mv.top();
  initChangeColor(4);
  ////////////////////////////////////////
  //////          Head
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0,-2.3 , 0.1)); 
   model = glm::scale(model, glm::vec3(1.0, 1.0 , 1.2));
   temp =model;
  if(head_flag == true) 
  {temp = glm::rotate(model, head, glm::vec3(0.0,1.0,0.0));}//		head_flag = false;	}

  modelview1 = view * temp; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview1[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  
  model=mv.top();
 

   ////////////////////////////////////////
  //////          neck
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0,-1.8 , 0.0)); 
  model = glm::scale(model,glm::vec3(0.4,0.4,0.4));

  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  

   model=mv.top();
 
    initChangeColor(3);
   ////////////////////////////////////////
  //////          eyes
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.25, -2.4 , 0.8)); 
  model = glm::scale(model,glm::vec3(0.20,0.25,0.10));
  temp = model;
   if(head_flag == true) 
 {
	 temp = glm::translate(temp, glm::vec3(-0.25, 2.4 , -0.8)); 
	 temp = glm::rotate(temp, head, glm::vec3(0.0,1.0,0.0));
     temp = glm::translate(temp, glm::vec3(0.25, -2.4 , 0.8)); }
  modelview1 = view * temp; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview1[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
	model=mv.top();

 
    initChangeColor(3);
   ////////////////////////////////////////
  //////          eyes2
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(-0.25, -2.4 , 0.8)); 
  model = glm::scale(model,glm::vec3(0.20,0.25,0.10));
  temp = model;
   if(head_flag == true) 
 {
	 temp = glm::translate(temp, glm::vec3(0.25, 2.4 , -0.8));
	 temp = glm::rotate(temp, head, glm::vec3(0.0,1.0,0.0));		
	 temp = glm::translate(temp, glm::vec3(-0.25, -2.4 , 0.8));}
  modelview1 = view * temp; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview1[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
	model=mv.top();
 
	 initChangeColor(3);
	 	model=mv.top();

 

   ////////////////////////////////////////
  //////          mouth
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0, -2.0 , 0.8)); 
  model = glm::scale(model,glm::vec3(0.20,0.25,0.10));
  temp = model;
  if(head_flag == true) 
 {
	  //temp = glm::scale(temp,glm::vec3(1,1,1));
	  temp = glm::translate(temp, glm::vec3(0.0, 2.0 , -0.8)); 
	  temp = glm::rotate(temp, head, glm::vec3(0.0,1.0,0.0));		
	  temp = glm::translate(temp, glm::vec3(0.0, -2.0 , 0.8));
	  //temp = glm::scale(temp,glm::vec3(0.20,0.25,0.10));
	  head_flag = false;	}
	 
  modelview1 = view * temp; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview1[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
	model=mv.top();
 
	 initChangeColor(4);
  ////////////////////////////////////////
  //////          shoulder
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0,-1.2 , 0.0)); 
  model = glm::scale(model,glm::vec3(4.4,0.8,0.4));

  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);


  model=mv.top();
  model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  mv.push(model);
   model=mv.top();

  
  ////////////////////////////////////////
  //////          2nd box left hand
  ///////////////////////////////////////
  
   model = glm::translate(model, glm::vec3(2.1, 0.2, 0.0)); 
   model = glm::rotate(model, l_upper_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 
  
    if(hand_move_flag == true)
  {		
	model = glm::rotate(model, -hand_move , glm::vec3(1.0f, 0.0f, 0.0f)); 
  }
  modelview = view * model; 

  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);


  model = mv.top();  
  mv.pop();
  mv.push(model);

   initChangeColor(5);
  ////////////////////////////////////////
  //////          2nd box left hand lower arm
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(l_lower_arm_x, l_lower_arm_y, 0.0)); 
 // model = glm::rotate(model, 2.57f , glm::vec3(0.0f, 1.0f, 2.0f)); 
 model = glm::rotate(model, l_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 
 if(hand_move_flag == true)
  {		 
   model = glm::rotate(model, -l_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f));
   model = glm::translate(model, glm::vec3(-l_lower_arm_x, -l_lower_arm_y, 0.0)); 
   
   model = glm::rotate(model, -hand_move, glm::vec3(1.0f, 0.0f, 0.0f));
  

   
   model = glm::translate(model, glm::vec3(l_lower_arm_x, l_lower_arm_y, 0.0));
   model = glm::rotate(model, l_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 

     if(hand_move_forward == true)	
   { model = glm::rotate(model, -0.04f, glm::vec3(1.0f, 0.0f, 0.0f));	
     model = glm::translate(model, glm::vec3(0.0, 0.2, 0.0));
   }
 }

  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 

  model = mv.top();  
  mv.pop();
  mv.push(model);

   initChangeColor(4);
  ////////////////////////////////////////
  //////          3rd box right hand
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(-2.1, 0.2, 0.0)); 
  model = glm::rotate(model, r_upper_arm , glm::vec3(0.0f, 0.0f, 1.0f));
  

     if(hand_move_flag == true)
  {	
  model = glm::rotate(model, hand_move , glm::vec3(1.0f, 0.0f, 0.0f)); 

  }
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 


   initChangeColor(5);
  ///////////////////////////////////////////////////////
  //////          3rd box right hand lower arm
  //////////////////////////////////////////////////////
  
  model = mv.top();  
  mv.pop();
  mv.push(model);
    model = glm::translate(model, glm::vec3(r_lower_arm_x, r_lower_arm_y, 0.0)); 
    model = glm::rotate(model, r_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 

  if(hand_move_flag == true)
  {		
	  hand_move_flag = false;
	  printf("hand move angle: %f",hand_move);
  if (hand_move>=0.6)
  {
	  hand_move_forward = false;
  }
  else if(hand_move <=-0.6)
  {
	  hand_move_forward = true;
  }  
   model = glm::rotate(model, -r_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f));
   model = glm::translate(model, glm::vec3(-r_lower_arm_x, -r_lower_arm_y, 0.0)); 
   
   model = glm::rotate(model, hand_move, glm::vec3(1.0f, 0.0f, 0.0f));
  
   
   model = glm::translate(model, glm::vec3(r_lower_arm_x, r_lower_arm_y, 0.0));
   model = glm::rotate(model, r_lower_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 

   if(hand_move_forward == true)	
   { model = glm::rotate(model, 0.04f, glm::vec3(1.0f, 0.0f, 0.0f));	
     model = glm::translate(model, glm::vec3(0.0, 0.2, 0.0));
   }

  }



  //model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 
   initChangeColor(6);
  ///////////////////////////////////////////////
  //////          4th box left upper leg
  //////////////////////////////////////////////

  model=mv.top();
    mv.pop();
  mv.push(model);

  model = glm::translate(model, glm::vec3( 0.7, 2.2, 0.0)); 
  model = glm::rotate(model, legangle1 , glm::vec3(1.0f, 0.0f, 0.0f)); 
  //model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 

  
  ////////////////////////////////////////
  //////          4th box left lower leg
  ///////////////////////////////////////
  model=mv.top();
    mv.pop();
  mv.push(model);

  model = glm::translate(model, glm::vec3( 0.7, 3.2, 0.0)); 
  model = glm::rotate(model, lowerlegangle1 , glm::vec3(1.0f, 0.0f, 0.0f)); 
  //model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 
  
  ////////////////////////////////////////
  //////          5th box right leg
  ///////////////////////////////////////
  model = mv.top();
  model = glm::translate(model, glm::vec3( -0.7, 2.2, 0.0)); 
  model = glm::rotate(model, legangle2 , glm::vec3(1.0f, 0.0f, 0.0f)); 
  //model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glColor3f(0.0f,0.5f,1.0f); 

  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 
  
  /////////////////////////////////////////////
  //////          5th box right lower leg
  /////////////////////////////////////////////
  model = mv.top();
  model = glm::translate(model, glm::vec3( -0.7, 3.2, 0.0)); 
  model = glm::rotate(model, lowerlegangle2 , glm::vec3(1.0f, 0.0f, 0.0f)); 
  //model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  modelview = view * model; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glColor3f(0.0f,0.5f,1.0f); 
 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0); 

 mv.pop();
  model=mv.top();
//  mv.pop();
  
  //mv.pop();
  //model=mv.top();
  //mv.pop();
  //glutSwapBuffers(); 

} 
//////////////////////////////////////////////////////////////////////////////////
//
//    GLUT stuff 
//
//*

void mymouse(int button, int state, int x, int y)
{
  if (state == GLUT_DOWN) {
    press_x_cube = x; press_y_cube = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode_cube = XFORM_ROTATE; 
	 else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode_cube = XFORM_SCALE; 
  }
  else if (state == GLUT_UP) {
	  xform_mode_cube = XFORM_NONE; 
  }

//cylidnrt
    if (state == GLUT_DOWN) {
    press_x_cyl = x; press_y_cyl = y; 
    if (button == GLUT_LEFT_BUTTON)
      xform_mode_cyl = XFORM_ROTATE; 
	 else if (button == GLUT_RIGHT_BUTTON) 
      xform_mode_cyl = XFORM_SCALE; 
  }
  else if (state == GLUT_UP) {
	  xform_mode_cyl = XFORM_NONE; 
  }

}

void mymotion(int x, int y)
{
 

	    if (xform_mode_cyl==XFORM_ROTATE) {
      x_angle_rot += (x - press_x_cyl)/5.0; 
      if (x_angle_rot > 180) x_angle_rot -= 360; 
      else if (x_angle_rot <-180) x_angle_rot += 360; 
      press_x_cyl = x; 
	   rotate_x = true;
	   rotate_y = true;
      y_angle_rot += (y - press_y_cyl)/5.0; 
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

void mykey(unsigned char key, int x, int y)
{
	if(key == 'h')  
	{printf("r_upper_arm: %f, lower_arm_y: %f ", r_upper_arm , r_lower_arm_y);
		if (r_arm_up == true)
		{
		 
		 if(r_upper_arm < 1.6) 
		 { 
			r_upper_arm += 0.1f;
			r_lower_arm += 0.1f;
			r_lower_arm_x -= 0.08f;
			r_lower_arm_y -= 0.08f;			
		 }
		 else
			{r_lower_arm += 0.1f;	r_lower_arm_y -= 0.06f;
			}
		 if( r_lower_arm_y <= -1.04)
			 r_arm_up=false;
		}
	else
		{
			if( r_lower_arm_y <= -0.14)
				{
					r_lower_arm -= 0.1f;	
					r_lower_arm_y += 0.06f;
				}
			else if (r_lower_arm_y == -0.14)	r_lower_arm_y = 0.0f;
			else
			{  
				r_upper_arm -= 0.1f;
				r_lower_arm -= 0.1f;
				r_lower_arm_x += 0.08f;
				r_lower_arm_y += 0.08f;			
			}
			if(r_lower_arm_y >= 1.2f)
				r_arm_up = true;
		}
	}

	if(key == 'g')  
	{//printf("l_upper_arm: %f, l_lower_arm_y: %f ", l_upper_arm , l_lower_arm_y);
		if (l_arm_up == true)
		{
		 
		 if(l_upper_arm > -1.6) 
		 {  printf("here");
			l_upper_arm -= 0.1f;
			l_lower_arm -= 0.1f;
			l_lower_arm_x += 0.08f;
			l_lower_arm_y -= 0.08f;			
		 }
		 else
			{
				//printf("not here");
				l_lower_arm -= 0.1f;	l_lower_arm_y -= 0.06f;
			}
		 if( l_lower_arm_y <= -1.04)
			 l_arm_up=false;
		}
	else
		{
			if( l_lower_arm_y <= -0.14)
				{
					l_lower_arm += 0.1f;	
					l_lower_arm_y += 0.06f;
				}
			else if (l_lower_arm_y == 0.14)	l_lower_arm_y = 0.0f;
			else
			{  
				l_upper_arm += 0.1f;
				l_lower_arm += 0.1f;
				l_lower_arm_x -= 0.08f;
				l_lower_arm_y += 0.08f;			
			}
			if(l_lower_arm_y >= 1.2f)
				l_arm_up = true;
		}
	}

	if(key == 'w')  
		{
			hand_move_flag = true;
			if(hand_move_forward == true) {hand_move += 0.05; }
			else {hand_move-= 0.05;}
			
			if(r_leg_front == true)
			{
			 legangle1  += 0.05f;	
			 lowerlegangle1 += 0.05f;	

				 if(legangle2 <= -0.45)				
					 {//	printf("inside");
							 if( l_lowerleg_front == false) lowerlegangle2 -= 0.05f;
							 else	{ lowerlegangle2 += 0.05f; 
										if(lowerlegangle2 >= -0.45 )
										{legangle2+=0.05;	lowerlegangle2+= 0.05;	 l_lowerleg_front = false; r_leg_front = false; }
									}
							 if(l_lowerleg_front == false && lowerlegangle2 < -0.8) l_lowerleg_front = true;

					 }
				 else								//halfway left leg going back, make static upper, move lower leg
					 {  
						 legangle2 -= 0.05f; 
						 lowerlegangle2 -= 0.05f;
					 }
			}
			else		//r_leg_front = false
			{
			 		
			 legangle2 += 0.05f;
			 lowerlegangle2 += 0.05f;

			  if(legangle1 <= -0.45)				
					 {	printf("inside");
							 if( r_lowerleg_front == false) lowerlegangle1 -= 0.05f;
							 else	{ lowerlegangle1 += 0.05f; 
										if(lowerlegangle1 >= -0.45 )
										{legangle1+=0.05;	lowerlegangle1+= 0.05;	 r_lowerleg_front = false; r_leg_front = true; }
									}
							 if(r_lowerleg_front == false && lowerlegangle1 < -0.8) r_lowerleg_front = true;

					 }
				 else								//halfway left leg going back, make static upper, move lower leg
					 {  
						 legangle1 -= 0.05f; 
						 lowerlegangle1 -= 0.05f;
					 }
			}
		if(legangle1 >= 2)
			r_leg_front = false;
		if(legangle1 <= -0.2 && r_leg_front == false)
			r_leg_front = true;
		 z_coord += 0.00005f; 
		 forwardmove = true;
		 //printf("legangle2 %f, lowerlegangle2 %f \n",legangle2,lowerlegangle2);
		 //samescale = false;
		}
	if(key == 'a')	{if(x_angle_body<0.6) {x_angle_body += 0.01f; x_angle_body_flag = true;}}
	if(key == 'd')  {if(x_angle_body>-0.6) {x_angle_body -= 0.01;	x_angle_body_flag = true;}}
	if(legangle1>3) legangle1 = 0.0f;
	if(legangle2>1.74) legangle2 = 0.0f;
	if(key == 'z') if(head< 0.6) {head+= 0.02; head_flag = true;}
	if(key == 'x') if(head> -0.6) {head-= 0.02; head_flag = true;}

//	if (key == 's') draw_line = !draw_line; 
	if (key == 'q') exit(0); 
	glutPostRedisplay(); 
	
}



///////////////////////////////////////////////////////////////////////////////
//
void display()
{
	
	glClearColor(0,0,0,1); 
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
	
	display_cube();
	
//	display_cyl();
	glutSwapBuffers();
}


int main(int argc, char** argv) 
{ 
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 

  glutCreateWindow("shader cube"); 
  
  glutDisplayFunc(display);
  //glutDisplayFunc(display_cube); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 

  // initialize GLEW 

  //  GLenum err = glewInit(); 
  //  if ( err != GLEW_OK)  printf(" Error initializing GLEW! \n"); 
  //  else printf("Initializing GLEW succeeded!\n"); 
  GLenum err = glewInit(); 

  if ( err != GLEW_OK)  printf(" Error initializing GLEW! \n"); 
  else printf("Initializing GLEW succeeded!\n"); 

  // define the discretion level for the cylinder 
  int nslices, nstacks; 
  nslices = 60; 
  nstacks = 20; 

  InitCylinder(nslices, nstacks, 1.0, 0.0, 1.0); 

  glewInit(); 

  InitGeometry();
  InitVBO_cyl(nslices,nstacks);

  nslices = 40; 
  nstacks = 40; 

  which_sph_color = 0; 


  //  which_sph_color = 2; // 0:T, 1:N, 2:B
  InitSphere(nslices, nstacks, 1.0, 0.0, 1.0); 
  InitVBO_sphere(nslices, nstacks); 

  InitVBO_cube(); 
  glutMainLoop(); 
} 
