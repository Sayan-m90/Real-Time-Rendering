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


GLuint programObject;
GLuint SetupGLSL(char*);
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
Vertexsphere *sphverts;   // sphere vertices
GLuint *sp_indices; 

GLuint sphere_vboHandle[1];   // a VBO that contains interleaved positions and sph_colors 
GLuint sphere_indexVBO; 

///////////////////////////////////////////////////////////////
// Define Light Properties -  Ia, Id, Is, and light position 
////////////////////////////////////////////////////////////////

GLfloat light_ambient[4] = {1,1,1,1};  //Ia 
GLfloat light_diffuse[4] = {0.4,0.8,0.8,1};  //Id
GLfloat light_specular[4] = {1,1,1,1};  //Is
GLfloat light_pos [4] = {4, -3.0, 8.8, 1};

/////////////////////////////////////////////////////////////////
// Define Default Material Properties -  Ka, Kd, Ks, Shininess 
////////////////////////////////////////////////////////////////

GLfloat mat_ambient[4] = {0.1,0.1,0.1,1};  //Ka 
GLfloat mat_diffuse[4] = {0.8,0.8,0,1};  //Kd
GLfloat mat_specular[4] = {1,1,1,1};  //Ks
GLfloat mat_shine[1] = {10}; 

////////////////////////////////////////////////
typedef struct 
{
  float location[4]; 
  float normal[4]; 
//  float cube_color[4]; 
} Vertexcube; 
//
Vertexcube cverts[24];       // cube vertices 



GLubyte tindices[24];   // cube face indices: 6 faces, 6 indices each face 

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
float shoulder_move = 0.0;
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
// create VBO objects and send the triangle vertices/cube_colors to the graphics card
// 
void InitCube()
{ float size = 1.0;
	

	cverts[0].location[0] = cverts[1].location[0] = cverts[2].location[0] = cverts[3].location[0] = -0.5*size;
	cverts[4].location[0] = cverts[5].location[0] = cverts[14].location[0] = cverts[15].location[0] = -0.5*size;
	cverts[17].location[0] = cverts[19].location[0] = cverts[21].location[0] = cverts[23].location[0] = -0.5*size;


	cverts[0].location[1] = cverts[2].location[1] = cverts[4].location[1] = cverts[6].location[1] = -0.5*size;
	cverts[8].location[1] = cverts[10].location[1] = cverts[12].location[1] = cverts[14].location[1] = -0.5*size;
	cverts[16].location[1] = cverts[17].location[1] = cverts[18].location[1] = cverts[19].location[1] = -0.5*size;


	cverts[0].location[2] = cverts[1].location[2] = cverts[10].location[2] = cverts[11].location[2] = -0.5*size;
	cverts[12].location[2] = cverts[13].location[2] = cverts[14].location[2] = cverts[15].location[2] = -0.5*size;
	cverts[18].location[2] = cverts[19].location[2] = cverts[22].location[2] = cverts[23].location[2] = -0.5*size;

	
	cverts[6].location[0] = cverts[7].location[0] = cverts[8].location[0] = cverts[9].location[0] = 0.5*size;
	cverts[10].location[0] = cverts[11].location[0] = cverts[12].location[0] = cverts[13].location[0] = 0.5*size;
	cverts[16].location[0] = cverts[18].location[0] = cverts[20].location[0] = cverts[22].location[0] = 0.5*size;


	cverts[1].location[1] = cverts[3].location[1] = cverts[5].location[1] = cverts[7].location[1] = 0.5*size;
	cverts[9].location[1] = cverts[11].location[1] = cverts[13].location[1] = cverts[15].location[1] = 0.5*size;
	cverts[20].location[1] = cverts[21].location[1] = cverts[22].location[1] = cverts[23].location[1] = 0.5*size;

	cverts[2].location[2] = cverts[3].location[2] = cverts[4].location[2] = cverts[5].location[2] = 0.5*size;
	cverts[6].location[2] = cverts[7].location[2] = cverts[8].location[2] = cverts[9].location[2] = 0.5*size;
	cverts[16].location[2] = cverts[17].location[2] = cverts[20].location[2] = cverts[21].location[2] = 0.5*size;
	
	for (int i = 0; i < 24; i++)
	{
		cverts[i].location[3] = 1.0;
		cverts[i].normal[3] = 0.0;
		
		tindices[i] = i;
	}

	for (int i = 0; i < 6; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			float x, y, z;
			if (i == 0 || i == 2)
			{
				y = 0.0; z = 0.0;
				if (i == 0)
					x = -1.0;
				else
					x = 1.0;
			}

			else if (i == 1 || i ==3)
			{
				x = 0.0; y = 0.0;
				if (i == 1)
					z = 1.0;
				else
					z = -1.0;
			}
			else 
			{
				x = 0.0; z = 0.0;
				if (i == 4)
					y = -1.0;
				else
					y = 1.0;
			}
			cverts[i * 4 + j].normal[0] = x;
			cverts[i * 4 + j].normal[1] = y;
			cverts[i * 4 + j].normal[2] = z;		
		}
	}
	


}

void InitVBO_cube() 
{InitCube();

  glGenBuffers(1, &cube_vboHandle[0]);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(Vertexcube)*24, cverts, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &cube_indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint)*24, tindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

} 


 

bool first = true;


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


void draw_cube(float* local2clip, float* local2eye, float* world2eye, float color[3],  GLuint c0, GLuint c1, GLuint c2, GLuint m1, GLuint m2, GLuint m3, GLuint m4)
{

  glBindBuffer(GL_ARRAY_BUFFER, cube_vboHandle[0]); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, cube_indexVBO);

  glVertexAttribPointer(c0,4,GL_FLOAT, GL_FALSE, sizeof(Vertexcube),(char*) NULL+0);  // position 
  glVertexAttribPointer(c2,4,GL_FLOAT, GL_FALSE, sizeof(Vertexcube),(char*) NULL+16); // normal 

  glm::mat4 normal_matrix = glm::inverse((glm::mat4)*local2eye);
  normal_matrix = glm::transpose(normal_matrix);

  glUniformMatrix4fv(m1, 1, GL_FALSE, (float*) local2clip);   // pass the local2clip matrix
  glUniformMatrix4fv(m2, 1, GL_FALSE, (float*) local2eye);   // pass the local2eye matrix
  glUniformMatrix4fv(m3, 1, GL_FALSE, (float*) &normal_matrix[0][0]);   // pass the local2eye matrix
  glUniformMatrix4fv(m4, 1, GL_FALSE, (float*) world2eye);   // pass the local2eye matrix

  glDrawElements(GL_TRIANGLE_STRIP, 24, GL_UNSIGNED_BYTE, (char*)NULL+0); 

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


//////////////////////////////////////////////////////////////

#define SetMaterialColor(d, r, g, b, a)  glUniform4f(d, r, g, b, a); 

/////////////////////////////////////////////////////////////



void display() 
{ 
	glClearColor(0,0,0,1); 
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT); 
	
  glEnable(GL_DEPTH_TEST);    // need depth test to correctly draw 3D objects 
   glUseProgram(programObject);
   glLightModelf(GL_LIGHT_MODEL_TWO_SIDE, GL_TRUE); 
  /////////////////////////////////////////////
  ////////////////////////    Shader 
  /////////////////////////////////////////////

  // glUseProgram(programObject);

  GLuint c0 = glGetAttribLocation(programObject, "position");
  GLuint c1 = glGetAttribLocation(programObject, "color");
  GLuint c2 = glGetAttribLocation(programObject, "normal");
  GLuint m1 = glGetUniformLocation(programObject, "local2clip");
  GLuint m2 = glGetUniformLocation(programObject, "local2eye");
  GLuint m3 = glGetUniformLocation(programObject, "normal_matrix");
  GLuint m4 = glGetUniformLocation(programObject, "world2eye");  

  GLuint Ia = glGetUniformLocation(programObject, "light_ambient");
  GLuint Id = glGetUniformLocation(programObject, "light_diffuse");
  GLuint Is = glGetUniformLocation(programObject, "light_specular");
  GLuint Lpos = glGetUniformLocation(programObject, "light_pos");

  GLuint Ka = glGetUniformLocation(programObject, "mat_ambient");
  GLuint Kd = glGetUniformLocation(programObject, "mat_diffuse");
  GLuint Ks = glGetUniformLocation(programObject, "mat_specular");
  GLuint Shine = glGetUniformLocation(programObject, "mat_shine"); 

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


  
  /////////////////////////////////////////////
  ///////////////////////////   End of Shader
  ////////////////////////////////////////////

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
	//glutPostRedisplay();


	//////////////////////////////////////////////////
  glm::mat4 mvp;
  glm::mat4 mv;
  float toradian = 197.74*3.141/180;
  glm::mat4 projection = glm::perspective(toradian,1.0f,.1f,100.0f); 
  //glMatrixMode(GL_PROJECTION);
  //glLoadMatrixf(&projection[0][0]); 

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
	  x_angle_rot = 0.0f;
  }

  if(rotate_y == true) 
	{model = glm::rotate(model, y_angle_rot, glm::vec3(0.0f, 1.0f, 0.0f)); 
	 model2 = glm::rotate(model2, y_angle_rot, glm::vec3(0.0f, 1.0f, 0.0f)); 
	 rotate_y = false;
	 y_angle_rot = 0.0f;
    }
   
   if(samescale==false) 
	  { model = glm::scale(model, glm::vec3(scale_size_cyl, scale_size_cyl, scale_size_cyl)); 
		model2 = glm::scale(model, glm::vec3(scale_size_cyl, scale_size_cyl, scale_size_cyl)); 
		samescale = true;
		}
   if(forwardmove == true) 
   {model = glm::translate(model, glm::vec3(0.0, 0.0, z_coord)); forwardmove = false; z_coord = 0;}
   if(x_angle_body_flag == true)
   {
	   model = glm::rotate(model, x_angle_body , glm::vec3(0.0f, 1.0f, 0.0f)); 
	   x_angle_body_flag = false;
   }
   
   st.push(model);


  if (draw_line) glPolygonMode(GL_FRONT_AND_BACK,GL_LINE); 
	  else glPolygonMode(GL_FRONT_AND_BACK,GL_FILL); 
		 float color[3];

	    glm::mat4 lightM = model2; 
  lightM = glm::translate(lightM, glm::vec3(light_pos[0], light_pos[1], light_pos[2]-6.0f)); 
  lightM = glm::scale(lightM, glm::vec3(0.3, 0.3, 0.3));

 mvp = projection*view*lightM;
 mv = view*lightM;

  
  SetMaterialColor(Kd, 0, 0, 0, 1); 
  draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3, m4);



	  //////////////////////////////////////
	  ////////////////////   Sphere
	  //////////////////////////////////////
  
  model = glm::translate(model, glm::vec3(0.0, -3.0, -0.1f));   
  model = glm::scale(model, glm::vec3(0.2,0.2,0.2)); 

  glm::mat4 temp = glm::mat4(1.0f); 
  temp = model;
   if(head_flag == true) 
   {temp = glm::rotate(model, head, glm::vec3(0.0,1.0,0.0));	}//	head_flag = false;	}

//   modelview = view*model;
 //  modelview1 = view * temp; 
//  glMatrixMode(GL_MODELVIEW); 
 // glLoadMatrixf(&modelview1[0][0]); 
     mvp = projection*view*model;
     mv = view * model;

	 color[0] = 1.0;   color[1] = 0.7;   color[2] = 0.3;
	 SetMaterialColor(Kd, 0.7, 0, 0.7, 1); 
     draw_sphere(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 


  //glDrawElements(GL_TRIANGLES, sphindices, GL_UNSIGNED_INT, (char*) NULL+0); 

  //glDisableClientState(GL_VERTEX_ARRAY); // enable the Vertexsphere array on the client side
  //glDisableClientState(GL_COLOR_ARRAY); // enable the sph_color array on the client side
  
  ////////////////////////////////////////////////////// 
  ////////////////////    Cylinder    /////////////////
  ////////////////////////////////////////////////////
  model=st.top();
  
   model = glm::rotate(model, 1.57f, glm::vec3(1.0f, 0.0f, 0.0f)); 
   model = glm::scale(model, glm::vec3(1.0f, 1.0f, 3.0f)); 
   model = glm::translate(model,glm::vec3(0.0f,0.0f,-0.8f));
   
   

   model1 = model;
   model1 = glm::scale(model1, glm::vec3(1.1f, 1.0f, 0.4f));

 // modelview = view * model; 
  //glMatrixMode(GL_MODELVIEW); 
 
  
   color[0] = 1.0;   color[1] = 1.0;   color[2] = 0.0;
     mvp = projection*view*model;
  mv = view * model;
  SetMaterialColor(Kd, 0.7, 0, 0.7, 1); 
   draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 

   mvp = projection*view*model1;
   mv = view * model1;
   SetMaterialColor(Kd, 0.7, 0.4, 0.9, 1); 
   draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  
   
   model3 = glm::rotate(model2, 1.57f, glm::vec3(1.0f, 0.0f, 0.0f)); 
   model3 = glm::translate(model3, glm::vec3(-8.5, 3.2 , -3.5)); 
  model3 = glm::scale(model3, glm::vec3(0.3f, 0.1f, 2.1f)); 
  glm::mat4  modell = model;
  mvp = projection*view*model3;
   mv = view * model3;
   SetMaterialColor(Kd, 0.4, 0.9, 0.4, 1); 
   draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3,m4);

   model3  = glm::scale(model3, glm::vec3(0.09f, 0.09f, 1.5f)); 
   //model3 = glm::rotate(model2, 0.0014f, glm::vec3(1.0f, 1.0f, 0.0f)); 
   mvp = projection*view*model3;
   mv = view * model3;
   SetMaterialColor(Kd, 0.9, 0.4, 0.4, 1); 
   draw_cylinder(mvp, mv, &view[0][0], color, c0, c1, c2, m1, m2, m3,m4);

  ///////////////////////////////////////////
  ///////////// Cubes          /////////////
  /////////////////////////////////////////
  

  /////////////////////////////////////////
  //////          base    ////////////////
  ///////////////////////////////////////
  //model=st.top();
  model3 = glm::translate(model2, glm::vec3(0.0,8.0 , 0.0)); 
  model3 = glm::scale(model3, glm::vec3(50.0f, 0.20f, 10.0f)); 

  /*modelview = view * model3; 
  glMatrixMode(GL_MODELVIEW); 
  glLoadMatrixf(&modelview[0][0]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*) NULL+0);
  */
  mvp = projection*view*model3;
   mv = view * model3;
	
	// color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
		 SetMaterialColor(Kd, 0.9, 0.78, 0.49, 1); 
	 draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

  model=st.top();

  //initChangeColor(0);
  ////////////////////////////////////////
  //////          base2
  ///////////////////////////////////////
  //model=st.top();
  model3 = glm::translate(model2, glm::vec3(-1.7,8.0 , 0.0)); 
  model3 = glm::scale(model3, glm::vec3(10.0f, 0.20f, 10.0f)); 


    mvp = projection*view*model3;
     mv = view * model3;
	
	 color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	 
   // draw_cube(&mvp[0][0], &mv[0][0], color, c0, c1, c2, m1, m2, m3); 

  model=st.top();

  ////////////////////////////////////////
  //////          back
  ///////////////////////////////////////
  //model3 = glm::translate(model2, glm::vec3(0.0,0.0 , -10.0)); 
  model3 = glm::scale(model2, glm::vec3(50.0f, 30.0f, 0.2f)); 


    mvp = projection*view*model3;
     mv = view * model3;
	
	 color[0] = 1.0;   color[1] = 0.7;   color[2] = 0.3;
SetMaterialColor(Kd, 0.99, 0.23, 0.31, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

  //model=st.top();

  model=st.top();
 
  
  //initChangeColor(8);

  ////////////////////////////////////////
  //////          books
  ///////////////////////////////////////

  ////////////////////////// one /////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-5.4,2.0 , 3.2)); 
  model3 = glm::scale(model3, glm::vec3(0.3f, 1.2f, 0.1f)); 

 // modelview = view * model3; 
 //glMatrixMode(GL_MODELVIEW); 
 // glLoadMatrixf(&modelview[0][0]); 

 
   mvp = projection*view*model3;
     mv = view * model3;
	
	 color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	 SetMaterialColor(Kd, 0.0, 0.55, 0.3, 1); 
       draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  ////////////////////		two		//////////////////////////
  //model=st.top();
 //initChangeColor(7);
  model3 = glm::translate(model2, glm::vec3(-6.3,1.6 , 3.2)); 
  model3 = glm::scale(model3, glm::vec3(0.5f, 2.0f, 0.1f)); 
	
	mvp = projection*view*model3;
    mv = view * model3;
	//color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	SetMaterialColor(Kd, 0.5, 0.55, 0.0, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  ///initChangeColor(9);
  ///////////////////		three		/////////////////////////////
  //model=st.top();
  model3 = glm::translate(model2, glm::vec3(-7.0,2.0 , 3.2)); 
  model3 = glm::scale(model3, glm::vec3(0.8f, 1.2f, 0.1f)); 

	mvp = projection*view*model3;
    mv = view * model3;
	color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	SetMaterialColor(Kd, 0.99, 0.13, 0.9, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

  //model=st.top();
  ////////////////////////////////////////
  ////////////  TV	    /////////////////
  //////////////////////////////////////
	
	model3 = glm::translate(model2, glm::vec3(-7.0,-1.97 , 1.2)); 
	model3 = glm::scale(model3, glm::vec3(4.8f, 3.2f, 0.1f)); 

	mvp = projection*view*model3;
    mv = view * model3;
	color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	SetMaterialColor(Kd, 0.1, 0.1, 0.1, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  

	model3 = glm::translate(model2, glm::vec3(-6.85,-1.95 , 1.25)); 
	model3 = glm::scale(model3, glm::vec3(4.0f, 2.4f, 0.1f)); 

	mvp = projection*view*model3;
    mv = view * model3;
	color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	SetMaterialColor(Kd, 0.6, 0.6, 0.6, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

	////////////////////////////////////////
  //////         boxes
  ///////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(7.5, 5.2 , 3.5)); 
  model3 = glm::scale(model3, glm::vec3(2.0f, 2.0f, 0.1f)); 
  

   mvp = projection*view*model3;
   mv = view * model3;
	
	 color[0] = 0.0;   color[1] = 0.7;   color[2] = 0.3;
	SetMaterialColor(Kd, 0.99, 0.13, 0.9, 1);  
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  model=st.top();
  //initChangeColor(2);
  ////////////////////////////////////////
  //////          table
  ///////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-8, 3.2 , 2.5)); 
  model3 = glm::scale(model3, glm::vec3(4.8f, 0.4f, 1.8f)); 
  mvp = projection*view*model3;
  mv= view * model3; 
  SetMaterialColor(Kd, 0.49, 0.31, 0.0, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
 
  //////////////////////////////////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-7,5.0 , 2.3)); 
  model3 = glm::scale(model3, glm::vec3(1.0f, 2.3f, 1.5f)); 

  //modelview = view * model3; 
   mvp = projection*view*model3;
   mv= view * model3; 
   SetMaterialColor(Kd, 0.13, 0.1, 0.0, 1); 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
 
  ////////////////////////////////////////////////////////////////////////
  model3 = glm::translate(model2, glm::vec3(-11, 5.0 , 2.3)); 
  model3 = glm::scale(model3, glm::vec3(1.0f, 2.3f, 1.5f)); 

 // modelview = view * model3; 
   mvp = projection*view*model3;
  mv= view * model3; 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
 
  //////////////////////////////////////////////////////////////
  ///////////////////		chair		///////////////////////
  ////////////////////////////////////////////////////////////
  //initChangeColor(2);
  model3 = glm::translate(model2, glm::vec3(8.0, 2.0, 2.9)); 
  model3 = glm::scale(model3, glm::vec3(1.5f, 4.0f, 0.05f)); 

  //modelview = view * model3; 
   mvp = projection*view*model3;
  mv= view * model3; 
   
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
 
  /////////////////////////////////////////////////
 
	model3 = glm::translate(model2, glm::vec3(8.0, 3.0, 3.26)); 
    model3 = glm::scale(model3, glm::vec3(1.5f, 1.6f, 0.05f)); 
    mvp = projection*view*model3;
    mv= view * model3; 
	SetMaterialColor(Kd, 0.13, 0.1, 0.0, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
 // modelview = view * model3; 
 
 // glLoadMatrixf(&modelview[0][0]); 
 

    /////////////////////////////////////////////////
  //initChangeColor(3);
  model3 = glm::translate(model2, glm::vec3(8.0, 2.0, 3.05)); 
  model3 = glm::scale(model3, glm::vec3(1.8f, 0.30f, 0.5f)); 
   mvp = projection*view*model3;
  mv= view * model3; 
SetMaterialColor(Kd, 0.49, 0.31, 0.0, 1); 
    draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  //modelview = view * model3; 
 
 // glLoadMatrixf(&modelview[0][0]); 
 


  model=st.top();
  //initChangeColor(4);
  ////////////////////////////////////////
  //////          Head
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0,-2.3 , 0.1)); 
   model = glm::scale(model, glm::vec3(1.0, 1.0 , 1.2));
   temp =model;
  if(head_flag == true) 
  {temp = glm::rotate(model, head, glm::vec3(0.0,1.0,0.0));}//		head_flag = false;	}
   mvp = projection*view*temp;
  
   mv= view * temp; 
   SetMaterialColor(Kd, 0.69, 0.51, 0.0, 1); 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

  
  model=st.top();
 

   ////////////////////////////////////////
  //////          neck
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0,-1.8 , 0.0)); 
  model = glm::scale(model,glm::vec3(0.4,0.4,0.4));

//  modelview = view * model; 
    mvp = projection*view*model;
   mv= view * model; 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
  

   model=st.top();
 SetMaterialColor(Kd, 0.0, 0.0, 0.00, 1); 
    //initChangeColor(3);
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
//  modelview1 = view * temp; 
   mvp = projection*view*temp;
   mv= view * temp; 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 

 
    //initChangeColor(3);
   ////////////////////////////////////////
  //////          eyes2
  ///////////////////////////////////////
   model = st.top();
  model = glm::translate(model, glm::vec3(-0.25, -2.4 , 0.8)); 
  model = glm::scale(model,glm::vec3(0.20,0.25,0.10));
  temp = model;
   if(head_flag == true) 
 {
	 temp = glm::translate(temp, glm::vec3(0.25, 2.4 , -0.8));
	 temp = glm::rotate(temp, head, glm::vec3(0.0,1.0,0.0));		
	 temp = glm::translate(temp, glm::vec3(-0.25, -2.4 , 0.8));}
	 //modelview1 = view * temp; 
   mvp = projection*view*temp;
   mv= view * temp; 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
	 //initChangeColor(3);
	 	model=st.top();

 

   ////////////////////////////////////////
  //////          mouth
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(0.0, -2.0 , 0.8)); 
  model = glm::scale(model,glm::vec3(0.20,0.25,0.10));
  temp = model;
  if(head_flag == true) 
 {
	  temp = glm::translate(temp, glm::vec3(0.0, 2.0 , -0.8)); 
	  temp = glm::rotate(temp, head, glm::vec3(0.0,1.0,0.0));		
	  temp = glm::translate(temp, glm::vec3(0.0, -2.0 , 0.8));
	  head_flag = false;	}
	 
 // modelview1 = view * temp; 
   mvp = projection*view*temp;
   mv= view * model; 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3,m4); 
	model = st.top();
	 //initChangeColor(4);
  ////////////////////////////////////////
  //////          shoulder
  ///////////////////////////////////////
	SetMaterialColor(Kd, 0.69, 0.51, 0.0, 1); 
  if(hand_move_forward == true && hand_move_flag == true)
	  { shoulder_move += 0.01f;
		model = glm::translate(model,glm::vec3(0.0,shoulder_move,0.0));
		}
else if(hand_move_forward == false && hand_move_flag == true)
	 { shoulder_move -= 0.01f;
  	   model = glm::translate(model,glm::vec3(0.0,shoulder_move,0.0));
	}
 
   model = glm::translate(model,glm::vec3(0.0,-1.1,0.0));
   model = glm::scale(model,glm::vec3(4.9,0.8,1.0));
 
 // modelview = view * model; 
   mvp = projection*view*model;
   mv= view * model; 
   draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 
  //model=st.top();
//  model = glm::scale(model, glm::vec3(0.8, 2.0, 0.8)); 
  //st.push(model);
  model=st.top();

  
  ////////////////////////////////////////
  //////          2nd box left hand
  ///////////////////////////////////////
  
   model = glm::translate(model, glm::vec3(2.1, 0.2, 0.0)); 
   model = glm::rotate(model, l_upper_arm , glm::vec3(0.0f, 0.0f, 1.0f)); 
  
    if(hand_move_flag == true)
  {		
	model = glm::rotate(model, -hand_move , glm::vec3(1.0f, 0.0f, 0.0f)); 
  }
  
     mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 0.3); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 

  //model=st.top();



  model = st.top();  
  st.pop();
  st.push(model);

   //initChangeColor(5);
  //////////////////////////////////////////////////////
  //////          2nd box left hand lower arm
  ////////////////////////////////////////////////////
  model = glm::translate(model, glm::vec3(l_lower_arm_x, l_lower_arm_y, 0.0)); 

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


 mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 0.3); 
     draw_cube(&mvp[0][0], &mv[0][0],&view[0][0],  color, c0, c1, c2, m1, m2, m3, m4); 

  model = st.top();  
  st.pop();
  st.push(model);

   //initChangeColor(4);
  ////////////////////////////////////////
  //////          3rd box right hand
  ///////////////////////////////////////
  model = glm::translate(model, glm::vec3(-2.1, 0.2, 0.0)); 
  model = glm::rotate(model, r_upper_arm , glm::vec3(0.0f, 0.0f, 1.0f));
  

     if(hand_move_flag == true)
  {	
  model = glm::rotate(model, hand_move , glm::vec3(1.0f, 0.0f, 0.0f)); 

  }

	 mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 

   //initChangeColor(5);
  ///////////////////////////////////////////////////////
  //////          3rd box right hand lower arm
  //////////////////////////////////////////////////////
  
  model = st.top();  
  st.pop();
  st.push(model);
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

     mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 1); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 

   //initChangeColor(6);
  ///////////////////////////////////////////////
  //////          4th box left upper leg
  //////////////////////////////////////////////

  model=st.top();
//  st.pop();
//  st.push(model);

  model = glm::translate(model, glm::vec3( 0.7, 3.5, 0.0)); 

   if(leg_move_flag == true)
  {		
	model = glm::rotate(model, -leg_move , glm::vec3(1.0f, 0.0f, 0.0f)); 
  }

     mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 1); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 
  
  //////////////////////////////////////////////
  //////          4th box left lower leg
  /////////////////////////////////////////////
  //initChangeColor(5);
  model=st.top();
   // st.pop();
 // st.push(model);

  model = glm::translate(model, glm::vec3( 0.7, 4.2, 0.0)); 

   if(leg_move_flag == true)
  {		 

   model = glm::translate(model, glm::vec3(0.0, -1.0, 0.0)); 
   
   model = glm::rotate(model, -leg_move, glm::vec3(1.0f, 0.0f, 0.0f));
  

   
   model = glm::translate(model, glm::vec3(0.0, 1.0 , 0.0));

 }
  model = glm::scale(model,glm::vec3(1.0,2.5,1.0));
     mvp = projection*view*model;
     mv = view * model;
	 color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 1); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0],  color, c0, c1, c2, m1, m2, m3, m4); 
  
  ////////////////////////////////////////
  //////          5th box right leg
  ///////////////////////////////////////
  model = st.top();
  //initChangeColor(6);
  model = glm::translate(model, glm::vec3( -0.7, 3.5, 0.0)); 
 
   if(leg_move_flag == true)
  {		
	model = glm::rotate(model, leg_move , glm::vec3(1.0f, 0.0f, 0.0f)); 

	if(leg_move_forward == true)	
   { model = glm::rotate(model, -0.04f, glm::vec3(1.0f, 0.0f, 0.0f));	
     model = glm::translate(model, glm::vec3(0.0, -0.2, 0.0));
   }
  }

     mvp = projection*view*model;
     mv = view * model;
      color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	// SetMaterialColor(Kd, 1.0, 0, 0.0, 1); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4); 
  //initChangeColor(5);
  /////////////////////////////////////////////
  //////          5th box right lower leg
  /////////////////////////////////////////////

	 model = st.top();
  //initChangeColor(6);
  model = glm::translate(model, glm::vec3( -0.7, 4.2, 0.0)); 
 
   if(leg_move_flag == true)
  {		
	  leg_move_flag = false;
	  printf("leg move angle: %f",hand_move);
  if (leg_move>=0.6)
  {
	  leg_move_forward = false;
  }
  else if(leg_move <=-0.6)
  {
	  leg_move_forward = true;
  }  

   model = glm::translate(model, glm::vec3(0.0, -1.0, 0.0)); 
   
   model = glm::rotate(model, leg_move, glm::vec3(1.0f, 0.0f, 0.0f));
  
   
   model = glm::translate(model, glm::vec3(0.0, 1.0, 0.0));
 
   if(leg_move_forward == false)	
   { model = glm::rotate(model, -0.04f, glm::vec3(1.0f, 0.0f, 0.0f));	
     model = glm::translate(model, glm::vec3(0.0, -0.2, 0.0));
   }

  }

     model = glm::scale(model,glm::vec3(1.0,2.5,1.0));
     mvp = projection*view*model;
     mv = view * model;
      color[0] = 1.0;   color[1] = 0.0;   color[2] = 0.3;
	 //SetMaterialColor(Kd, 1.0, 0, 0.0, 1); 
     draw_cube(&mvp[0][0], &mv[0][0], &view[0][0], color, c0, c1, c2, m1, m2, m3, m4);



//st.pop();
  model=st.top();
  glDisableClientState(GL_VERTEX_ARRAY); 
  glutSwapBuffers();

}	// end of function display

//////////////////////////////////////////////////////////////////////////////////
//
//    GLUT stuff 
//
//*

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

//cylidnrt
  //  if (state == GLUT_DOWN) {
  //  
  //  if (button == GLUT_LEFT_BUTTON)
  //    
	 //else if (button == GLUT_RIGHT_BUTTON) 
  //    
  //}
  //else if (state == GLUT_UP) {
	 //
  //}

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
	{
		if (l_arm_up == true)
		{
		 
		 if(l_upper_arm > -1.6) 
		 {  //printf("here");
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
			leg_move_flag = true;
			if(leg_move_forward == true) {leg_move += 0.05; }
			else {leg_move-= 0.05;}
			

		 z_coord += 0.005f; 
		 forwardmove = true;

		}
		if(key == 's')  
		{
			hand_move_flag = true;
			if(hand_move_forward == true) {hand_move += 0.05; }
			else {hand_move-= 0.05;}
			leg_move_flag = true;
			if(leg_move_forward == true) {leg_move += 0.05; }
			else {leg_move-= 0.05;}
			

			z_coord = z_coord - 0.005f; 
		 forwardmove = true;

		}
	if(key == 'a')	{if(x_angle_body<0.6) {x_angle_body = 0.01f;	x_angle_body_flag = true;}}
	if(key == 'd')  {if(x_angle_body>-0.6) {x_angle_body = -0.01;	x_angle_body_flag = true;}}
	//if(legangle1>3) legangle1 = 0.0f;
	//if(legangle2>1.74) legangle2 = 0.0f;
	if(key == 'z') if(head< 0.5) {head+= 0.02; head_flag = true;} else{head = 0.0;}
	if(key == 'x') if(head> -0.5) {head-= 0.02; head_flag = true;}else{head = 0.0;}
	head_flag = true;
		if (key == '+') 
          mat_shine[0] += 1; 	
		if (key == '-') 
          mat_shine[0] -= 1; 	
        if (key == '3') 
	  light_pos[0] += 0.5; 
        if (key == '4') 
	  light_pos[0] -= 0.5; 
        if (key == '5') 
	  light_pos[1] += 0.5; 
        if (key == '6') 
	  light_pos[1] -= 0.5; 
        if (key == '7') 
	  light_pos[2] += 0.5; 
        if (key == '8') 
	  light_pos[2] -= 0.5;
		if(key == ' ')
			light_ambient[0]=light_ambient[1]=light_ambient[2]=light_ambient[0]+1;
	if (key == 'q') exit(0); 
	glutPostRedisplay(); 
	
}



///////////////////////////////////////////////////////////////////////////////
//

//int times,timed;
void animation(void)
{

 timed = time(NULL);
 if(timed-times >= 3)
	 {
		float temp;
		//temp = 0.2;
		temp = light_diffuse[0];
		 light_diffuse[0] = light_diffuse[1];  //Id
		 light_diffuse[1] = temp;
		 times=time(NULL);
	 }
	// xTargetRotate += 0.005;
	glutPostRedisplay();
}

int main(int argc, char** argv) 
{ 
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 

  glutCreateWindow("shader cube"); 
  times=time(NULL);
  glutDisplayFunc(display);
  //glutDisplayFunc(display); 
  glutMouseFunc(mymouse); 
  glutMotionFunc(mymotion);
  glutKeyboardFunc(mykey); 

  // initialize GLEW 
 // glutIdleFunc(animation);
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
  model = glm::translate(model,glm::vec3(0,1.0,0));
  InitCylinder(nslices, nstacks, 1.0, 0.0, 1.0); 

  glewInit(); 

 // InitGeometry();
  InitVBO_cyl(nslices,nstacks);

  nslices = 40; 
  nstacks = 40; 

  which_sph_color = 0; 


  //  which_sph_color = 2; // 0:T, 1:N, 2:B
  InitSphere(nslices, nstacks, 1.0, 0.0, 1.0); 
  InitVBO_sphere(nslices, nstacks); 
 // programObject = glCreateProgram();
  InitVBO_cube(); 
  programObject = SetupGLSL("robot");
  glutMainLoop(); 
} 
