
#include <GL/glew.h> 
#include <GL/glut.h> 
#include <GL/gl.h>
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <stack> 
#include <math.h> 
#include <time.h>
#include <Windows.h>

#ifdef __APPLE__
#include <GLUT/glut.h> 
#include <OpenGL/gl.h> 
#else 
#include <GL/glut.h> 
#include <GL/glew.h> 
#include <GL/gl.h> 
#endif

#include<glm/glm.hpp> 
#include<glm/gtx/transform2.hpp>
#include<glm/gtx/transform.hpp> 


using namespace std; 

float vertices[] = { 0.08f, 0.08f, 0.08f,  -0.08f, 0.08f, 0.08f,  -0.08f,-0.08f, 0.08f,   0.08f,-0.08f, 0.08f,   // v0,v1,v2,v3 (front)
                        0.08f, 0.08f, 0.08f,   0.08f,-0.08f, 0.08f,   0.08f,-0.08f,-0.08f,   0.08f, 0.08f,-0.08f,   // v0,v3,v4,v5 (right)
                        0.08f, 0.08f, 0.08f,   0.08f, 0.08f,-0.08f,  -0.08f, 0.08f,-0.08f,  -0.08f, 0.08f, 0.08f,   // v0,v5,v6,v1 (top)
                       -0.08f, 0.08f, 0.08f,  -0.08f, 0.08f,-0.08f,  -0.08f,-0.08f,-0.08f,  -0.08f,-0.08f, 0.08f,   // v1,v6,v7,v2 (left)
                      -0.08f,-0.08f,-0.08f,   0.08f,-0.08f,-0.08f,   0.08f,-0.08f, 0.08f,  -0.08f,-0.08f, 0.08f,   // v7,v4,v3,v2 (bottom)
                        0.08f,-0.08f,-0.08f,  -0.08f,-0.08f,-0.08f,  -0.08f, 0.08f,-0.08f,   0.08f, 0.08f,-0.08f }; // v4,v7,v6,v5 (back)
 
GLubyte tindices[]  = { 0, 1, 2,   2, 3, 0,      // front
                       4, 5, 6,   6, 7, 4,      // right
                       8, 9,10,  10,11, 8,      // top
                      12,13,14,  14,15,12,      // left
                     16,17,18,  18,19,16,      // bottom
                      20,21,22,  22,23,20 };    // back


GLuint vboHandle[1];   // a VBO that contains interleaved positions and colors 
GLuint indexVBO; 
int times,timed;
bool won = false;
bool lost =  false;
float angle1=0, angle2=0;
float angle3=0, angle4=0;
float gun_angle = 0;
bool bullet = false;

glm::mat4 modelM      = glm::mat4(1.0f); 
glm::mat4 modelm      = glm::mat4(1.0f); 
glm::mat4 modeltar    = glm::mat4(1.0f); 
glm::mat4 modeltars   = glm::mat4(1.0f); 
glm::mat4 modelnozzle = glm::mat4(1.0f); 
vector<glm::mat4> mat_list;
stack<glm::mat4> mat_stack; 
stack<glm::mat4> tar_stack; 

void InitVBO() 
{
  glGenBuffers(1, vboHandle);   // create an interleaved VBO object
  glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);   // bind the first handle 

  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*72, vertices, GL_STATIC_DRAW); // allocate space and copy the position data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 
  glGenBuffers(1, &indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLubyte)*36, tindices, GL_STATIC_DRAW);  // load the index data 

  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,0);  // clean up 

} 

void drawBullets(glm::mat4 m, float color[3])
{
	printf("inside bullet function\n");
	glLoadMatrixf(&m[0][0]); 
	glColor3f(color[0],color[1],color[2]); 
	glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*)NULL+0); 

}
void draw_square(glm::mat4 m, float color[3]) 
{
  glLoadMatrixf(&m[0][0]); 
  glColor3f(color[0],color[1],color[2]); 
  glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*)NULL+0); 
} 

void draw_triangle(glm::mat4 m, float color[3]) 
{
	glLoadMatrixf(&m[0][0]); 
    glColor3f(color[0],color[1],color[2]);
    glDrawElements(GL_TRIANGLES, 36, GL_UNSIGNED_BYTE, (char*)NULL+0); 
} 

/////////////////////////////////////////////////////////////

vector <glm::mat4>::const_iterator mi;
  float targety = -0.4;
  float targetx =  0.9;
  float color[3]; 
  bool col = true;
  float anglex;
  float colort[3];
  float bullety = 0.0f;
  float xTargetRotate = 0.0f;
  

void displayMainBox()
{
	if( won == true)
			{
				glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT); 
				glClearColor(.18,0,.18,1); 
				char menu[] = "YOU WON !!!!";
				int len = strlen(menu);
				glColor3f(0.5f, 0.5f, 0.2f);
				glRasterPos2f(-1.0f, -0.2f);
				for (int p = 0; p < len; p++) 
					glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_24, menu[p]);
				printf("you won!");
				glutSwapBuffers();
				getchar();
				exit(1);
			}
	glClear(GL_COLOR_BUFFER_BIT| GL_DEPTH_BUFFER_BIT); 
    glClearColor(.18,0,.18,1); 
	glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array 
    	
			
	glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]); 
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO); 

  
  glVertexPointer(3,GL_FLOAT, 0,0);
  modeltar =  glm::mat4(1.0f);
  modeltars=  glm::mat4(1.0f);
  //////////////////////////////////////////
  //////////////////     target   /////////
  ////////////////////////////////////////
  colort[0] = 1.0; colort[1] = 0.47; colort[2] = 0.0;
  modeltar=glm::scale(modeltar, glm::vec3(1.0f-angle1, 1.0f-angle2, 1.0f));				
  modeltar=glm::translate(modeltar, glm::vec3(targetx, targety, 1.0f));		
  tar_stack.push(modeltar);
  modeltar=glm::rotate(modeltar,xTargetRotate,glm::vec3(0.0f, 0.0f, 1.0f));
  //modeltar=glm::rotate(modeltar,xTargetRotate,glm::vec3(1.0f, 0.0f, 0.0f));
  modeltar=glm::rotate(modeltar,xTargetRotate,glm::vec3(0.0f, 1.0f, 0.0f));
  draw_triangle(modeltar, colort);
  

  //////////////////////////////////////////
  ////////////////// inner target /////////
  ////////////////////////////////////////
  modeltar = tar_stack.top();
  tar_stack.pop();
  float scalext  = 0.6f-angle1, scaleyt = 0.6f-angle2;
  if (0.6f-angle1 < 0) scalext = 0;
  if (0.6f-angle2 < 0) scaleyt = 0;
  colort[0] = .80; colort[1] = 0.27; colort[2] = 0.0;
  //modeltar = glm::translate(modeltar, glm::vec3(targetx, targety, 1.0f));					//inner target
  modeltar = glm::scale(modeltar, glm::vec3(scalext, scaleyt, 1.0f));				
  modeltar = glm::rotate(modeltar,-xTargetRotate,glm::vec3(0.0f, 0.0f, 1.0f));
  //modeltar = glm::rotate(modeltar,-xTargetRotate,glm::vec3(1.0f, 0.0f, 0.0f));
  modeltar = glm::rotate(modeltar,-xTargetRotate,glm::vec3(0.0f, 1.0f, 0.0f));
  draw_triangle(modeltar, colort);

  color[0] = 0; color[1] = 0; color[2] = 1;
  
  
  modelm = modelM;	//copy matrix position for bullet

  //////////////////////////////////////////////
  ///////////////   Nozzle    /////////////////
  ////////////////////////////////////////////
 
  mat_stack.push(modelM);
  modelM = glm::translate(modelM, glm::vec3(-.55f, 0.0f, 0.0f));
  
  modelM = glm::scale(modelM, glm::vec3(3.5f, 0.6f, 1.0f));
  modelnozzle = modelM;

  
  
 
 modelM = mat_stack.top();  mat_stack.pop(); 
 mat_stack.push(modelM);

  //////////////////////////////
  //    Big Gun Handle  ///////
  /////////////////////////////
  
  color[0] = 1; color[1] = 0; color[2] = 0;				
  modelM = glm::translate(modelM, glm::vec3(-.75f, -0.15f, 0.0f));
  modelM = glm::rotate(modelM, angle1, glm::vec3(0.0f, 0.0f, 1.0f));
  modelM = glm::shearX3D(modelM,0.4f,0.5f);
  
  mat_stack.push(modelM);  
  modelM = glm::scale(modelM, glm::vec3(0.4f, 1.5f, 1.0f)); 
  draw_square(modelM, color);
  modelM = mat_stack.top();  mat_stack.pop(); 
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
  
  //////////////////////////////
  //    Small Gun Handle //////
  /////////////////////////////
  color[0] = 1; color[1] = 0; color[2] = 0;								
  modelM = glm::translate(modelM, glm::vec3(-0.55f, -0.1f, 0.0f));
  modelM = glm::rotate(modelM, -angle2, glm::vec3(0.0f, 0.0f, -1.0f));
  
  mat_stack.push(modelM);  
  modelM = glm::scale(modelM, glm::vec3(0.4f, 1.0f, 1.0f)); 
  draw_square(modelM, color);
  modelM = mat_stack.top();  mat_stack.pop(); 
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
  color[0] = 0; color[1] = 1; color[2] = 1;
  draw_square(modelnozzle,color);		//main nozzle
  ///////////////////////////////
  //////// Extra parts//////////
  /////////////////////////////
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
  color[0] = 0.65; color[1] = 0; color[2] = 0;				
  modelM = glm::translate(modelM, glm::vec3(-.55f, 0.0f, 0.0f));
  modelM = glm::scale(modelM, glm::vec3(2.5f, 0.3f, 1.0f));
  draw_square(modelM,color);			//pipe inside nozzle
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
  
  modelM = glm::translate(modelM, glm::vec3(-.85f, 0.0f, 0.0f));
  modelM = glm::scale(modelM, glm::vec3(0.5f, 0.3f, 0.5f));
  
  draw_square(modelM,color);		//back stub of gun
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
  
  modelM = glm::translate(modelM, glm::vec3(-.25f, 0.0f, 0.0f));
  modelM = glm::scale(modelM, glm::vec3(0.9f, 0.4f, 0.5f));
  color[0] = 0; color[1] = 1; color[2] = 1;
  draw_square(modelM,color);		//silencer, small box in front of nozzle
  modelM = mat_stack.top();  mat_stack.pop(); 
  mat_stack.push(modelM);
}
void display() 
{ 
  
if(bullet == true)
  { 
	anglex = 0;
	for (float i = 0 ; i<=2 ; i+= 0.05)
	{   
			if(i>targetx && abs(bullety-targety) <= 0.15 )
					{glutSwapBuffers();  break;	}
		displayMainBox();  
		modelm = glm::scale(modelm, glm::vec3(0.6f, 0.3f, 1.0f));	
		modelm = glm::translate(modelm,glm::vec3(-.25f, 0.0f, 0.0f));
		anglex += 1;						//used to rotate the bullet
		color[0] = 1; color[1] = 1; color[2] = 1;		
		printf("inside bullet\n");
		modelm = glm::translate(modelm, glm::vec3(i, 0.0f , 0.0f));
		 modelm = glm::rotate(modelm, anglex, glm::vec3(1.0f, -0.2f, 0.0f));
		 color[0] = 1; color[1] = 1; color[2] = 0;
		
		glEnableClientState(GL_VERTEX_ARRAY);
		
			drawBullets(modelm,color);
		Sleep(50+gun_angle*4);
		glDisableClientState(GL_VERTEX_ARRAY); 
		//glDisableClientState(GL_COLOR_ARRAY);
		glutSwapBuffers(); 
  
	} // end of for loop
	bullet = false;
	if(abs(bullety-(targety)) <= 0.15)// || (bullety - ceil(targety*10)/10)< 0.1 )
		{   
			won = true; 
			
		}
	else printf("bullety:%f, targety:%f",bullety,ceil(targety*10)/10,floor(targety*10)/10);
		printf(" ****\n"); 
		for (mi=mat_list.begin(); mi!=mat_list.end(); mi++) 
			{
			printf(" hello bullet!\n"); 
			//draw_square((*mi), color); 
			} 
  
		glDisableClientState(GL_VERTEX_ARRAY); 
		glutSwapBuffers(); 
	}	// end of if 
	else 
	{
		displayMainBox();
		//printf(" **--**\n"); 
		for (mi=mat_list.begin(); mi!=mat_list.end(); mi++) 
			{
			printf(" hello non-bullet!\n"); 
			//draw_square((*mi), color); 
			} 
  
		glDisableClientState(GL_VERTEX_ARRAY); 
		glutSwapBuffers();
	}
}


///////////////////////////////////////////////////////////////////////////////
//////////////////////////// decide which key to press  //////////////////////
/////////////////////////////////////////////////////////////////////////////
void mykey(unsigned char key, int x, int y)
{
	float d_angle = 0.05; 

	if (key == 'q') exit(1); 
	if (key == ' ')
	{
								//bullets
		bullet = true;
    
	}
	if (key == 'a') 
		{
			modelM = glm::rotate(modelM, d_angle, glm::vec3(1.0f, 1.0f, 0.0f)); 
			gun_angle += d_angle;
			if(gun_angle>= 6.2) gun_angle = 0;
			printf("gun_angle %f\n",gun_angle);
		}
	if (key == 'w') 
	  { modelM = glm::translate(modelM, glm::vec3(0.0f, 0.1f, 0.0f)); 
		bullety += 0.1f;
		printf("inside w bullety:%f\n",bullety);
		}
	//modelM =  translate(modelM, 0,.1,0);
	if (key == 's') 
	  {modelM = glm::translate(modelM, glm::vec3(0.0f, -0.1f, 0.0f)); 
		bullety -= 0.1f;
	}
	//modelM =  translate(modelM, 0,-.1,0);
	if (key == 'c') {
	  modelM =  glm::mat4(1.0f);
	  angle1 = angle2 = angle3 = angle4 = 0; 
	}

	if (key == '1') {
          angle1 += 0.05;
          printf(" hello!\n"); 
        }
        if (key == '2') 
          angle2 += 0.05;

	
	if (key == 'p')  {
	  glm::mat4 pm = glm::scale(modelM, glm::vec3(0.5f, 0.5f, 1.0f)); 
	  mat_list.push_back(pm); 
	} 

	glutPostRedisplay(); 
}


void animation(void)
{
	
 timed = time(NULL);
 if(timed-times >= 5)
	 {
		 targety = ((double) rand() / (RAND_MAX));
		 targetx = (((float) rand() / RAND_MAX) * 0.5) + 0.5;
		 times=time(NULL);
	 }
	 xTargetRotate += 0.005;
	 display(); 
}

int main(int argc, char** argv) 
{ 
	 
times=time(NULL);
  glutInit(&argc, argv); 
  glutInitDisplayMode(GLUT_RGB|GLUT_DOUBLE|GLUT_DEPTH); 
  glutInitWindowSize(600,600); 

  glutCreateWindow("Shooting Game"); 
  glutDisplayFunc(display); 
  glutKeyboardFunc(mykey); 
  glutIdleFunc(animation);
  //glutReshapeFunc(reshape);
  mat_list.clear(); 
  
#ifdef __APPLE__
#else
  glewInit(); 
#endif


  InitVBO();
  	glutMainLoop(); 
  
  printf("out of main loop");
  getchar();
} 
