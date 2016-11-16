#ifdef __APPLE__
#include <GLUT/glut.h> 
#include <OpenGL/gl.h> 
#else
#include <GL/glew.h> 
#include <GL/glut.h> 
#include <GL/gl.h> 
#endif

#include<glm/glm.hpp>
#include <glm/gtx/transform.hpp>
#include <stdio.h>
#include <stdlib.h>
#include<vector>

int width, height, w,h, border=30, flagv1, flagv2, ratio, sp;;
std::vector<float> v1, v2, v3, v4, v5, vtempj, vtempk ;
int oneboxw, oneboxh,pts,size,ptGlobal=0;
float max_v1 = -10000, min_v1 = 10000;
float max_v2 = -10000, min_v2 = 10000;
float max_v3 = -10000, min_v3 = 10000;
float max_v4 = -10000, min_v4 = 10000;
float min_vk = -10000, min_vj = -10000;
float *vertices; 
GLubyte *pindices; 

GLuint vboHandle[1]; 
GLuint indexVBO; 

float scale_v1, scale_v2, scale_v3, scale_v4, scale_vk, scale_vj; 



void InitVBO() 
{
  glGenBuffers(1, &vboHandle[0]);   // create a VBO handle
  glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]);   // bind the handle to the current VBO 
  glBufferData(GL_ARRAY_BUFFER, size*2*4, vertices, GL_STATIC_DRAW); // allocate space and copy the data over
  glBindBuffer(GL_ARRAY_BUFFER, 0);   // clean up 

  glGenBuffers(1, &indexVBO); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO); 
  glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLubyte)*size, pindices, GL_STATIC_DRAW); 
  glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);  //clean up 
  
} 

void display()
{
	 int i, j, k,t;
	 size = v1.size();
	 printf("size = %d\n", size);
	 glClearColor(0.6, 0.16, 0.16, 1); 
	 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
	 float oneboxw=(width-5*border)/4.0;		//height and width of smaller box
     float oneboxh= (height-5*border)/4.0;
	 vertices = new float[size*2]; 
	 pindices = new GLubyte[size]; 
    
	

  for( j = 0 ; j < 4 ; j++ )			//j and k loops are for the rows and columns controlling 4 boxes in each row and each column
	{   
		switch (j)				// assigns to temp which variable will be plotted along j axis
			{ 
			case 0:
				vtempj = v1;
				scale_vj = scale_v1;
				min_vj = min_v1;
				break;
			case 1:
				vtempj = v2;
				scale_vj = scale_v2;
				min_vj = min_v2;
				break;
			case 2:
				vtempj = v3;
				scale_vj = scale_v3;
				min_vj = min_v3;
				break;
			case 3:
				vtempj = v4;
				scale_vj = scale_v4;
				min_vj = min_v4;
				break;
			default:
				break;
			}


		for ( k = 0;k <4 ; k++)			//j and k loops are for the rows and columns controlling 4 boxes in each row and each column
		{ 
			int w = (k+1)*border + k*oneboxw;		//drawing bondary boxes and inner black coloured box using glscissor
			int h = (j+1)*border + j*oneboxh;
			glScissor(w-2, h-2, oneboxw+4, oneboxh+4);
			glEnable(GL_SCISSOR_TEST);
			glClearColor(0,0.9,0.9,1);
			glClear(GL_COLOR_BUFFER_BIT);
			
			glScissor(w, h, oneboxw,oneboxh);
			glEnable(GL_SCISSOR_TEST);
			glClearColor(0.1,0.1,0.1,1);
			glClear(GL_COLOR_BUFFER_BIT);
			glDisable(GL_SCISSOR_TEST);
			 
			glFlush();
			
		glViewport(w, h-border, oneboxw, border);	//label x-axis
			char menu[20];
			if(j == 0)
				strcpy(menu,"sepal-length");
			else if (j == 1)
				strcpy(menu,"sepal-width");
			else if (j == 2)
				strcpy(menu,"petal-length");
			else strcpy(menu,"petal-width");

			glColor3f(0.0f, 1.0f, 0.0f);
			glRasterPos2f(-0.4, -0.3);
		    int len = strlen(menu);
			for (int p = 0; p < len; p++) 
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, menu[p]);
			
			glViewport(w-border, h, border, oneboxh);		//label y-axis
			
			if(k == 0)
				strcpy(menu,"sl");
			else if (k == 1)
				strcpy(menu,"sw");
			else if (k == 2)
				strcpy(menu,"pl");
			else strcpy(menu,"pw");

			glColor3f(0.0f, 1.0f, 0.0f);
			glRasterPos2f(-0.3, -0.1);
		    len = strlen(menu);
			for (int p = 0; p < len; p++) 
			glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, menu[p]);
			
			glFlush();
			
			glViewport(w, h, oneboxw, oneboxh);	//draw horizontal and vertical tick-mark lines
			
			glEnable(GL_ALPHA_TEST);				//for making round points
			//glAlphaFunc(GL_NOTEQUAL, 0);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable( GL_POINT_SMOOTH );
			
			glBegin(GL_LINES);
			glEnable(GL_BLEND);
			glColor4f(0, 0.9, 0.9,0.5);
			for (int m=0;m<11;m++)		//ticks on each viewport
			{
				
				glVertex2f((2*(float)m/10)-1, -1); 
				glVertex2f((2*(float)m/10)-1, +1);
				
				glVertex2f(-1,(2*(float)m/10)-1); 
				glVertex2f(1,(2*(float)m/10)-1);
				
			}
			glDisable(GL_BLEND);
			glDisable(GL_POINT_SMOOTH);
			glBlendFunc(GL_NONE, GL_NONE);
			glDisable(GL_BLEND);
			glEnd();
			glFlush();

			switch (k)				// assigns to temp which variable will be plotted along k axis
			{ 
			case 0:
				vtempk = v1;
				scale_vk = scale_v1;
				min_vk = min_v1;
				break;
			case 1:
				vtempk = v2;
				scale_vk = scale_v2;
				min_vk = min_v2;
				break;
			case 2:
				vtempk = v3;
				scale_vk = scale_v3;
				min_vk = min_v3;
				break;
			case 3:
				vtempk = v4;
				scale_vk = scale_v4;
				min_vk = min_v4;
				break;
			default:
				break;
			}

		for(int g=0;g<3;g++)		// for three species of flowers in dataset
			{
			for ( i=g*(size/3), t=0; i<g*(size/3)+(size/3); t++,i++)	//taking each species of flowers, divide by 3 as there are three species
			{   
			

				float x = (vtempk[i]-min_vk)*2/scale_vk - 1 ;
				float y = (vtempj[i]-min_vj)*2/scale_vj - 1 ;
				vertices[2*t] = x;		//store to pass in VBO
				vertices[2*t+1] = y; 
				pindices[t] = t; 
			
			}	//end of i loop
			pts = ptGlobal;
			if(g == 0)			//choose color and size for each species
				{ glColor3f(1, 0, 0); pts=pts+(g*2);}		//size of each species is a function of the input size
			else if (g == 1)
					{glColor3f(0, 1, 0);	 pts=pts+(g*2);}
			else 	{glColor3f(0, 0, 1);	 pts=pts+(g*2);}
			glPointSize(pts);
			#ifdef __APPLE__
			#else
				glewInit(); 
			#endif
			glEnable(GL_ALPHA_TEST);				//for making round points
			glAlphaFunc(GL_NOTEQUAL, 0);
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
			glEnable( GL_POINT_SMOOTH );
			
			InitVBO();							  //draw points using VBO
			glEnableClientState(GL_VERTEX_ARRAY); // enable the vertex array 
			glBindBuffer(GL_ARRAY_BUFFER, vboHandle[0]); 
			glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, indexVBO);
			glVertexPointer(2,GL_FLOAT, 0,0);
			glViewport(w, h, oneboxw, oneboxh);
			glDrawElements(GL_POINTS, size, GL_UNSIGNED_BYTE, (char*)NULL+0); 
			glDisableClientState(GL_VERTEX_ARRAY); 
			glDisable(GL_POINT_SMOOTH);
			glBlendFunc(GL_NONE, GL_NONE);
			glDisable(GL_BLEND);
			
	
			glFlush();
			}//end of g
    
		}//end of k
	}//end of j
	
}	//end of function

void handleResize(GLint w,GLint h)
{
	glMatrixMode(GL_PROJECTION);	// Select The Projection Matrix
    glLoadIdentity();		// Reset The Projection Matrix
	glClearColor(0.6, 0.16, 0.16, 1); 
	glClear(GL_COLOR_BUFFER_BIT);
	height =h;
	width= w;
}

int main(int argc, char** argv) 
{
	 float f1, f2, f3, f4;
	 float *d1, *d2, *d3, *d4;
	 int border = 20, cnt = 0, i;;
	 char h1[20],h2[20],h3[20],h4[20],h5[20],type[20],species[20];
	 width = atoi(argv[1]);
	height = atoi(argv[2]);
	   ptGlobal = atoi(argv[4]);
  printf(" width %d  height %d plot-point size %d file %s \n", width, height, pts, argv[3]); 
  
  FILE* fp = fopen(argv[3], "r");
  
   
  fscanf(fp, "%20[^,\n],%20[^,\n],%20[^,\n],%20[^,\n],%20[^,\n]\n", h1, h2, h3, h4, h5);		//read file header
  
  printf("File Headers:\n h1: %s \n h2: %s \n h3: %s \n h4: %s \n h5: %s\n",h1,h2,h3,h4,h5);

  while (fscanf(fp, "%19[^ ,],%g,%g,%g,%g \n", species, &f1, &f2, &f3, &f4)==5) 
  {
	
	if (strcmp(species,"setosa")==0)
		sp  = 1;
	if (strcmp(species,"versicolor")==0)
		sp = 2;
	if (strcmp(species,"virginica")==0)
		{sp = 3;
		}
	
    if (f1 < min_v1) min_v1 = f1;
    if (f1 > max_v1) max_v1 = f1;
    
    if (f2 < min_v2) min_v2 = f2;
    if (f2 > max_v2) max_v2 = f2;

    if (f3 < min_v3) min_v3 = f3;
    if (f3 > max_v3) max_v3 = f3;


    if (f4 < min_v4) min_v4 = f4;
    if (f4 > max_v4) max_v4 = f4;

    
    v1.push_back(f1);
    v2.push_back(f2);
    v3.push_back(f3);
    v4.push_back(f4); 
	v5.push_back(sp);
    cnt++; 
  }
  scale_v1 = max_v1 - min_v1;
  scale_v2 = max_v2 - min_v2;
  scale_v3 = max_v3 - min_v3;
  scale_v4 = max_v4 - min_v4; 
 
  glutInit(&argc, argv);
  glutInitDisplayMode(GLUT_RGB|GLUT_SINGLE);
  glutInitWindowSize(width,height);
  glutCreateWindow("Scatter plot");
  glClearColor(0.6, 0.16, 0.16, 1); 

  glutDisplayFunc(display);
  glutReshapeFunc(handleResize);

  glutMainLoop();
}
