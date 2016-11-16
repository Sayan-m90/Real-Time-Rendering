uniform int use_texture;

attribute vec4 position; 
attribute vec4 color;
attribute vec4 normal; 
attribute vec2 texCoord;

varying vec3 N;
varying vec3 L;
varying vec3 R;
varying vec3 V;

varying vec2 FtexCoord;
varying vec4 pos_in_eye;


uniform mat4 local2clip;
uniform mat4 local2eye;
uniform mat4 normal_matrix;
uniform mat4 world2eye; 

uniform vec4 light_ambient;
uniform vec4 light_diffuse;
uniform vec4 light_specular;
uniform vec4 light_pos;

uniform vec4 mat_ambient;
uniform vec4 mat_diffuse;
uniform vec4 mat_specular;
uniform float mat_shine; 



void main(){


       gl_Position = local2clip * position;
	   vec4 Vpos =  local2eye * position; 
	   pos_in_eye = Vpos;
	   N = normalize(vec3(normal_matrix * normal)); 
	   vec4 Lpos =  world2eye * light_pos;
	   L = normalize(vec3(Lpos - Vpos)); 
	   R = normalize(reflect(-L, N)); 
	   V = normalize(vec3(-Lpos));
	   FtexCoord = texCoord;
		
}