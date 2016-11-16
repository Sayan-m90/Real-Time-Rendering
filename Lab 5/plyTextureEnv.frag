
uniform int use_texture;

varying vec2 FtexCoord;

uniform vec4 light_ambient;
uniform vec4 light_diffuse;
uniform vec4 light_specular;
uniform vec4 light_pos;

uniform vec4 mat_ambient;
uniform vec4 mat_diffuse;
uniform vec4 mat_specular;
uniform float mat_shine;

uniform mat4 local2clip;
uniform mat4 local2eye;
uniform mat4 normal_matrix;
uniform mat4 world2eye;
uniform mat4 eye2world;

varying vec3 N;
varying vec3 L;
varying vec3 R;
varying vec3 V;

varying vec4 pos_in_eye;

uniform sampler2D Tex1;
uniform sampler2D Tex2;
uniform sampler2D Tex3;
uniform sampler2D Tex4;
uniform sampler2D Tex5;
uniform sampler2D Tex6;
uniform sampler2D Tex7;
uniform samplerCube cubeMap;

void main()
{

	 vec4 ambient = light_ambient * mat_ambient;
	 vec4 Lightning;
	 vec4 Lightning2;

	 float NdotL; 
       if (dot(N,L) <0.0) NdotL = 0.0; 
       else NdotL = dot(N, L); 

	   vec4 diffuse = light_diffuse * mat_diffuse * NdotL;
  

      float RdotV; 
       RdotV = dot(R, V); 

       if (NdotL == 0.0) {RdotV = 0.0; }
       if (RdotV <0.0) {RdotV = 0.0;}
	  vec4 specular = light_specular * mat_specular * pow(RdotV,mat_shine);   
	  Lightning = diffuse+ambient+specular;
	  Lightning2 = diffuse+ambient+(2*specular);
	  vec4 texcolor;
	  vec4 env_color = vec4(1,0,0,1);
	  vec3 ref, view_vector;
	  if(use_texture == 0)
	  { 
	  view_vector = normalize(vec3(vec4(0,0,0,1)-pos_in_eye));
		ref = normalize(reflect(-view_vector,N));
		ref = vec3(eye2world*vec4(ref,0));
		env_color = textureCube(cubeMap, ref);
		gl_FragColor = env_color;
		gl_FragColor = env_color+(2*Lightning2);
	  }  
	  else if(use_texture ==1)
	  {

		texcolor = texture2D(Tex1,FtexCoord);
		gl_FragColor = texcolor;
	  }
	  else if(use_texture ==3)
	  {
	  texcolor = texture2D(Tex3,FtexCoord);
		gl_FragColor = texcolor;
	  }
	  else if(use_texture == 2)
	  {
		view_vector = normalize(vec3(vec4(0,0,0,1)-pos_in_eye));
		ref = normalize(reflect(-view_vector,N));
		ref = vec3(eye2world*vec4(ref,0));
		env_color = textureCube(cubeMap, ref);
		gl_FragColor = env_color;
	  }
	  else if(use_texture == 8)
	  {
	  view_vector = normalize(vec3(vec4(0,0,0,1)-pos_in_eye));
		ref = normalize(reflect(-view_vector,N));
		ref = vec3(eye2world*vec4(ref,0));
		env_color = textureCube(cubeMap, ref);
		gl_FragColor = env_color*Lightning2;
	  }
	   else if(use_texture == 4)
	  {
	  texcolor = texture2D(Tex4,FtexCoord);
		gl_FragColor = texcolor;
	  }
	   else if(use_texture == 5)
	  {
	  texcolor = texture2D(Tex5,FtexCoord);
		gl_FragColor = texcolor;
	  }
	   else if(use_texture == 6)
	  {
	  texcolor = texture2D(Tex6,FtexCoord);
		gl_FragColor = texcolor;
	  }
	   else if(use_texture == 7)
	  {
	  texcolor = texture2D(Tex7,FtexCoord);
		gl_FragColor = texcolor;
	  }
}





