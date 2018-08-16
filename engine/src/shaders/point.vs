#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec3 color;

//uniform vec3 cam_left;
//uniform vec3 cam_up;
uniform mat4 view;

varying vec4 g_color;

void main(void) {
  //vec2 tmp = (sq_pos - vec2(0.5, 0.5)) * pos.w;
  //vec3 p_pos = pos.xyz + cam_left * tmp.x + cam_up * tmp.y;

  //f_color = color;
  //gl_Position = view * ;

  g_color = vec4(color, 0.3);
  gl_Position = view * vec4(pos, 1.0);
}
