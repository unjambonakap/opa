#version 430 core

layout(location = 0) in vec2 sq_pos;
layout(location = 1) in vec4 pos1;
layout(location = 2) in vec4 pos2;
layout(location = 3) in vec3 color;

uniform vec3 cam_left;
uniform vec3 cam_up;
uniform mat4 view;

varying vec3 f_color;

void main(void) {
  vec3 dir = pos2.xyz - pos1.xyz;
  vec3 orth = cross(cam_left, cam_up);
  vec3 orth_dir = normalize(cross(dir, orth));
  vec3 p_pos = pos1.xyz + dir * sq_pos.x + orth_dir * sq_pos.y * 0.1;

  f_color = color;
  gl_Position = view * vec4(p_pos, 1.0);
}
