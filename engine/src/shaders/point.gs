#version 440 core
#define M_PI 3.1415926535897932384626433832795
const int N_TRS = 8;

layout (points) in;
layout (triangle_strip, max_vertices = N_TRS*3) out;

in vec4 g_color[];
out vec4 f_color;
uniform float r_fact;

void draw_point(vec4 pos){
  int i;
  float step = 2 * M_PI / N_TRS;
  float r = r_fact * (2 - pos.z / pos.w);
  for (i = 0; i < N_TRS; i++) {
    gl_Position = pos + vec4(r*cos(i*step), r*sin(i*step), 0, 0);
    EmitVertex();

    gl_Position = pos + vec4(r*cos((i+1)*step), r*sin((i+1)*step), 0, 0);
    EmitVertex();

    gl_Position = pos;
    EmitVertex();
    EndPrimitive ();
  }
}

void main() {
  f_color = g_color[0];
  draw_point(gl_in[0].gl_Position);
}
