#version 130


attribute vec3 pos;
attribute vec3 x_color;

uniform mat4 view;
varying vec3 f_color;


void main(void) {
    vec4 tmp = vec4(pos, 1.0);
    f_color=x_color;
    gl_Position=view*tmp;
}


