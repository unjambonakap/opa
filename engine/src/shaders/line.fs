#version 430 core

out vec4 FragColor;
varying vec3 f_color;

void main(void) { FragColor = vec4(f_color, 1); }
