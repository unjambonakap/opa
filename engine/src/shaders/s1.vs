#version 430 core

layout(location = 0) in vec3 pos;
layout(location = 1) in vec2 texpos;

varying vec2 f_texpos;

uniform mat4 view;

void main(void) {
    vec4 tmp = vec4(pos, 1.0);
    f_texpos = texpos;
    gl_Position = view * tmp;
}
