#version 430


out vec4 FragColor;
varying vec2 f_texpos;
uniform sampler2D texture1;

void main(void) {
    FragColor=texture2D(texture1,f_texpos);
}

