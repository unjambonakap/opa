#pragma once

#include <opa/math/game/base.h>
#include <opa/math/game/mesh_util.h>

OPA_NM_MATH_GAME

void triangularize_sphere(FaceCollection *fc, double radius, int ntr);

void add_octahedron(FaceCollection *fc);
void add_dodecahedron(FaceCollection *fc);
void add_icosahedron(FaceCollection *fc);

OPA_NM_MATH_GAME_END
