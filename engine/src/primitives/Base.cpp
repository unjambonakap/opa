#include "primitives/Base.h"
#include "Drawer.h"

using namespace opa::math::game;
using namespace opa::utils;
OPA_NAMESPACE_DECL2(opa, engine)

Square::Square() {}

void TrianglesSceneObj::init(Game *game) {
  SceneObj::init(game);
  auto x = new TriangleContainerIsec;
  x->init(TriangleContainerIsecParams(&this->tr()));
  auto y = new TriangleContainerDrawer;
  y->init(TriangleContainerDrawer::Params(&this->tr()));

  isec().register_isec(IsecPtr(x));
  drawer().register_drawer(DrawerPtr(y));

}

void Square::init(Game *game, float side, TexturePtr texture,
                  const std::vector<TexUV> &texpos) {
  SceneObj::init(game);
  auto x = new TriangleContainerIsec;
  x->init(TriangleContainerIsecParams(&this->tr()));
  auto y = new TriangleContainerDrawer;
  y->init(TriangleContainerDrawer::Params(&this->tr()));

  isec().register_isec(IsecPtr(x));
  drawer().register_drawer(DrawerPtr(y));

  std::vector<Pos2> coords;
  coords.pb(Pos2(0, 0));
  coords.pb(Pos2(1, 0));
  coords.pb(Pos2(0, 1));
  coords.pb(Pos2(1, 1));

  REP (i, 2) {
    TexturedTriangle &tr = this->tr().triangles().get_new2();
    tr.init(texture);
    REP (j, 3) {
      if (texpos.size() == 0) {
        tr.data()[j].texpos.x = coords[i + j].x;
        tr.data()[j].texpos.y = 1 - coords[i + j].y;
      } else
        tr.data()[j].texpos = texpos[i + j];
      tr.data()[j].pos =
        glm::vec3(0, coords[j + i] - glm::vec2(0.5, 0.5)) * side;
    }
    tr.refresh();
  }
}
void Cube::init(Game *game, float side, TexturePtr *textures) {
  SceneObj::init(game);

  REP (i, 6) {
    Square *square = new Square;
    square->init(game, side, textures[i], {});
    square->status().weak = 1;
    square->fixed() = true;

    float v = OPA_BITSIGN(i & 1) * side / 2;
    if (i & 4) {
      square->pos_obj().pos().z = v;
      square->pos_obj().rot() = quat_look_at(vec_z, vec_x);
    } else if (i & 2) {
      square->pos_obj().pos().y = v;
      square->pos_obj().rot() = quat_look_at(vec_y, vec_x);
    } else {
      square->pos_obj().pos().x = v;
      square->pos_obj().rot() = quat_look_at(vec_x, vec_y);
    }
    m_game->register_obj(SceneObjPtr(square), this);
  }
}

void Icosahedron::init(Game *game, float side, TexturePtr *textures) {
  SceneObj::init(game);
  auto x = new TriangleContainerIsec;
  x->init(TriangleContainerIsecParams(&this->tr()));
  auto y = new TriangleContainerDrawer;
  y->init(TriangleContainerDrawer::Params(&this->tr()));

  isec().register_isec(IsecPtr(x));
  drawer().register_drawer(DrawerPtr(y));

  Pos pl, ph;
  pl = -vec_z * side;
  ph = -pl;

  float step = 2 * PI / 5;
  int ns = 5;
  float ang2 = PI / 6;
  std::vector<Pos> vl(ns), vh(ns);
  REP (i, ns) {
    vh[i] = glm::vec3(cos(i * step), sin(i * step), sin(ang2));
    vh[i].xy() *= cos(ang2);
    vh[i] *= side;

    vl[i] =
      glm::vec3(cos(i * step + step / 2), sin(i * step + step / 2), -sin(ang2));
    vl[i].xy() *= cos(ang2);
    vl[i] *= side;
  }

  std::vector<TexturedTriangle *> trs;
  std::vector<IdType> ids = this->tr().triangles().get_new_vec(20);
  REP (i, ids.size()) {
    TexturedTriangle *tr = &this->tr().triangles().get(ids[i]);
    trs.pb(tr);
    tr->init(textures[i]);
    tr->data()[0].texpos = Pos2(0, 0);
    tr->data()[1].texpos = Pos2(1, 0);
    tr->data()[2].texpos = Pos2(0, 1);
  }
  vl.pb(vl[0]);
  vh.pb(vh[0]);

  REP (i, 5) {
    trs[i]->data()[0].pos = vl[i];
    trs[i]->data()[1].pos = vh[i + 1];
    trs[i]->data()[2].pos = vl[i + 1];

    trs[5 + i]->data()[0].pos = vh[i];
    trs[5 + i]->data()[1].pos = vl[i];
    trs[5 + i]->data()[2].pos = vh[i + 1];

    trs[10 + i]->data()[0].pos = vh[i];
    trs[10 + i]->data()[1].pos = ph;
    trs[10 + i]->data()[2].pos = vh[i + 1];

    trs[15 + i]->data()[0].pos = vl[i];
    trs[15 + i]->data()[1].pos = pl;
    trs[15 + i]->data()[2].pos = vl[i + 1];
  }

  REP (i, 20)
    trs[i]->refresh();
}

OPA_NAMESPACE_DECL2_END
