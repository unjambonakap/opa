#include "Sphere.h"
#include "procgen/Texture.h"

OPA_NAMESPACE_DECL2(opa, engine)

Sphere::Sphere() {}

static const float eps = 1e-9;

void Sphere::init(Game *game, float radius) {
  m_mesh_cache.init(3);

  SceneObj::init(game);
  auto x = new TriangleContainerDrawer;
  x->init(TriangleContainerDrawer::Params(&tr()));
  drawer().register_drawer(DrawerPtr(x));

  auto y = new SphereIsec;
  y->init(SphereIsecParams(radius));
  isec().register_isec(IsecPtr(y));

  m_r = radius;
}

void Sphere::build_strip(TexturedStrip *strip, float zl, float zh, int ntr) {
  strip->init(ntr, m_texture);

  auto data = strip->data();
  float angle = 2 * PI / (ntr / 2);
  float rl = sqrt(m_r * m_r - zl * zl);
  float rh = sqrt(m_r * m_r - zh * zh);

  for (int i = 0; i < ntr / 2; ++i) {
    float xg1 = cos(i * angle);
    float yg1 = sin(i * angle);
    // printf("got %f %f %f %f\n", xg1, yg1, xg2, yg2);

    data[i * 2].pos = glm::vec3(xg1 * rl, yg1 * rl, zl);
    data[i * 2 + 1].pos = glm::vec3(xg1 * rh, yg1 * rh, zh);

    REP (j, 2)
      data[i * 2 + j].texpos = glm::vec2(0.f, 0.f);
  }
  data[ntr] = data[0];
  data[ntr + 1] = data[1];
  strip->refresh();
}

void Sphere::build(u32 npoly) {
  check_init();

  int nstrips = sqrt(npoly);
  int ntr_per_strip = nstrips + 1 & ~1;
  ntr_per_strip *= 4;

  t_buf = new u8[128 * 128 * 3];
  UniformTexturePtr build_tex(new UniformTexture());
  build_tex->init(
    UniformTexture::UniformTextureParameters(Color(0xff, 0xff, 0xff)));
  m_texture = build_tex;

  OPA_ASSERT0(nstrips > 5);
  nstrips = nstrips + 1 & ~1;

  float r2 = m_r * m_r;
  float area = 4 * PI * r2;

  float curz = m_r;
  float area_per_strip = area / nstrips;
  float coeff = 2 * PI * m_r;

  if (1) {
    float mul = 2 * m_r / nstrips;
    REP (i, nstrips) {
      auto &strip = tr().strips().get_new2();
      build_strip(&strip, i * mul - m_r, (i + 1) * mul - m_r, ntr_per_strip);
    }
  }
}

opa::math::game::Mesh *Sphere::build_mesh(double precision) const {
  // TODO: proper one
  return tr().to_tr_collection().to_mesh();
}

OPA_NAMESPACE_DECL2_END
