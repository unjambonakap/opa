#include <opa/math/game/mesh_primitives.h>

OPA_NM_MATH_GAME

PointVec compute_sphere_points(double radius, double z, int n) {
  PointVec res(n);
  double angle = 2 * PI / n;
  double r = sqrt(radius * radius - z * z);
  REP (i, n) {
    double xg1 = cos(i * angle);
    double yg1 = sin(i * angle);
    // printf("got %f %f %f %f\n", xg1, yg1, xg2, yg2);

    res[i] = glm::vec3(xg1 * r, yg1 * r, z);
  }
  return res;
}

void triangularize_sphere_strip(FaceCollection *fc, double radius, int ntr,
                                double zl, double zh) {

  double angle = 2 * PI / (ntr / 2);
  double rl = sqrt(radius * radius - zl * zl);
  double rh = sqrt(radius * radius - zh * zh);
  PointVec ll, lh;
  ll = compute_sphere_points(radius, zl, ntr / 2);
  lh = compute_sphere_points(radius, zh, ntr / 2);
  ll.push_back(ll[0]);
  lh.push_back(lh[0]);
  REP (i, ntr / 2) {
    fc->push(Triangle{ ll[i], lh[i], ll[i + 1] });
    fc->push(Triangle{ lh[i], ll[i + 1], lh[i + 1] });
  }
}

void triangularize_sphere_fan(FaceCollection *fc, double radius, int ntr,
                              double zfix, double zfan) {
  double angle = 2 * PI / (ntr / 2);
  PointVec tb = compute_sphere_points(radius, zfan, ntr);
  Pos pfix(0, 0, zfix);
  tb.push_back(tb[0]);
  REP (i, ntr) { fc->push(Triangle{ pfix, tb[i], tb[i + 1] }); }
}

void triangularize_sphere(FaceCollection *fc, double radius, int ntr) {
  int nstrips = sqrt(ntr);
  int ntr_per_strip = nstrips + 1 & ~1; // make it even

  double r2 = radius * radius;
  double area = 4 * PI * r2;

  double curz = radius;
  double area_per_strip = area / nstrips;
  double coeff = 2 * PI * radius;

  if (1) {
    double mul = 2 * radius / nstrips;
    FOR (i, 1, nstrips - 1) {
      triangularize_sphere_strip(fc, radius, ntr_per_strip, i * mul - radius,
                                 (i + 1) * mul - radius);
    }
    triangularize_sphere_fan(fc, radius, ntr_per_strip, -radius, -radius + mul);
    triangularize_sphere_fan(fc, radius, ntr_per_strip, radius, radius - mul);
  }
}

Pos get_octahedron_point(int id) { return vec_tb[id / 3] * OPA_BITSIGN(id); }

void add_octahedron(FaceCollection *fc) {
  REP (i, 8)
    fc->push(Triangle{ get_octahedron_point(i & 1),
                       get_octahedron_point(GETB(i, 1) + 2),
                       get_octahedron_point(GETB(i, 2) + 4) });
}

PointVec icosahedron_vertices() {
  double phi = (1 + sqrt(5)) / 2;
  PointVec tb;
  REP (i, 4)
    REP (j, 3) {
      Pos cur;
      cur[j] = 0;
      cur[(j + 1) % 3] = OPA_BITSIGN(GETB(i, 0)) * 1;
      cur[(j + 2) % 3] = OPA_BITSIGN(GETB(i, 1)) * phi;
      tb.push_back(cur);
    }
  return tb;
}

void add_dodecahedron(FaceCollection *fc) { OPA_CHECK0(false); }
void add_icosahedron(FaceCollection *fc) {
  Mesh mesh;
  MeshBuilder::PolyhedraFromCloud(&mesh, icosahedron_vertices());
  fc->push(*mesh.to_faces());
}

OPA_NM_MATH_GAME_END
