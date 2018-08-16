#pragma once

#include <fplll/fplll.h>
#include <opa/algo/sat2.h>
#include <opa/math/co/base.h>
#include <opa/math/common/fast_gf2.h>
#include <opa/math/common/matrix_utils.h>
#include <opa/or/grid_search.h>

OPA_MATH_CO

class CombHelper {
public:
  CombHelper(int n = -1) { this->init(n); }
  void init(int n) {
    if (n == -1) return;
    fact.resize(n + 1);
    fact[0] = 1;
    REP (i, n)
      fact[i + 1] = fact[i] * (i + 1);
  }

  bignum cnk(int n, int k) const {
    if (n < k) return 0;
    if (k < 0) return 0;
    return fact[n] / fact[k] / fact[n - k];
  }

  bignum ank(int n, int k) const {
    if (n < k) return 0;
    if (k < 0) return 0;
    return fact[n] / fact[n - k];
  }

  std::vector<Z> fact;
};

class Hyperplane {
public:
  static constexpr bool kLinear = true;
  Hyperplane(const QModule_t &w, const Q &level) { this->init(w, level); }

  Hyperplane(const std::vector<QModule_t> &points, bool linear = !kLinear) {
    OPA_CHECK0(!points.empty());
    int n = points[0].size();
    Matrix<Q> tmat(&QF, n, n - 1 + n);
    REP (i, n)
      tmat(i, n - 1 + i) = QF.getE();

    if (!linear) {
      OPA_CHECK_EQ0(points.size(), n);
      FOR (i, 1, n)
        tmat.set_col(i - 1, (points[i] - points[0]));
    } else {
      OPA_CHECK_EQ0(points.size(), n - 1);
      tmat.set_cols(points);
    }

    tmat.row_echelon(n, n - 1);
    OPA_CHECK(tmat.get_submatrix(n - 1, 0, 1, n - 1).is_null(), tmat);

    QModule_t w(tmat.get_submatrix(n - 1, n - 1).tovec());
    Q level = QF.getZ();
    if (!linear) {
      level = w.dot(points[0]);
    }
    this->init(w, level);
  }

  bool contains_check_not_on(const QModule_t &pt) const {
    OPA_CHECK(!is_on(pt), pt);
    return contains(pt);
  }

  void toggle() {
    level = -level;
    w = -w;
  }

  bool reorient_to_contain(const QModule_t &pt) {
    Q tmp = pt.dot(w);
    OPA_CHECK(tmp != level, pt, w, level);
    if (tmp < level) {
      toggle();
      return true;
    }
    return false;
  }

  bool reorient_to_contain(const std::vector<QModule_t> &pts) {
    bool oriented = false;
    for (auto &pt : pts) {
      int status = this->contains_status(pt);
      if (status == 0) continue;
      if (status == -1) {
        if (oriented) return false;
        toggle();
      }
      oriented = true;
    }
    return true;
  }

  void init(const QModule_t &w, const Q &level) {
    this->w = w;
    this->level = level;
  }

  int contains_status(const QModule_t &x) const {
    Q v = x.dot(w);
    // OPA_DISP0(x, w, v, level, v == level, v > level);
    return v == level ? 0 : v > level ? 1 : -1;
  }

  bool contains(const QModule_t &x) const {
    Q v = x.dot(w);
    return v >= level;
  }

  bool is_on(const QModule_t &x) const { return x.dot(w) == level; }
  int dim() const { return w.size(); }
  OPA_DECL_COUT_OPERATOR2(Hyperplane, a.w, a.level);

  QModule_t w;
  Q level;
};

struct PolyhedraSimplicialDecomposition {
  std::vector<int> ray_ids;
  OPA_DECL_COUT_OPERATOR2(PolyhedraSimplicialDecomposition, a.ray_ids);
};

struct PolyhedraCone {
  QModule_t vertex;
  std::vector<ZModule_t> rays;
  std::vector<PolyhedraSimplicialDecomposition> simplicial_decomposition;
  OPA_DECL_COUT_OPERATOR2(PolyhedraCone, a.vertex, a.rays,
                          a.simplicial_decomposition);
};

struct Polytope_ConeRepr {
  std::vector<PolyhedraCone> cones;
  OPA_DECL_COUT_OPERATOR2(Polytope_ConeRepr, a.cones);
};

class Polyhedra_VertexRepr {
public:
  void init(const std::vector<QModule_t> &vertices) {
    this->vertices = vertices;
    n = vertices.size();
    d = vertices[0].size();
    OPA_DISP0(n, d, vertices);

    // Compute faces
    OR::SubsetSelector<int>(utils::range(0, n).tb(),
                            [&](const std::vector<int> &sel) {
                              this->maybe_add_repr(sel);
                              return OR::SubsetSelectorRes{ false, true };
                            },
                            d);

    // Compute cones
    std::vector<BitVec> vertex_to_incident_faces(n, BitVec(faces.size()));
    REP (fid, faces.size()) {
      for (auto &v : faces[fid]) vertex_to_incident_faces[v].set(fid, 1);
    }

    std::vector<std::vector<int> > vertex_to_ray_vertices(n);
    cone_repr.cones.resize(n);
    REP (i, n) {
      cone_repr.cones[i].vertex = vertices[i];
      FOR (j, i + 1, n) {
        BitVec nd =
          vertex_to_incident_faces[i].andz(vertex_to_incident_faces[j]);

        int count = nd.count_set();
        std::vector<QModule_t> incident_hps;
        REP (fid, faces.size()) {
          if (nd.get(fid)) {
            incident_hps.push_back(faces_hp[fid].w);
          }
        }
        Matrix<Q> incident_hps_mat(&QF, incident_hps.size(), d);
        incident_hps_mat.set_rows(incident_hps);
        int row_span = incident_hps_mat.row_echelon();
        if (row_span != d - 1) continue;

        ZModule_t ray = force_zmodule(vertices[j] - vertices[i]);
        vertex_to_ray_vertices[i].push_back(j);
        vertex_to_ray_vertices[j].push_back(i);
        cone_repr.cones[i].rays.push_back(ray);
        cone_repr.cones[j].rays.push_back(-ray);
      }
      OPA_CHECK(cone_repr.cones[i].rays.size() >= d, cone_repr.cones[i]);
    }

    REP (i, n) {
      std::vector<pii> ray_vertices;
      REP (j, vertex_to_ray_vertices[i].size())
        ray_vertices.emplace_back(vertex_to_ray_vertices[i][j], j);
      build_simplicial_decomposition(i, ray_vertices);
    }
  }

  void build_simplicial_decomposition(int vid,
                                      const std::vector<pii> &ray_vertices) {
    if (ray_vertices.size() == d) {
      std::vector<int> ray_ids;
      for (auto &ray_vertex : ray_vertices)
        ray_ids.push_back(ray_vertex.second);
      cone_repr.cones[vid].simplicial_decomposition.emplace_back(
        PolyhedraSimplicialDecomposition{ ray_ids });
      return;
    }
    // OPA_DISP("===== start decomp ", vid, ray_vertices.size());

    std::vector<pii> left, right;
    OR::SubsetSelector<pii>(
      ray_vertices,
      [&](const std::vector<pii> &split_face_vids) {
        std::vector<int> maybe_face_ids;
        maybe_face_ids.push_back(vid);
        for (auto &e : split_face_vids) maybe_face_ids.push_back(e.first);
        std::sort(ALL(maybe_face_ids));

        // These vertices form a face of the cone. can't
        // use it to split
        if (this->seen.count(maybe_face_ids))
          return OR::SubsetSelectorRes{ false, true };
        std::vector<QModule_t> splitv;
        for (auto &e : split_face_vids)
          splitv.push_back(lift_zmodule(cone_repr.cones[vid].rays[e.second]));
        Hyperplane split_hp(splitv, Hyperplane::kLinear);

        for (auto &e : ray_vertices) {
          int status = split_hp.contains_status(
            lift_zmodule(cone_repr.cones[vid].rays[e.second]));
          if (status >= 0) right.push_back(e);
          if (status <= 0) left.push_back(e);
        }
        if (left.size() < ray_vertices.size() &&
            right.size() < ray_vertices.size()) {
          return OR::SubsetSelectorRes{ true, true };
        }
        left.clear();
        right.clear();
        return OR::SubsetSelectorRes{ false, true };

      },
      d - 1);

    OPA_CHECK0(!left.empty() && !right.empty());
    this->build_simplicial_decomposition(vid, left);
    this->build_simplicial_decomposition(vid, right);
  }

  void maybe_add_repr(const std::vector<int> &vertex_sel) {
    std::vector<QModule_t> hp_vertices;
    for (const auto &face : faces) {
      if (utils::contained_sorted(vertex_sel, face)) return;
    }

    for (auto &v : vertex_sel) hp_vertices.push_back(vertices[v]);
    Hyperplane hp(hp_vertices);

    std::vector<int> rem =
      utils::remove_sorted(utils::range(0, n).tb(), vertex_sel);
    // OPA_DISP0(rem, vertex_sel, hp_vertices);

    std::vector<QModule_t> rem_vertices;
    std::vector<int> final_vertex_sel = vertex_sel;
    for (auto &id : rem) {
      if (hp.is_on(vertices[id]))
        final_vertex_sel.push_back(id);
      else
        rem_vertices.push_back(vertices[id]);
    }
    if (!hp.reorient_to_contain(rem_vertices)) return;
    maybe_add_face(final_vertex_sel, hp);
  }

  void maybe_add_face(const std::vector<int> &vertex_sel,
                      const Hyperplane &hp) {
    std::vector<int> lst = vertex_sel;
    std::sort(ALL(lst));
    if (seen.insert(lst).second) {
      faces_hp.push_back(hp);
      faces.push_back(lst);
    }
  }

  bool contains(const QModule_t &pt) const {
    return std::all_of(ALL(faces_hp), std::bind(&Hyperplane::contains,
                                                std::placeholders::_1, pt));
  }

  std::vector<QModule_t> vertices;
  std::vector<std::vector<int> > faces;
  std::vector<Hyperplane> faces_hp;
  std::set<std::vector<int> > seen;
  int n, d;
  Polytope_ConeRepr cone_repr;
};

template <class U> struct LexCompare {
  bool operator()(const U &a, const U &b) const {
    return std::lexicographical_compare(ALL(a), ALL(b));
  }
};

class Polyhedra_RelRepr {
public:
  void init() {
    if (is_init) return;
    is_init = true;
    OPA_CHECK0(rel_ineq.size() > 0);
    n = rel_ineq.size();
    d = rel_ineq[0].dim();

    OPA_TRACE("Finding vertex of polyhedron. ", n, d);
    // TODO: compute faces directly here, instead of redoing work in
    // VertexRepr
    OR::SubsetSelector<int>(utils::range(0, n).tb(),
                            [&](const std::vector<int> &sel) {
                              this->maybe_add_vertex(sel);
                              return OR::SubsetSelectorRes{ false, true };
                            },
                            d);
    OPA_TRACE("Found vertex repr > >", rel_ineq, vertices);
    std::sort(ALL(vertices), LexCompare<QModule_t>());
    utils::make_unique(vertices);
    vertex_repr.init(vertices);
  }

  void maybe_add_vertex(const std::vector<int> &hps_sel) {
    Matrix<Q> mat(&QF, d, d);
    std::vector<Q> pt;
    REP (i, hps_sel.size()) {
      mat.set_row(i, rel_ineq[hps_sel[i]].w);
      pt.emplace_back(rel_ineq[hps_sel[i]].level);
    }
    Matrix<Q> imat;
    if (!mat.invert(&imat)) return;
    QModule_t ipt(imat.eval(pt));

    for (auto &hp : rel_ineq)
      if (!hp.contains(ipt)) return;
    vertices.push_back(ipt);
    return;
  }

  const Polyhedra_VertexRepr &get_vertex_repr() {
    init();
    return vertex_repr;
  }

  std::vector<QModule_t> vertices;
  std::vector<Hyperplane> rel_ineq;
  int n, d;
  Polyhedra_VertexRepr vertex_repr;

  bool is_init = false;
};

Polyhedra_VertexRepr
compute_poly_proj(const Polyhedra_VertexRepr &polyhedra,
                  const std::vector<Hyperplane> &proj_planes) {

  std::vector<Hyperplane> hps;
  Polyhedra_VertexRepr cur = polyhedra;
  const int dim = polyhedra.d;
  REP (i, proj_planes.size()) {
    const int cur_dim = dim - i;
    auto proj_mat = expand_unimodular(&QF, proj_planes[i].w.elems, cur_dim - 1);
    Polyhedra_RelRepr next_repr;
    OPA_DISP0(proj_mat, proj_planes[i]);
    for (auto &hp : polyhedra.faces_hp) {
      auto tmp = proj_mat.evalT(hp.w);
      Q level_shift = tmp.back() * proj_planes[i].level;
      tmp.pop_back();
      Hyperplane new_hp(tmp, hp.level - level_shift);
      next_repr.rel_ineq.push_back(new_hp);
    }
    cur = next_repr.get_vertex_repr();
  }
  return cur;
}

class Simplex {
public:
  Simplex(const std::vector<QModule_t> &pts) {
    n = pts[0].elems.size();
    this->pts = pts;
    OPA_CHECK_EQ0(n + 1, pts.size());

    REP (i, n + 1) {
      std::vector<QModule_t> lst;
      REP (j, n + 1)
        if (i != j) lst.push_back(pts[j]);
      hps.emplace_back(lst);
      hps.back().reorient_to_contain(pts[i]);
    }
  }

  bool contains(const QModule_t &x) const {
    return std::all_of(
      ALL(hps), std::bind(&Hyperplane::contains, std::placeholders::_1, x));
  }

  int n;
  std::vector<QModule_t> pts;
  std::vector<Hyperplane> hps;
};

class Cone {
public:
  Cone() {}
  void init(const std::vector<QModule_t> &basis) { this->basis = basis; }
  std::vector<QModule_t> basis;
};

class SimplicialCone {
public:
  SimplicialCone() {}
  SimplicialCone &operator=(const SimplicialCone &peer) {
    planes = peer.planes;
    n = peer.n;
    m = peer.m;
    rays = peer.rays;
    iq = peer.iq.clone();
    diag = peer.diag;
    mat = peer.mat.clone();
    coord_map = peer.coord_map.clone();
    parallelepiped_count = peer.parallelepiped_count;
    return *this;
  }

  void init(const std::vector<ZModule_t> &rays) {
    Matrix<Z> mat(&Ring_Z, rays[0].elems.size(), rays.size());
    REP (i, rays.size())
      mat.setCol(i, rays[i].elems);
    this->init(mat);
  }

  void init(const std::vector<QModule_t> &rays) {
    std::vector<ZModule_t> zrays;
    for (auto &ray : rays) {
      zrays.push_back(force_zmodule(ray));
    }
    this->init(zrays);
  }

  void init(const Matrix<Z> &mat) {
    this->mat = mat.clone();
    Matrix<Q> qmat = mat.lift(QF);
    OPA_CHECK(qmat.left_inverse(&coord_map), mat);
    n = mat.n;
    m = mat.m;
    REP (i, m)
      rays.emplace_back(mat.get_col(i).tovec());

    OPA_CHECK(m <= n, m, n, mat);

    Matrix<bignum> p, d, q;
    // OPA_DISP0("COMPUTE SNF >> ", mat);
    snf(mat, &d, &p, &q);
    // puts("DONE");
    OPA_CHECK(q.invert(&iq), q);

    parallelepiped_count = 1;
    REP (i, m) {
      diag.push_back(d.get(i, i).abs());
      parallelepiped_count *= diag.back();
    }
  }

  void compute_shift_vertex(const QModule_t &vertex, ZModule_t *z_vertex,
                            QModule_t *q_shift) const {
    *z_vertex = ZModule_t(vertex.size());
    *q_shift = QModule_t(vertex.size());
    REP (i, vertex.size()) {
      (*z_vertex)[i] = QF.integer_part(vertex[i]);
      (*q_shift)[i] = QF.rational_part(vertex[i]);
    }
  }

  QModule_t get_coords(const ZModule_t &tb) const {
    return get_coords(lift_zmodule(tb));
  }
  QModule_t get_coords(const QModule_t &tb) const { return coord_map.eval(tb); }

  bool in_scone(const ZModule_t &tb) const {
    auto x = get_coords(tb);
    return std::all_of(ALL(x.elems), std::bind(&Q::operator>=, QF.getZ(),
                                               std::placeholders::_1));
  }

  bool in_cone(const ZModule_t &tb) const {
    auto x = get_coords(tb);
    return std::all_of(ALL(x.elems), std::bind(&Q::operator<=, QF.getZ(),
                                               std::placeholders::_1));
  }

  ZModule_t remap_parallelepiped(const ZModule_t &tb,
                                 const std::vector<bool> &exclude_face,
                                 const QModule_t &qs_coord) const {
    QModule_t coords = get_coords(tb) - qs_coord;

    std::vector<Z> rem;
    REP (i, coords.size()) {
      rem.push_back(QF.integer_part(coords[i]));
      if (exclude_face[i] && QF.is_integer(coords[i])) rem[i] -= 1;
    }

    ZModule_t remv(mat.eval(rem));
    ZModule_t rmp = tb - remv;

    // OPA_DISP("remapping >> ", coords, get_coords(tb), rmp, rem);
    {
      QModule_t final_coords = get_coords(rmp);
      REP (i, final_coords.size()) {
        if (exclude_face[i]) {
          OPA_CHECK(final_coords[i] > qs_coord[i] &&
                      final_coords[i] <= qs_coord[i] + QF.getE(),
                    final_coords, rmp, coords, rem, remv);
        } else {
          OPA_CHECK(final_coords[i] >= qs_coord[i] &&
                      final_coords[i] < qs_coord[i] + QF.getE(),
                    final_coords, rmp, coords, rem, remv);
        }
      }
    }
    return rmp;
  }

  std::vector<ZModule_t> list_parallelepiped() const {
    std::vector<bool> exclude_face(m, false);
    return list_parallelepiped(QModule_t::zero(n), exclude_face);
  }

  std::vector<ZModule_t>
  list_parallelepiped(const QModule_t &qshift,
                      const std::vector<bool> &exclude_face) const {

    std::vector<std::vector<bignum> > ranges;
    REP (i, n) {
      ranges.push_back(
        utils::Range<bignum>::StepRange(0, utils::get_or(diag, i, bignum(1)), 1)
          .tb());
    }

    std::vector<bool> mod_exclude_face;
    QModule_t qs_coord = get_coords(qshift);

    std::vector<ZModule_t> points;
    OR::CrossProdGen<bignum>(ranges, [&](const std::vector<bignum> &v) {
      auto lat_point = iq.eval(v);
      auto lat_mod = remap_parallelepiped(lat_point, exclude_face, qs_coord);
      points.push_back(lat_mod);
    });
    return points;
  }

  Matrix<Z> find_nw_matrix() const {
    Matrix<Q> mq = mat.lift(QF).inverse().transpose();
    Matrix<Z> W = mq.cmul(QF.import(parallelepiped_count)).project(QF);
    OPA_DISP0(W);

    // fplll considers row lattices
    fplll::ZZ_mat<mpz_t> fmat;
    fmat.resize(m, n);
    REP (i, n)
      REP (j, m)
        fmat(i, j) = *(mpz_t *)W(i, j).get_internal();

    int prec = 0;
    int status = fplll::bkz_reduction(
      fmat, fplll::BKZ_SLD_RED, fplll::BKZ_DEFAULT, fplll::FT_DEFAULT, prec);
    OPA_CHECK0(status == fplll::RED_SUCCESS);
    OPA_DISP0(W, mat);

    Matrix<Z> res(&Ring_Z, m, n);
    REP (i, m) {
      REP (j, n)
        fmat(i, j).get_mpz(*(mpz_t *)res(j, i).get_internal());
    }
    return res;
  }

  bool is_on_cone_boundary(const ZModule_t &pt) const {
    return std::any_of(
      ALL(planes),
      std::bind(&Hyperplane::is_on, std::placeholders::_1, lift_zmodule(pt)));
  }

  std::pair<bool, ZModule_t> find_w() const {

    Matrix<Z> wmat = find_nw_matrix();
    REP (i, m) {
      ZModule_t l(wmat.get_col(i));
      ZModule_t w(mat.eval(l));
      gcd_simplify(Ring_Z, w.elems);

      if (this->in_scone(w.elems)) {
        w = -w;
      }
      OPA_DISP("Trying ", w);
      REP (j, n) {
        Matrix<Z> nmat = mat.clone();
        nmat.set_col(j, w);
        OPA_DISP0(j, hnf_det(nmat));
        // if (hnf_det(nmat) == 0) goto bad;
      }
      // OPA_CHECK(!is_on_cone_boundary(w), mat);
      // if (is_on_cone_boundary(w)) continue;

      return { true, w };
    }
    OPA_DISP0(wmat);
    return { false, ZModule_t() };
  }

  std::vector<bool> compute_exclude_face(const ZModule_t &w) const {
    compute_planes();
    std::vector<bool> res;
    QModule_t qw = lift_zmodule(w);
    for (const auto &plane : planes) {
      res.push_back(plane.contains(qw));
    }
    return res;
  }

  void compute_planes() const {
    if (!planes.empty()) return;
    OPA_CHECK_EQ0(m, n);
    REP (i, m) {
      std::vector<QModule_t> lst;
      REP (j, m)
        if (j != i) lst.emplace_back(lift_zmodule(mat.get_col(j).tovec()));
      Hyperplane hp(lst, Hyperplane::kLinear);
      hp.reorient_to_contain(lift_zmodule(this->rays[i]));
      planes.push_back(hp);
    }
  }

  OPA_DECL_COUT_OPERATOR2(SimplicialCone, a.n, a.m, a.rays,
                          a.parallelepiped_count, a.diag, a.mat, a.planes);

  mutable std::vector<Hyperplane> planes;
  int n, m;
  std::vector<ZModule_t> rays;
  Matrix<Z> iq;
  std::vector<Z> diag;
  Matrix<Z> mat;
  Matrix<Q> coord_map;
  Z parallelepiped_count;
};

/*
class Lattice {
public:
  Lattice() {}
  static Lattice FromElems(const std::vector<ZModule_t> &elems) {
    Matrix<Z> mat;
    int n = elems.size();
    int m = elems[0].elems.size();
    mat.initialize(&Ring_Z, m, n);
    REP (i, n) { mat.setCol(i, elems[i].elems); }

    Matrix<Z> hnf_mat;
    hnf(mat, &hnf_mat, nullptr, nullptr);

    int lattice_n = 0;
    REP (i, n) {
      bool is_z = true;
      REP (j, m)
        if (hnf_mat(j, i) != 0) {
          is_z = false;
          break;
        }
      if (is_z) break;
      ++lattice_n;
    }

    Lattice res;
    res.full_dim = lattice_n == m;
    res.lat = hnf_mat.get_submatrix(0, 0, m, lattice_n).clone();
    if (res.full_dim) {
      res.det = 1;
      REP (i, m)
        res.det *= res.lat(i, i);
    }
    return res;
  }

  std::pair<bool, std::vector<Z> > expl(const ZModule_t &e) const {
    std::vector<Z> res;
    std::vector<Z> cur = e.elems;
    bool ok = reduce_hnf_vec(lat, cur, &res);
    return { ok, res };
  }

  bool full_dim;
  Matrix<Z> lat;
  Z det;
};
*/

class TaylorBase {
public:
  TaylorBase(Q startv) {
    this->cur = startv;
    termlist.push_back(get());
  }
  virtual ~TaylorBase() {}

  virtual Q _next() = 0;
  void next() {
    cur = _next();
    termlist.push_back(get());
    ++curpw;
  }

  Q get() const { return cur; }
  Q get_next() {
    next();
    return get();
  }

  P_Q get_as_poly(int maxterm) {
    while (termlist.size() <= maxterm) next();
    return Q_x.import_vec(termlist);
  }
  std::vector<Q> termlist;

  Q cur;
  int curpw = 1;
};

class TaylorExp : public TaylorBase {
public:
  // taylor expansion of e^(v*x), about x = 0
  TaylorExp(Z v) : TaylorBase(QF.getE()) { this->v = v; }

  virtual Q _next() override {
    Q res = cur;
    res.p *= v;
    res.q *= bignum(curpw);
    return QF.reduce(res);
  }

  Z v;
};

class Taylor_x_over_1_eminusx : public TaylorBase {
public:
  Taylor_x_over_1_eminusx(int maxn) : TaylorBase(QF.getE()), maxn(maxn) {
    ch.init(maxn + 3);
    clist.push_back(1);
  }

  virtual Q _next() override {
    int sgn = 1;

    Z nextc = 0;
    FOR (j, 1, curpw + 1) {
      Z term = ch.cnk(curpw + 1, j + 1) *
               (ch.fact[curpw] / ch.fact[curpw + 1 - j]) * clist[curpw - j];
      nextc += term * sgn;
      sgn *= -1;
    }
    clist.push_back(nextc);
    Q res = QF.import(nextc, ch.fact[curpw] * ch.fact[curpw + 1]);
    return res;
  }
  std::vector<Z> clist;
  CombHelper ch;
  int maxn;
};
enum BarvinokErrType {
  NOT_INTEGER = 1,
  BAD_VECTOR = 2,
  TOO_MANY_RETRIES = 3,
};

ZModule_t generate_lambda_rand(int dim) {
  bignum gen_bound = bignum(2).lshift(30);
  ZModule_t cnd;
  REP (i, dim)
    cnd.elems.push_back(gen_bound.rand_signed());
  return cnd;
}

class BarvinokResidueHelper {
public:
  struct ConeCtx {
    P_Q dirs_poly;
    Z vertex_dot;
    P_Q points_poly;
    Q resv;
  };

  struct Result {
    BarvinokErrType err;
    bool ok;
    bool fatal;
    Q qv;
    Z v;
    ZModule_t bad_dir;
    OPA_DECL_COUT_OPERATOR2(Result, a.ok, a.err, a.v, a.qv, a.bad_dir);
  };

  void reset(int dim) {
    this->dim = dim;
    taylor_dirs_base = Taylor_x_over_1_eminusx(dim).get_as_poly(dim);
    res.ok = true;
    res.fatal = false;
    cone_ctx.resv = QF.getZ();
    data = ResidueData();
  }

  struct ResidueCone {
    std::vector<ZModule_t> points;
    ZModule_t vertex;
    std::vector<ZModule_t> dirs;
    int sgn;
    OPA_DECL_COUT_OPERATOR2(ResidueCone, a.vertex, a.dirs, a.points);
  };

  struct ResidueData {
    std::vector<ResidueCone> cones;
    OPA_DECL_COUT_OPERATOR2(ResidueData, a.cones);
  };

  ResidueCone &active() { return data.cones.back(); }

  void start_cone(const SimplicialCone &cone, int sgn,
                  const ZModule_t &vertex) {
    if (online_mode) {
      OPA_CHECK0(lambda_set);
      std::vector<ZModule_t> dirs;
      REP (i, cone.m)
        dirs.emplace_back(cone.mat.get_col(i));
      cone_ctx.dirs_poly = compute_dirs_poly(dirs, sgn);
      cone_ctx.points_poly = Q_x.getZ();
      cone_ctx.vertex_dot = vertex.dot(selected_lambda);

    } else {
      data.cones.emplace_back();
      active().vertex = vertex;
      active().sgn = sgn;
      REP (i, cone.m)
        active().dirs.emplace_back(cone.mat.get_col(i));
    }
  }

  void finish_cone() {
    if (!res.ok) return; // Poor man's exceptions. Maybe reconsider.
    if (online_mode) {
      Q qs = compute_final_res(cone_ctx.points_poly, cone_ctx.dirs_poly);
      OPA_DISP0("BEFORE >> ", cone_ctx.resv);
      cone_ctx.resv = cone_ctx.resv + qs;
      OPA_DISP0("FINISHING CONE >> ", cone_ctx.resv, qs);
    }
  }

  void add_point(const ZModule_t &pt) {
    if (online_mode) {
      cone_ctx.points_poly = cone_ctx.points_poly + get_point_poly(pt);
    } else {
      active().points.push_back(pt);
    }
  }

  bool is_good_lambda(const ZModule_t &cnd) const {
    for (auto &cone : data.cones) {
      for (auto &dir : cone.dirs) {
        if (dir.dot(cnd) == 0) return false;
      }
    }
    return true;
  }

  ZModule_t find_lambda_rand() const {
    while (true) {
      ZModule_t cnd = generate_lambda_rand(dim);
      if (is_good_lambda(cnd)) return cnd;
    }
  }

  P_Q get_point_poly(const ZModule_t &pt) const {
    Z val = cone_ctx.vertex_dot + pt.dot(selected_lambda);
    // OPA_DISP0(pt, val);
    return TaylorExp(val).get_as_poly(dim);
  }

  P_Q compute_dirs_poly(const std::vector<ZModule_t> &dirs,
                        int cone_sgn) const {
    Z denom = cone_sgn;
    P_Q final_p = Q_x.getE();
    // why is denom out?
    for (const auto &dir : dirs) {
      Z val = dir.dot(selected_lambda);
      if (val == 0) {
        res.ok = false;
        res.bad_dir = dir;
        res.err = BarvinokErrType::BAD_VECTOR;
        return Q_x.getZ();
      }
      denom = denom * -val;

      P_Q tmp = taylor_dirs_base;
      tmp.mulx(QF.import(-val));
      final_p = Q_x.mulmoddeg(final_p, tmp, dim + 1);
    }
    return Q_x.mulc(final_p, QF.importq(denom));
  }

  Q compute_final_res(const P_Q &points_poly, const P_Q &dirs_poly) const {
    return Q_x.mulmoddeg(points_poly, dirs_poly, dim + 1).get_safe(dim);
    // not the fastest way to get [n], does not matter
  }

  Q compute_cone(const ResidueCone &cone) const {
    OPA_CHECK0(lambda_set);
    OPA_CHECK_EQ0(cone.dirs.size(), dim);

    P_Q points_poly = Q_x.getZ();
    for (const auto &pt : cone.points) {
      points_poly = points_poly + get_point_poly(pt);
    }

    P_Q dirs_poly = compute_dirs_poly(cone.dirs, cone.sgn);
    OPA_DISP0(points_poly, dirs_poly, selected_lambda, cone_ctx.vertex_dot);
    return compute_final_res(points_poly, dirs_poly);
  }

  void set_lambda(const ZModule_t &lambda) {
    selected_lambda = lambda;
    lambda_set = true;
  }

  bool push_cone(const SimplicialCone &cone, int sgn,
                 const std::vector<ZModule_t> &cone_points,
                 const ZModule_t &vertex) {
    start_cone(cone, sgn, vertex);

    for (const auto &e : cone_points) {
      if (!res.ok) break;
      add_point(e);
    }
    finish_cone();
    return res.ok;
  }

  const Result &compute() {
    if (!online_mode) {
      if (!lambda_set) {
        set_lambda(find_lambda_rand());
      }
      Q val = QF.getZ();
      for (auto &cone : data.cones) {
        if (!res.ok) break;
        cone_ctx.vertex_dot = cone.vertex.dot(selected_lambda);
        val = val + compute_cone(cone);
      }
      cone_ctx.resv = val;
    }

    res.qv = cone_ctx.resv;

    if (!res.ok)
      ;
    else if (!QF.is_integer(res.qv)) {
      res.ok = false;
      res.fatal = true;
      res.err = BarvinokErrType::NOT_INTEGER;
    } else {
      res.v = QF.to_base_or_fail(res.qv);
    }
    return res;
  }

  mutable Result res;

  ConeCtx cone_ctx;

  P_Q taylor_dirs_base;
  bool online_mode = true;
  int dim;
  bool lambda_set = false;
  ZModule_t selected_lambda;
  ResidueData data;
};

class Barvinok : opa::utils::Initable {
public:
  typedef std::vector<bool> ExcludeFace_t;
  struct RecData {
    QModule_t vertex;
    SimplicialCone cone;
    ExcludeFace_t exclude_face;
    int sgn;
    int depth;
    OPA_DECL_COUT_OPERATOR2(RecData, a.vertex, a.cone, a.exclude_face, a.sgn,
                            a.depth);
  };

  struct Result {
    bool ok;
    BarvinokErrType err;
    BarvinokResidueHelper::Result residue_res;
    Z val;
    OPA_DECL_COUT_OPERATOR2(Result, a.ok, a.err, a.val, a.residue_res);
  };

  struct Configuration {
    bignum enumerate_ub;
    int max_depth = INT_MAX;
    int max_try = INT_MAX;
    bool online_mode = true;
    bool fixed_lambda = false;
    bool debug = false;
    ZModule_t lambda;
  };

  struct DescentData {
    std::vector<std::pair<bignum, bignum> > transitions;
    bool ok = true;
    BarvinokResidueHelper::Result residue_res;
    std::vector<DescentData> children;
    RecData rd;
    std::vector<int> hp_status;
    ZModule_t w;

    void reset() {
      transitions.clear();
      ok = true;
    }
  };

  bool should_enumerate_cone(const RecData &data) const {
    if (data.cone.parallelepiped_count <= conf.enumerate_ub) return true;
    if (data.depth >= conf.max_depth) return true;
    return false;
  }

  bool rec(const RecData &data, DescentData &descent_data) {
    if (conf.debug) {
      descent_data.rd = data;
    }

    if (should_enumerate_cone(data)) {
      QModule_t qs_shift;
      ZModule_t zvertex;
      data.cone.compute_shift_vertex(data.vertex, &zvertex, &qs_shift);

      auto lst = data.cone.list_parallelepiped(qs_shift, data.exclude_face);
      OPA_DISP("Leaf enumeration >> ", data.sgn, data.exclude_face, data.cone,
               zvertex, qs_shift);
      if (!residue_helper.push_cone(data.cone, data.sgn, lst, zvertex)) {
        descent_data.residue_res = residue_helper.res;
        return false;
      }

    } else {
      bool found;
      ZModule_t w;
      std::tie(found, w) = data.cone.find_w();
      OPA_CHECK(found, data.cone);
      int m = data.cone.m;

      std::vector<int> hp_status;
      data.cone.compute_planes();

      REP (i, m) {
        hp_status.push_back(
          data.cone.planes[i].contains_status(lift_zmodule(w)));
        OPA_DISP0(
          data.cone.planes[i].contains_status(lift_zmodule(data.cone.rays[i])));
      }
      OPA_DISP("Splitting ", data.cone, w, hp_status, data.exclude_face,
               data.cone.planes);
      // a,b,false: exclude face when drop ray a, then ray b
      // a,b,true: exclude face face intersection missing (ray a and w)
      algo::Sat2<std::tuple<int, int, bool> > sat2;
      std::vector<int> next_sgns(m);
      std::vector<bool> exclude_w(m);

      REP (ignorev, m) {
        if (hp_status[ignorev] == 0) continue;
        int nsgn = hp_status[ignorev] * data.sgn;
        next_sgns[ignorev] = nsgn;
        exclude_w[ignorev] =
          (data.exclude_face[ignorev] ^ (hp_status[ignorev] == -1));

        REP (j, m) {
          if (j == ignorev || hp_status[j] == 0) continue;

          sat2 += sat2.term(ignorev, j, true) ==
                  (sat2.term(ignorev, j, false) || exclude_w[ignorev]);
          if (j < ignorev) {
            // face without ray ignorev,j should be:
            // ++ / -- -> only 1
            // +- -> both set or both unset
            sat2 +=
              (sat2.term(ignorev, j, false) ^ sat2.term(j, ignorev, false)) ==
              (hp_status[ignorev] == hp_status[j]);

            // Number of time the face-face intersection is counted should not
            // change.
            bool target = !(data.exclude_face[j] || data.exclude_face[ignorev]);
            sat2 += (!sat2.term(j, ignorev, true) * (s8)hp_status[ignorev] +
                     !sat2.term(ignorev, j, true)) *
                      (s8)hp_status[j] ==
                    target;
          }
        }
      }
      OPA_CHECK0(sat2.compute());

      if (conf.debug) {
        descent_data.w = w;
        descent_data.hp_status = hp_status;
      }

      REP (ignorev, m) {
        if (hp_status[ignorev] == 0) {
          OPA_DISP("Skipping ", ignorev, data.cone.rays[ignorev], w);
          continue;
        }
        RecData ndata;
        std::vector<ZModule_t> nrays;
        REP (i, m)
          if (i != ignorev) nrays.push_back(data.cone.rays[i]);
        nrays.push_back(w);

        ndata.cone.init(nrays);
        OPA_DISP("NCONE >>> ", ndata.cone.parallelepiped_count);
        ndata.depth = data.depth + 1;
        ndata.vertex = data.vertex;

        ndata.sgn = next_sgns[ignorev];
        ndata.exclude_face.resize(m, false); // include by default
        ndata.exclude_face[m - 1] = exclude_w[ignorev];

        REP (i, m) {
          if (i == ignorev) continue;
          int id = i - (i >= ignorev);
          if (hp_status[i] == 0) {
            // w is on face i, so its new face without i should respect
            // data.exclude_face.
            ndata.exclude_face[id] = data.exclude_face[i];
          } else {
            ndata.exclude_face[id] = sat2.get2(ignorev, i, false);
          }
        }

        descent_data.transitions.emplace_back(data.cone.parallelepiped_count,
                                              ndata.cone.parallelepiped_count);
        descent_data.children.emplace_back();
        DescentData &nd = descent_data.children.back();
        if (!this->rec(ndata, nd)) return false;
      }
    }
    return true;
  }

  std::vector<RecData> prepare_cone(const PolyhedraCone &cone) {

    int dim = cone.rays[0].elems.size();

    std::vector<SimplicialCone> scones;

    ZModule_t Y;
    std::vector<ZModule_t> normed_rays = cone.rays;
    for (auto &ray : normed_rays) {
      Z gcd = gcd_simplify(Ring_Z, ray.elems);
      if (gcd < 0)
        for (auto &e : ray.elems) e *= -1;
    }
    OPA_DISP("Normed rays >> ", normed_rays, cone.rays);

    std::vector<PolyhedraSimplicialDecomposition> decompositions =
      cone.simplicial_decomposition;
    if (decompositions.empty()) {
      decompositions.push_back(PolyhedraSimplicialDecomposition{
        utils::range(0, normed_rays.size()).tb() });
    }

    bool first = true;
    std::vector<bool> toggle_exclude_face;
    int pos = 0;

    std::vector<RecData> data_to_call;
    for (auto &decomp : decompositions) {
      OPA_CHECK0(decomp.ray_ids.size() == dim);
      Matrix<Z> mat(&Ring_Z, dim, dim);
      REP (i, decomp.ray_ids.size()) {
        mat.setCol(i, normed_rays[decomp.ray_ids[i]].elems);
      }
      RecData data;
      data.cone.init(mat);

      if (first) {
        first = false;
        Y = ZModule_t(mat.eval(mat.ones(dim, 1)));
        toggle_exclude_face = data.cone.compute_exclude_face(Y);
      }

      data.exclude_face = data.cone.compute_exclude_face(Y);
      REP (i, toggle_exclude_face.size())
        data.exclude_face[i] = data.exclude_face[i] ^ toggle_exclude_face[i];
      data.vertex = cone.vertex;
      data.sgn = 1;
      data.depth = 0;
      data_to_call.push_back(data);
    }
    return data_to_call;
  }

  bool add_cone(const PolyhedraCone &cone) {
    std::vector<RecData> data_to_call = prepare_cone(cone);
    descent_data.reset();
    for (auto &data : data_to_call) {
      descent_data.children.emplace_back();
      DescentData &nd = descent_data.children.back();
      // if (pos++ != 1) continue;
      if (!this->rec(data, nd)) return false;
    }
    return true;
  }

  void init(const Configuration &conf) {
    this->conf = conf;
    utils::Initable::init();
  }

  void reset_residue_helper(int dim) {
    residue_helper.reset(dim);
    residue_helper.online_mode = conf.online_mode;
    if (conf.fixed_lambda) {
      residue_helper.set_lambda(conf.lambda);
    } else {
      residue_helper.set_lambda(generate_lambda_rand(dim));
    }
  }

  Result try_solve(const Polytope_ConeRepr &polytope) {
    utils::Initable::check_init();
    int max_try = !conf.online_mode || conf.fixed_lambda ? 1 : conf.max_try;

    Result fres;
    fres.ok = true;
    int dim = polytope.cones[0].vertex.size();

    REP (i, max_try) {
      reset_residue_helper(dim);
      for (const auto &cone : polytope.cones) {
        if (!this->add_cone(cone)) break;
      }
      auto res = residue_helper.compute();
      if (!res.ok) {
        OPA_DISP0(res);
        if (res.fatal) {
          fres.err = res.err;
          fres.residue_res = res;
          fres.ok = false;
          return fres;
        }
        continue;
      }
      fres.val = res.v;
      fres.residue_res = res;
      return fres;
    }

    fres.ok = false;
    fres.err = BarvinokErrType::TOO_MANY_RETRIES;
    return fres;
  }

  Z solve(const Polytope_ConeRepr &polytope) {
    Result res = this->try_solve(polytope);
    OPA_CHECK(res.ok, res);
    return res.val;
  }

  Z solve_simplex(const Simplex &simplex) {
    Polytope_ConeRepr polytope;
    int m = simplex.pts.size();
    REP (i, m) {
      PolyhedraCone cone;

      cone.vertex = simplex.pts[i];
      REP (j, m)
        if (i != j)
          cone.rays.push_back(force_zmodule(simplex.pts[j] - simplex.pts[i]));
      polytope.cones.push_back(cone);
    }
    return this->solve(polytope);
  }

  Configuration conf;
  BarvinokResidueHelper residue_helper;
  DescentData descent_data;
};

OPA_MATH_CO_END
