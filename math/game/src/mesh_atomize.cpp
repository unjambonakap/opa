#include <opa/math/game/mesh_atomize.h>
#include <opa/math/game/point_cloud.h>
#include <opa/math/game/quat.h>
#include <opa/or/adaptative_search.h>

#include <lemon/concepts/graph.h>
#include <lemon/list_graph.h>
#include <lemon/matching.h>
#include <lemon/smart_graph.h>
#include <opa/algo/graph.h>
#include <opa/math/common/float.h>
#include <random>
#include <valarray>

std::mt19937 rng(0);

OPA_NAMESPACE_DECL3(opa, math, game)
namespace atom {
const double eps = 1e-6;

enum ScorePriority {
  Normal = 0,
  ExtendDegenerateBox = 1,
  FixDegenerateBox = 2,
  Free = 3,
};

template <typename Field> struct FieldCompare {
  Field field;
  bool ord = false;
  template <class T> bool operator()(const T &a, const T &b) const {
    OPA_DISP0(a.dest_cluster, b.dest_cluster);
    return (a.*field < b.*field) ^ ord;
  }
};
template <typename Field>
FieldCompare<Field> field_compare(Field field, bool ord = false) {
  return FieldCompare<Field>{ field, ord };
}

template <typename Method> struct MemberCompare {
  Method method;
  bool ord = false;

  template <class T> bool operator()(const T &a, const T &b) const {
    return (a.*method() < b.*method()) ^ ord;
  }
};
template <typename Member>
MemberCompare<Member> member_compare(Member member, bool ord = false) {
  return MemberCompare<Member>{ member, ord };
}

class VariableMarkovChain {
public:
  VariableMarkovChain(int n, ProbHelper &helper) : n(n), helper(helper) {
    mat.resize(n);
    REP (i, n)
      mat[i] = std::vector<double>(n, 0.);
  }

  int process(double self_state_prob) {
    mat[state][state] = 0;
    if (helper.select_one(self_state_prob)) {
      this->change_state = false;
    } else {
      state = helper.select_random_id_double(mat[state]);
      this->change_state = true;
    }
    return state;
  }

  std::vector<std::vector<double> > mat;
  int n;
  ProbHelper &helper;
  bool change_state = true;
  int state = 0;
};

KAtomizer::KAtomizer(const Mesh &mesh, int K) : K(K) {
  this->mesh = mesh;
  this->mesh.whiten();
  OPA_CHECK0(mesh.is_connected());
  FaceCollection fc(true);
  mesh.to_faces(&fc);
  for (auto &f : fc.faces()) {
    triangles.push_back(vec_to_tr(f));
  }
  mesh_area = mesh.area();
  face_graph = mesh.compute_face_graph2();
}

BoxSpec KAtomizer::compute_box(const TrSet &tr_ids) const {
  PointVec tb = get_points(tr_ids);
  // OPA_DISP0(tb, tr_ids);
  BoxSpec res = compute_best_box(tb).expand(1.01);
  for (auto &pt : tb) {
    OPA_CHECK(res.dist2(pt) < 1e-8, res.str(), pt);
  }
  return res;
}

PointVec &KAtomizer::add_points(PointVec &res, const TrSet &tr_ids) const {
  for (auto &id : tr_ids) {
    res.push_back(std::get<0>(triangles[id]));
    res.push_back(std::get<1>(triangles[id]));
    res.push_back(std::get<2>(triangles[id]));
  }
  return res;
}

PointVec KAtomizer::get_points(const TrSet &tr_ids) const {
  PointVec res;
  add_points(res, tr_ids);
  return res;
}

void KAtomizer::add_face_to_cluster(ClusterData *cluster, int face_pos) {
  cluster->faces.push_back(face_pos);
  cluster->bfs_ctx.sources = { face_pos };
  face_graph->graph->bfs_one(cluster->bfs_ctx);
  cluster->box = compute_box(cluster->faces);
}

RoundData round_data_from_cluster(const ClusterData &cluster) {
  RoundData res;
  res.rot = rot_to_box_space(cluster.box);
  res.box = BoxAASpec::FromBoxSpace(cluster.box);
  return res;
}

void KAtomizer::setup_candidate(const GlobalRoundData &grd,
                                const ClusterData &cl, CandidateStruct &cnd) {
  const RoundData &rd = grd.rd[cnd.dest_cluster];
  double MIN_VAL = 1e-9;

  cnd.priority = ScorePriority::Normal;
  PointVec tr = tr_to_vec(triangles[cnd.face]);
  {
    bool is_free = true;
    REP (j, 3)
      if (!cl.box.in(tr[j])) {
        is_free = false;
        break;
      }
    if (is_free) {
      cnd.priority = ScorePriority::Free;
      return;
    }
  }

  BoxAASpec nbox = rd.box;
  REP (j, 3)
    nbox.update(rd.rot * (tr[j] - cl.box.corner));
  double area_change_ratio = nbox.area() - rd.box.area();

  std::vector<double> dims(3);
  REP (j, 3)
    dims[j] = nbox.get_range(j).get_length();
  double minv = std::max<double>(base_eps, *std::min_element(ALL(dims)));
  double maxv = *std::max_element(ALL(dims));
  double box_like_score = std::min<double>(1000., maxv / minv);

  std::vector<double> box_areas;
  box_areas.push_back(nbox.area());
  REP (j, K) {
    if (cnd.dest_cluster == j) continue;
    box_areas.push_back(grd.rd[j].box.area());
  }
  // TODO: maybe score disparity
  std::valarray<double> box_array(box_areas.data(), box_areas.size());
  double area_sum = box_array.sum();

  double nmax = std::max<double>(nbox.area(), grd.max_box_area_round);
  nmax = std::max(MIN_VAL, nmax);
  double all_box_score = nmax / grd.max_box_area_round;
  if (cnd.dest_cluster == 0) {
    OPA_DISP0(nbox.str(), nmax, grd.max_box_area_round, all_box_score);
  }

  if (cl.box.almost_empty()) {
    cnd.priority = ScorePriority::ExtendDegenerateBox;
    PointVec tb = get_points({ cnd.face });
    REP (k, 3)
      tb.push_back(cl.box.get(k));
    auto non_planar_quad = find_non_planar_point(tb);
    if (!non_planar_quad.empty()) {
      Dir plane_dir =
        get_plane_normal(tb[0], tb[non_planar_quad[0]], tb[non_planar_quad[1]]);
      double mx = 0;
      FOR (k, tb.size() - 3, tb.size())
        mx = std::max<double>(mx, std::abs(glm::dot(tb[k] - tb[0], plane_dir)));
      if (mx > 1e-2) {
        auto new_best_box = compute_best_box(tb);
        OPA_CHECK(!new_best_box.almost_empty(), tb, cnd.face, cl.box,
                  new_best_box, mx);
        cnd.priority = ScorePriority::FixDegenerateBox;
        OPA_DISP0(cnd.dest_cluster, mx, tb);
      }
    }
  }
  cnd.box_like_score = std::max(MIN_VAL, box_like_score);
  cnd.all_box_score = std::max(MIN_VAL, all_box_score);
  cnd.box_like_score = 1.;
  cnd.penalty = linearize(nmax / std::max<double>(MIN_VAL, nbox.area()), 2, 7, 10., 1);
  // cnd.all_box_score = 1.;
  cnd.area_change_ratio = std::max(MIN_VAL, area_change_ratio / nmax);
}

AtomizeAssignment KAtomizer::initial_assignment() {
  AtomizeAssignment res;
  res.cluster_data.resize(K);
  std::set<int> rem;
  REP (i, triangles.size())
    rem.insert(i);

  algo::FastGraph::DijkstraCtx ctx;
  ctx.edge_costs = face_graph->edge_costs;
  face_graph->graph->dijkstra(ctx); // no source, used to init dists

  AssignmentDebug debug_data;

  REP (i, K) {
    std::vector<double> dist_tmp = ctx.dists;
    for (auto &x : dist_tmp) x *= x;
    int id = prob_helper.select_random_id_double_topr(dist_tmp, 0.7);
    OPA_DISP0(id, dist_tmp[id], dist_tmp);

    add_face_to_cluster(&res.cluster_data[i], id);
    rem.erase(id);

    ctx.sources = { id };
    // dists[id] will be 0, so that it does not get selected aftewards
    face_graph->graph->dijkstra(ctx);
  }

  if (debug) {
    InitRoundDebug round_debug;
    round_debug.cluster_data = res.cluster_data;
    debug_data.init_round_debug.push_back(round_debug);
  }

  while (rem.size() > 0) {
    GlobalRoundData grd;

    REP (i, K) {
      grd.rd.push_back(round_data_from_cluster(res.cluster_data[i]));
      OPA_DISP0(res.cluster_data[i].box.str(),
                res.cluster_data[i].box.corners());
    }

    grd.max_box_area_round =
      std::max_element(ALL(grd.rd), [](const RoundData &a, const RoundData &b) {
        return a.box.area() < b.box.area();
      })->box.area();
    grd.max_box_area_round = std::max(grd.max_box_area_round, base_eps);

    std::vector<CandidateStruct> cnds;
    REP (i, K) {
      const RoundData &rd = grd.rd[i];
      auto &cl = res.cluster_data[i];
      for (auto &x : cl.bfs_ctx.waiting) {
        if (!rem.count(x)) continue;

        CandidateStruct cnd;
        cnd.dest_cluster = i;
        cnd.face = x;
        setup_candidate(grd, cl, cnd);
        cnds.push_back(cnd);
      }
    }

    int max_priority =
      utils::MinFinder<int, std::greater<int> >()
        .update<CandidateStruct, &CandidateStruct::priority>(cnds)
        .get();

    cnds.erase(std::remove_if(ALL(cnds),
                              [&](const CandidateStruct &cnd) {
                                return cnd.priority != max_priority;
                              }),
               cnds.end());

    for (auto &cnd : cnds) {
      cnd.score = (1. / cnd.area_change_ratio) * (1. / cnd.box_like_score) *
                  (1. / cnd.all_box_score) * (1. / cnd.penalty);
      // OPA_DISP0(cnd);
    }
    // std::sort(ALL(cnds), field_compare(&CandidateStruct::score, true));
    std::sort(ALL(cnds));
    std::map<int, int> seen_disp;
    for (auto &cnd : cnds) {
      if (++seen_disp[cnd.dest_cluster] >= 8) continue;
      OPA_DISP0(cnd);
    }
    for (auto &cnd : cnds) {
      OPA_CHECK(cnd.dest_cluster >= 0 && cnd.dest_cluster < K, cnd);
    }

    OPA_CHECK0(cnds.size() > 0);
    CandidateStruct selected_cnd =
      prob_helper.select_random<CandidateStruct, &CandidateStruct::score>(cnds);
    OPA_DISP0(selected_cnd, get_points({ selected_cnd.face }));
    rem.erase(selected_cnd.face);
    add_face_to_cluster(&res.cluster_data[selected_cnd.dest_cluster],
                        selected_cnd.face);
    // OPA_CHECK0(selected_cnd.dest_cluster != 2);
    if (selected_cnd.priority == 2) {
      OPA_CHECK0(
        res.cluster_data[selected_cnd.dest_cluster].box.almost_empty() ==
        false);
    }

    if (debug) {
      InitRoundDebug round_debug;
      round_debug.cluster_data = res.cluster_data;
      round_debug.selected_candidate = selected_cnd;
      round_debug.candidates = cnds;
      debug_data.init_round_debug.push_back(round_debug);
    }
  }

  res.debug_data = debug_data;
  return res;
}

double KAtomizer::get_assignment_extra_area(
  const AtomizeAssignment &assignment) const {
  double extra_area = 0;
  for (auto &cluster : assignment.cluster_data) {
    extra_area += get_cluster_covered_area(cluster) - cluster.box.area();
  }
  return extra_area;
}

double KAtomizer::get_cluster_covered_area(const ClusterData &cluster) const {
  auto cut_meshes = mesh_cut_box(mesh, cluster.box);
  double area = 0;
  for (auto &mesh : cut_meshes) {
    area += mesh->area();
  }
  return area;
}

double
KAtomizer::get_assignment_cost(const AtomizeAssignment &assignment) const {
  double sum_area = 0;
  for (auto &x : assignment.cluster_data) sum_area += x.box.area();
  return sum_area;
}

bool KAtomizer::check_box_and_face(const BoxSpec &box,
                                   const TrSet &faces) const {
  auto pts = get_points(faces);
  for (auto &pt : pts) {
    OPA_CHECK(box.in(pt), pt);
    OPA_CHECK(box.dist2(pt) < 1e-6, pt);
  }

  return true;
}

bool KAtomizer::is_cluster_valid(const ClusterData &cluster, int mode) const {
  if (mode == 0) {
    auto cut_meshes = mesh_cut_box(mesh, cluster.box);
    if (cut_meshes.size() == 1) return true;
    std::vector<BoxSpec> cut_boxes;
    for (auto &cut_mesh : cut_meshes)
      cut_boxes.push_back(compute_best_box(cut_mesh->point_cloud()));

    utils::ConstantChecker<int> mesh_id;
    for (auto &face : cluster.faces) {
      bool fd = false;
      PointVec tr = tr_to_vec(triangles[face]);
      REP (i, cut_boxes.size()) {
        if (cut_boxes[i].in(tr[0]) && cut_boxes[i].in(tr[1]) &&
            cut_boxes[i].in(tr[2])) {
          fd = true;
          if (!mesh_id.add(i)) return false;
          break;
        }
      }
      OPA_CHECK0(fd);
    }
  } else if (mode == 1) {
    return check_box_and_face(cluster.box, cluster.faces);

  } else
    OPA_CHECK(false, mode);
  return true;
}

bool KAtomizer::is_assignment_valid(const AtomizeAssignment &assignment,
                                    int mode) const {
  REP (i, K)
    if (!is_cluster_valid(assignment.cluster_data[i], mode)) return false;
  return true;
}

std::vector<BoxSpec> atomize_mesh_K(const Mesh &mesh, int K) {
  utils::MinFinderPair<double, std::vector<BoxSpec> > finder;
  KAtomizer atomizer(mesh, K);

  int nstep = 40;
  REP (step, nstep) {
    auto ans = atomizer.initial_assignment();
    double max_area = 0.;
    double area_sum = 0;
    std::vector<BoxSpec> res;
    for (auto &x : ans.cluster_data) {
      max_area = std::max(max_area, x.box.area());
      area_sum += x.box.area();
      res.push_back(x.box);
    }
    finder.update(0 * max_area * max_area + area_sum, res);
  }
  return finder.get();
}

std::vector<BoxSpec> atomize_mesh(const Mesh &mesh) {
  std::vector<BoxSpec> res;
  return res;
}

void KAtomizer::improve_assignment(AtomizeAssignment *sol,
                                   const KAtomizer_SearchParams &params) {
  OptimizationContext ctx;
  ctx.stop_criteria.reset(params.stop_criteria->clone());
  ctx.stop_criteria->reset();
  ctx.temp.T = params.T0;

  ctx.score = 0;
  REP (i, K)
    ctx.score += sol->cluster_data[i].box.area();
  OPA_DISP("Initial score: ", ctx.score);

  std::array<double, 3> iir_coeffs;
  {
    double n = 5;
    double prob_at_n = 0.5;
    iir_coeffs[0] = std::pow(prob_at_n, 1. / n);
  }

  {
    double n = 30;
    double prob_at_n = 0.5;
    iir_coeffs[1] = std::pow(prob_at_n, 1. / n);
  }

  {
    double n = 30;
    double prob_at_n = 0.5;
    iir_coeffs[2] = std::pow(prob_at_n, 1. / n);
  }
  OPA_DISP0(iir_coeffs);

  int nstate = 3;
  VariableMarkovChain mc(nstate, this->prob_helper);
  REP (i, nstate)
    REP (j, nstate)
      if (i != j) mc.mat[i][j] = 1;

  // Starting state is A
  mc.state = params.state;
  const double iir_status_start = 1.;
  double iir_status = iir_status_start;

  while (!ctx.stop_criteria->should_stop()) {
    if (!params.force_state) mc.process(iir_status);
    if (mc.change_state) {
      iir_status = iir_status_start;
    }

    OPA_DISP("On ", ctx.stop_criteria->round(), mc.state, ctx.score,
             iir_status);
    ctx.history.push_back(0);
    if (mc.state == 0) {
      this->do_search_a(sol, &ctx);
    } else if (mc.state == 1) {
      this->do_search_b(sol, &ctx);
    } else {
      this->do_search_c(sol, &ctx);
    }
    OPA_CHECK0(is_assignment_valid(*sol, 1));

    iir_status = iir_status * iir_coeffs[mc.state] + ctx.history.back();

    ctx.temp.next();
    OR::SearchRoundDataForStopping stop_data;
    stop_data.new_best = ctx.score;
    stop_data.point_dist = -1;
    ctx.stop_criteria->new_round(stop_data);
  }
  OPA_DISP("Final score >> ", ctx.score);
}

void KAtomizer::do_search_a(AtomizeAssignment *assignment,
                            OptimizationContext *ctx) {
  std::vector<double> box_score;
  REP (i, K)
    box_score.push_back(assignment->cluster_data[i].box.area());
  int selected_box = prob_helper.select_random_id_double(box_score);

  auto &source_cl = assignment->cluster_data[selected_box];
  ProbHelperBuilder<std::pair<int, int> > tr_selector;
  REP (face_pos, source_cl.faces.size()) {
    int face = source_cl.faces[face_pos];
    PointVec tr_points = get_points({ face });
    Pos center = gravity_center(tr_points);
    double dist = source_cl.box.dist_to_center(center);
    REP (j, K) {
      if (j != selected_box) {
        // use volume extend metric + axis outlier (search_b)
        double score =
          dist / (1. + assignment->cluster_data[j].box.dist2(center));
        tr_selector.push(score, face_pos, j);
      }
    }
  }
  // select K candidates, to study transfer
  int selected_tr;
  int target_box;
  std::tie(selected_tr, target_box) = tr_selector.get(prob_helper);

  auto &target_cl = assignment->cluster_data[target_box];
  TrSet n_target_set = target_cl.faces;
  TrSet n_source_set = source_cl.faces;
  n_target_set.push_back(selected_tr);
  BoxSpec n_target_box = compute_box(n_target_set);
  OPA_CHECK(n_target_box.in(get_points({ selected_tr })),
            get_points({ selected_tr }), n_target_box.str());

  REP (i, n_source_set.size()) {
    int face = n_source_set[i];

    bool can_transfer = n_target_box.in(get_points({ face }));
    if (face == selected_tr) OPA_CHECK0(can_transfer);
    if (can_transfer) {
      n_target_set.push_back(face);
      std::swap(n_source_set.back(), n_source_set[i]);
      n_source_set.pop_back();
      --i;
    }
  }

  if (n_source_set.size() > 0) {
    BoxSpec n_source_box = compute_box(n_source_set);
    BoxSpec n_target_box = compute_box(n_target_set);
    double old_score = source_cl.box.area() + target_cl.box.area();
    double new_score = n_source_box.area() + n_target_box.area();
    double diff_score = new_score - old_score;
    double rel_score = diff_score / old_score;
    OPA_DISP("candidate ", n_source_set.size(), n_target_set.size(), old_score,
             new_score);

    bool accept = ctx->temp.accept(rel_score / ctx->score);
    ctx->history.back() = accept;

    if (accept) {
      puts("ACCEPT");
      source_cl.faces = n_source_set;
      source_cl.box = n_source_box;
      target_cl.faces = n_target_set;
      target_cl.box = n_target_box;
    } else {
      diff_score = 0;
    }
    ctx->score += diff_score;
  }
}

void KAtomizer::do_search_b(AtomizeAssignment *assignment,
                            OptimizationContext *ctx) {

  struct SearchBCandidate {
    int a;
    int b;
    int axis;
    bool dir;
    std::vector<int> faces_in_both;
    std::vector<int> faces_in_a_only;
    std::array<std::vector<double>, 3> spread = {};
  };

  ProbHelperBuilder<SearchBCandidate> selector;
  REP (a, K) {
    const auto &ca = assignment->cluster_data[a];
    if (ca.box.almost_empty()) continue;
    REP (b, K) {
      if (a == b) continue;
      // transfer from a to b
      const auto &cb = assignment->cluster_data[b];
      std::vector<int> faces_not_in;
      SearchBCandidate cnd;
      cnd.a = a;
      cnd.b = b;

      for (auto &face : ca.faces) {
        if (cb.box.in(this->get_points({ face })))
          cnd.faces_in_both.push_back(face);
        else
          cnd.faces_in_a_only.push_back(face);
      }
      if (cnd.faces_in_both.empty()) continue;

      PointVec pts = get_points(cnd.faces_in_both);
      std::array<double, 3> axis_scores = { 0. };
      for (auto &pt : pts) {
        REP (j, 3) {
          double dist = ca.box.main_axis_dist_normed(pt, j);
          dist = clamp(dist, -1, 1) * (1 - base_eps);
          OPA_CHECK(std::abs(dist) < 1, dist, pt, ca.box);
          cnd.spread[j].push_back(dist);
          axis_scores[j] += common::Float(dist).atanh().to_double();
        }
      }
      OPA_DISP0(axis_scores);

      PointVec pts_not = get_points(cnd.faces_in_both);
      std::array<double, 3> axis_scores_not = { 0. };
      std::array<std::vector<double>, 3> spread_not = {};
      for (auto &pt : pts_not) {
        REP (j, 3) {
          double dist = ca.box.main_axis_dist_normed(pt, j) * 0.99;
          OPA_CHECK(std::abs(dist) < 1, pt, ca.box.str(), ca.box.dist2(pt),
                    ca.box.in(pt));
          spread_not[j].push_back(dist);
          axis_scores_not[j] += common::Float(dist).atanh().to_double();
        }
      }

      REP (j, 3) {
        REP (dir, 2) {
          cnd.axis = j;
          cnd.dir = dir;
          double v = axis_scores[j] * OPA_BITSIGN(dir);
          if (v < 0) continue;
          REP (k, 2) {
            v *=
              std::atan(1. / (1e-4 + std::abs(axis_scores_not[(j + k) % 3])));
          }
          OPA_DISP("Pushing candidate ", cnd.a, cnd.b, cnd.axis, cnd.dir, v,
                   axis_scores, axis_scores_not, spread_not);
          selector.push(v, cnd);
        }
      }
    }
  }
  if (!selector.can()) return;
  SearchBCandidate cnd = selector.get(prob_helper);
  const auto &axis_spread = cnd.spread[cnd.axis];
  double cutoff = axis_spread[prob_helper.random_id(axis_spread.size())] *
                  OPA_BITSIGN(cnd.dir);
  auto &ca = assignment->cluster_data[cnd.a];
  auto &cb = assignment->cluster_data[cnd.b];

  TrSet n_b_set = cb.faces;
  TrSet n_a_set = cnd.faces_in_a_only;

  for (auto &face : cnd.faces_in_both) {
    Pos fp = gravity_center(get_points({ face }));
    double u =
      ca.box.main_axis_dist_normed(fp, cnd.axis) * OPA_BITSIGN(cnd.dir);
    if (u >= cutoff)
      n_b_set.push_back(face);
    else
      n_a_set.push_back(face);
  }

  if (n_a_set.size() > 0) {
    BoxSpec n_a_box = compute_box(n_a_set);
    BoxSpec n_b_box = compute_box(n_b_set);
    double old_score = ca.box.area() + cb.box.area();
    double new_score = n_a_box.area() + n_b_box.area();
    double diff_score = new_score - old_score;
    double rel_score = diff_score / old_score;
    OPA_TRACE("candidate ", n_a_set.size(), n_b_set.size(), old_score,
              new_score);

    bool accept = ctx->temp.accept(rel_score / ctx->score);
    ctx->history.back() = accept;

    if (accept) {
      puts("ACCEPT");
      ca.faces = n_a_set;
      ca.box = n_a_box;
      cb.faces = n_b_set;
      cb.box = n_b_box;
    } else {
      diff_score = 0;
    }
    ctx->score += diff_score;
  }
}

void KAtomizer::do_search_c(AtomizeAssignment *assignment,
                            OptimizationContext *ctx) {

  struct SearchBCandidate {
    int a;
    int b;
    int axis;
    bool dir;
    std::vector<int> faces_in_both;
    std::vector<int> faces_in_a_only;
    std::array<std::vector<double>, 3> spread = {};
  };

  std::vector<std::vector<int> > face_sets;
  for (auto &cl : assignment->cluster_data) face_sets.push_back(cl.faces);
  std::vector<std::pair<pii, pii> > neighbours =
    face_graph->graph->find_set_neighbours(face_sets);
  ;

  std::vector<std::tuple<int, int, int> > transfer_list;
  for (auto &e : neighbours) {
    transfer_list.emplace_back(e.first.first, e.first.second, e.second.second);
    transfer_list.emplace_back(e.second.first, e.second.second, e.first.second);
  }
  opa::utils::make_unique(transfer_list);
  std::shuffle(ALL(transfer_list), rng);
  OPA_DISP("Transfer >> ", transfer_list.size());

  for (auto &cnd : transfer_list) {
    int target_tr, cluster_from, cluster_to;
    std::tie(target_tr, cluster_from, cluster_to) = cnd;
    auto &cl_from = assignment->cluster_data[cluster_from];
    auto &cl_to = assignment->cluster_data[cluster_to];

    BoxSpec n_box_to = expand_best_box(cl_to.box, { target_tr });
    check_box_and_face(n_box_to, {target_tr});
    TrSet n_trs_from;
    TrSet n_trs_to = cl_to.faces;
    for (auto &tr : cl_from.faces) {
      PointVec pts = get_points({ tr });
      if (!n_box_to.in(pts))
        n_trs_from.push_back(tr);
      else
        n_trs_to.push_back(tr);
    }
    OPA_DISP0(n_trs_from.size(), n_trs_to.size(), cl_from.faces.size(),
              cl_to.faces.size(), n_box_to.str(), cl_to.box.str(),
              cl_from.box.str(), cluster_from, cluster_to);
    OPA_DISP0(cl_to.box.area(), n_box_to.area());
    if (n_trs_from.size() < cl_from.faces.size() * 2 / 3) continue;
    BoxSpec n_box_from = compute_box(n_trs_from);

    double old_score = cl_to.box.area() + cl_from.box.area();
    double new_score = n_box_to.area() + n_box_from.area();
    double diff_score = new_score - old_score;
    bool accept = ctx->temp.accept(diff_score / ctx->score);
    OPA_DISP0("Got candidate ", old_score, new_score, n_trs_from.size(),
              n_trs_to.size());
    ctx->history.back() = accept;

    if (accept) {
      puts("ACCEPT");
      check_box_and_face(n_box_to, n_trs_to);
      check_box_and_face(n_box_from, n_trs_from);
      cl_to.faces = n_trs_to;
      cl_to.box = n_box_to;
      cl_from.faces = n_trs_from;
      cl_from.box = n_box_from;
      ctx->score += diff_score;
      return;
    } else {
      diff_score = 0;
    }
  }
}

AtomizeAssignment
KAtomizer::search(const KAtomizer_GlobalSearchParams &params) {
  AtomizeAssignment res;
  return res;
}

BoxSpec KAtomizer::expand_best_box(const BoxSpec &box,
                                   const TrSet &tr_ids) const {
  BoxSpec res;
  PointVec pts = get_points(tr_ids);
  utils::container_extend(pts, box.corners());
  return compute_best_box(pts).expand(1.01);
}

BoxSpec KAtomizer::reduce_best_box(const BoxSpec &box, const TrSet &orig_trs,
                                   const TrSet &remove_trs) const {
  BoxSpec res;
  std::set<int> to_remove(ALL(remove_trs));
  TrSet ntrs;
  for (auto &tr : orig_trs) {
    if (!to_remove.count(tr)) ntrs.push_back(tr);
  }
  return compute_box(ntrs);
}
}
OPA_NAMESPACE_DECL3_END
