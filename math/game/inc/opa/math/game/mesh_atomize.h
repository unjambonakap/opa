#pragma once

#include <opa/math/common/rng.h>
#include <opa/math/game/base.h>
#include <opa/math/game/conf.h>
#include <opa/math/game/mesh.h>
#include <opa/math/game/mesh_util.h>
#include <opa/math/game/point_cloud.h>
#include <opa/math/game/quat.h>
#include <opa/or/adaptative_search.h>

OPA_NAMESPACE_DECL3(opa, math, game)
namespace atom {

class RngHelper {
public:
  std::mt19937 rng;
  std::uniform_real_distribution<> unif_rng{ 0., 1. };
  RngHelper(int seed = 0) {
    if (seed == 0) rng.seed(opa::math::common::rng());
  }

  double unif() { return unif_rng(rng); }
};

double linearize(double x, double xl, double xh, double yl, double yh) {
  if (x<xl) return yl;
  if (x>xh) return yh;
  return (x-xl) * (yh - yl) / (xh - xl) + yl;
}

class ProbHelper {
public:
  RngHelper rng_helper;
  ProbHelper() {}

  bool select_one(double score) { return rng_helper.unif() < score; }

  int select_random_id_double_topr(const std::vector<double> &tb, double ratio) {
    std::vector<double> tmp = tb;
    std::sort(ALL(tmp));
    double mx = tmp[ratio * tb.size()];
    double tot = 0;
    for (auto &x : tb) {
      if (x >= mx) tot += x;
    }

    double rand_num = rng_helper.unif() * tot;
    OPA_DISP0(tot, rand_num);
    REP (i, tb.size()) {
      if (tb[i] < mx) continue;
      rand_num -= tb[i];
      if (rand_num <= 0) return i;
    }
    OPA_CHECK(false, tb);
    return -1;
  }

  int select_random_id_double(const std::vector<double> &tb) {
    double tot = 0;
    for (auto &x : tb) {
      tot += x;
    }

    double rand_num = rng_helper.unif() * tot;
    OPA_DISP0(tot, rand_num);
    REP (i, tb.size()) {
      rand_num -= tb[i];
      if (rand_num <= 0) return i;
    }
    OPA_CHECK(false, tb);
    return -1;
  }

  int random_id(int n) { return rng_helper.rng() % n; }

  template <class A, double A::*get_field>
  int select_random_id(const std::vector<A> &tb) {
    std::vector<double> vals;
    for (auto &x : tb) vals.push_back(x.*get_field);
    return select_random_id_double(vals);
  }

  template <class A, double A::*get_field>
  A select_random(const std::vector<A> &tb) {
    return tb[select_random_id<A, get_field>(tb)];
  }
};

template <class K> class ProbHelperBuilder {
public:
  typedef std::pair<double, K> DataType;
  std::vector<DataType> data;

  template <typename... Args> void push(double v, const Args &... args) {
    data.emplace_back(v, K(args...));
  }

  K get(ProbHelper &helper) const {
    return helper.select_random<DataType, &DataType::first>(data).second;
  }
  bool can() const { return data.size() > 0; }
};

typedef std::vector<int> TrSet;

struct ClusterData {
  TrSet faces;
  BoxSpec box;
  algo::FastGraph::BfsCtx bfs_ctx;
};

struct CandidateStruct {
  int face = -1;
  int dest_cluster = -1;

  double score = 1.;
  double box_like_score = 1.;
  double area_change_ratio = 1.;
  double all_box_score = 1.;
  double penalty = 1.;

  int priority = -1;
  bool operator<(const CandidateStruct &other)const {
    return  score < other.score;
  }

  OPA_DECL_COUT_OPERATOR2(CandidateStruct, a.face, a.dest_cluster, a.score,
                          a.box_like_score, a.area_change_ratio,
                          a.all_box_score, a.priority, a.penalty);
};

struct RoundData {
  Rot rot;
  BoxAASpec box;
};

struct InitRoundDebug {
  std::vector<ClusterData> cluster_data;
  std::vector<CandidateStruct> candidates;
  CandidateStruct selected_candidate;
};

struct GlobalRoundData {
  double max_box_area_round;
  std::vector<RoundData> rd;
};

struct AssignmentDebug {
  std::vector<InitRoundDebug> init_round_debug;
};

struct AtomizeAssignment {
  std::vector<ClusterData> cluster_data;
  AssignmentDebug debug_data;
};

struct KAtomizer_SearchParams {
  int nstep;
  opa::OR::SearchStoppingCriteria *stop_criteria = nullptr;
  double T0;
  int state = 0;
  bool force_state = false;
};

struct KAtomizer_GlobalSearchParams {
  int nstep;
  opa::OR::SearchStoppingCriteria *stop_criteria = nullptr;
};

class Temperature {
public:
  RngHelper rng;
  double accept(double v) {
    if (v < 0) return true;
    double prob = std::exp(-v / T);
    return (rng.unif() <= prob);
  }
  void next() { T *= decay; }

  double T = 1.;
  double decay = 0.99;
};

struct OptimizationContext {
  Temperature temp;
  SPTR(opa::OR::SearchStoppingCriteria) stop_criteria;
  double score;
  std::vector<int> history;
};

class KAtomizer {
public:
  KAtomizer(const Mesh &mesh, int K);
  BoxSpec compute_box(const TrSet &tr_ids) const;
  BoxSpec expand_best_box(const BoxSpec &box, const TrSet &tr_ids) const;
  BoxSpec reduce_best_box(const BoxSpec &box, const TrSet &orig_trs, const TrSet &remove_trs) const;

  PointVec get_points(const TrSet &tr_ids) const;
  PointVec &add_points(PointVec &res, const TrSet &tr_ids) const;
  void add_face_to_cluster(ClusterData *cluster, int face_pos);
  AtomizeAssignment initial_assignment();
  double get_assignment_extra_area(const AtomizeAssignment &assignment) const;
  double get_cluster_covered_area(const ClusterData &cluster) const;
  double get_assignment_cost(const AtomizeAssignment &assignment) const;
  bool is_cluster_valid(const ClusterData &cluster, int mode) const;
  bool check_box_and_face(const BoxSpec &box, const TrSet &faces) const;
  bool is_assignment_valid(const AtomizeAssignment &assignment, int mode) const;
  void setup_candidate(const GlobalRoundData &grd, const ClusterData &cl,
                       CandidateStruct &cnd);

  void improve_assignment(AtomizeAssignment *sol,
                          const KAtomizer_SearchParams &params);

  AtomizeAssignment search(const KAtomizer_GlobalSearchParams &params);

  AtomizeAssignment
  create_assignment_from_faces(const std::vector<TrSet> &repartition) {
    AtomizeAssignment res;
    res.cluster_data.resize(K);
    REP (i, K) {
      res.cluster_data[i].faces = repartition[i];
      res.cluster_data[i].box = compute_box(repartition[i]);
    }
    return res;
  }

  void do_search_a(AtomizeAssignment *assignment, OptimizationContext *ctx);
  void do_search_b(AtomizeAssignment *assignment, OptimizationContext *ctx);
  void do_search_c(AtomizeAssignment *assignment, OptimizationContext *ctx);

  SPTR(FaceGraph2) face_graph;
  Mesh mesh;
  double mesh_area;
  std::vector<Triangle> triangles;
  int K;
  ProbHelper prob_helper;
  bool debug = false;
};

std::vector<BoxSpec> atomize_mesh_K(const Mesh &mesh, int K);
std::vector<BoxSpec> atomize_mesh(const Mesh &mesh);
}

OPA_NAMESPACE_DECL3_END

/*
   Fitting mesh with simpler curves (atoms). Atoms can be:
    - cylinder ? might be expansive for intersection test. Area gain not so
  great?
    - bounding box
    - sphere

  Want to minimize a function of:
    - number of atoms / cost per atom type
    - area overestimation -> area mesh / area atoms
    - ideally should represent cost of intersection computation / distance
  finding in a given direction

  x Need to cover all triangles -> give a set of points to cover.
  x primitives to compute intersections between mesh - atom


  Can split triangles if too large -> remeshing algorithm
  Optimization algorithm: K partition of triangles.

  initial configuration?
    - random assignment
    - RAX creation:
      * select random first face.  prob of selecting other face for parition
  tied to distance to other partitions
      * creating assigment:
        * greedy assignment
            -> each face to one which gets smallest size increase
            -> keep all partitions of equal size -> min(max(volume))
        * grasp like: prob of choice for each

  Improving:
  simulated annearling
    - prob of swap depends on:
      - new cost for the two atoms
        - volume of the atoms -> maybe hyperbolic
        - how tight the representation is -> bound area vs real area
        - how cube like
      - how recently touched the triagnle was


  Need to compute:
    - DONE: sphere from a set of points
    - DONE: one of
      * worse/average AABB for set of points -> low res sampling on rotation
  half ball
      * best bounding box
    - TODO: area of mesh intersected with atom -> not possible here since
  triangles are used



  -> find atom position/type/size


   */
