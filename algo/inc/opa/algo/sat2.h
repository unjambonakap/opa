#pragma once

#include <glib/gtl/map_util.h>
#include <opa/algo/graph.h>
#include <opa/algo/graph_util.h>
#include <opa/utils/hash.h>
#include <opa_common.h>
#include <type_traits>

OPA_NAMESPACE_DECL2(opa, algo)
constexpr int BAD_NODE = -1;
template <class Key> class Sat2 {
public:
  class Term {
  public:
    enum OpType {
      OR,
      AND,
      SUB,
      XOR,
      ADD,
      MUL,
      EQ,
    };

    Term(Sat2<Key> *sat) { m_sat = sat; }
    void set(const Key &k) {
      literals.push_back(m_sat->get_id(k));
      vals = { 0, 1 };
    }

    void set(bool v) { vals = { v }; }
    void set(s8 v) { vals = { v }; }

    void set(const Term &a, const Term &b, OpType op_type) {
      OPA_CHECK(a.literals.size() + b.literals.size() <= 2,
                "only solving 2 sat");
      literals = a.literals;
      literals.insert(literals.end(), ALL(b.literals));
      REP (i, 1 << literals.size()) {
        std::map<int, int> mp;
        REP (j, literals.size())
          mp[literals[j]] = i >> j & 1;
        int v1 = a.get_val(mp);
        int v2 = b.get_val(mp);
        vals.push_back(eval(op_type, v1, v2));
      }
      OPA_DISP0(a.vals, b.vals, vals);
      simplify();
    }

    int get_val(const std::map<int, int> &mp) const {
      int entry = 0;
      REP (i, literals.size())
        entry |= mp.find(literals[i])->second << i;
      return vals[entry];
    }

    void simplify() {
      if (literals.size() == 2) {
        REP (i, 2) {
          int other = 1 << (i ^ 1);
          int cur = 1 << i;
          if (vals[0] == vals[cur] && vals[other] == vals[other | cur]) {
            literals = { literals[other] };
            vals = { vals[0], vals[other] };
            break;
          }
        }
      }

      if (literals.size() == 1) {
        if (vals[0] == vals[1]) {
          literals = {};
          vals = { vals[0] };
        }
      }
    }

    template <typename... Args> Term create(const Args &... args) const {
      Term res(m_sat);
      res.set(args...);
      return res;
    }

    void create_op(const Term &a, const Term &b, OpType op_type) {}

    Term operator!() const {
      Term res = *this;
      for (auto &v : res.vals) v = !v;
      return res;
    }

#define FORWARD_OP(op)                                                         \
  template <typename... Args> Term operator op(const Args &... args) const {   \
    Term other = create(args...);                                              \
    return other op * this;                                                    \
  }
    FORWARD_OP(||);
    FORWARD_OP(&&);
    FORWARD_OP (^);
    FORWARD_OP(+);
    FORWARD_OP(*);
    FORWARD_OP(==);

#define DECL_OP(op, typ)                                                       \
  Term operator op(const Term &a) const {                                      \
    Term res(m_sat);                                                           \
    res.set(*this, a, OpType::typ);                                            \
    return res;                                                                \
  }

    DECL_OP(||, OR);
    DECL_OP(&&, AND);
    DECL_OP (^, XOR);
    DECL_OP(+, ADD);
    DECL_OP(-, SUB);
    DECL_OP(*, MUL);
    DECL_OP(==, EQ);

    int eval(OpType op_type, int v1, int v2) const {
      switch (op_type) {
      case OpType::AND:
        return v1 && v2;
      case OpType::OR:
        return v1 || v2;
      case OpType::XOR:
        return v1 ^ v2;
      case OpType::ADD:
        return v1 + v2;
      case OpType::SUB:
        return v1 - v2;
      case OpType::EQ:
        return v1 == v2;
      case OpType::MUL:
        return v1 * v2;
      }
      OPA_CHECK0(false);
    }

    std::vector<pii> get_edges() const {
      std::vector<pii> edges;
      if (literals.size() == 0) {
        if (vals[0] == 0) edges.emplace_back(BAD_NODE, BAD_NODE);
      } else if (literals.size() == 1) {
        OPA_CHECK0(vals[1] == 0 || vals[1] == 1);
        edges.emplace_back(literals[0] ^ vals[1], literals[0] ^ vals[0]);
      } else if (literals.size() == 2) {
        OPA_DISP0(literals, vals);
        REP (i, 4) {
          if (vals[i]) continue;
          int da = i & 1;
          int db = i >> 1;
          edges.emplace_back(literals[0] ^ da ^ 1, literals[1] ^ db);
          edges.emplace_back(literals[1] ^ db ^ 1, literals[0] ^ da);
        }
      }
      return edges;
    }

    OPA_DECL_COUT_OPERATOR2(Term, a.vals, a.literals);

    std::vector<int> vals;
    std::vector<int> literals;

    Sat2<Key> *m_sat;
  };

  template <typename... Args> Term term(const Args &... args) {
    Term res = Term(this);
    res.set(Key(args...));
    return res;
  }

  Sat2<Key> &operator+=(const Term &term) {
    terms.push_back(term);
    return *this;
  }

  bool compute() {
    int n = m_rmp.size();
    // Edge semantic:  source => dest
    graph.reset(2 * n, Mode::DIGRAPH);
    for (auto &term : terms) {
      auto edges = term.get_edges();
      for (auto &edge : edges) {
        if (edge.first == BAD_NODE) return false;
        OPA_DISP0(edge);
        graph.adde(edge.first, edge.second, false);
      }
    }

    auto ccs = compute_digraph_connected_components(graph);
    std::vector<bool> vis_cc(res.size(), false);
    std::vector<bool> isset;
    isset.resize(2 * n, false);
    res.resize(n, false);

    for (const auto &cc : ccs) {
      OPA_DISP0(cc.nodes);
      if (isset[cc.nodes[0] >> 1]) continue;
      for (const auto &x : cc.nodes) {
        if (isset[x >> 1]) return false;
        isset[x >> 1] = 1;
        res[x >> 1] = x & 1;
      }
    }
    return true;
  }

  int get_id(const Key &k) { return m_rmp.get(k) * 2 + 1; }
  bool get(const Key &k) const { return res[m_rmp.get_or_die(k)]; }

  template <typename... Args> bool get2(const Args &... args) const {
    return res[m_rmp.get_or_die(Key(args...))];
  }

  std::vector<bool> res;

  FastGraph graph;
  std::vector<Term> terms;
  OPA_ACCESSOR_R(utils::Remapper<Key>, m_rmp, rmp);
  utils::Remapper<Key> m_rmp;
};

OPA_NAMESPACE_DECL2_END
