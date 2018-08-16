#pragma once
#include <opa_common.h>
#include <opa/crypto/la/base.h>

namespace asmjit {
class X86GpVar;
}

OPA_NAMESPACE(opa, crypto, la)
class CipherEdge;

class CipherNode : public opa::utils::Initable {
public:
  enum class Type : int { Input, Key, Mid, Output };

  virtual void init(CipherBlock *block) {
    opa::utils::Initable::init();
    this->block = block;
  }
  bool is_input() const { return type == Type::Input; }
  bool is_output() const { return type == Type::Output; }
  bool is_key() const { return type == Type::Key; }
  bool is_mid() const { return type == Type::Mid; }

  CipherBlock *block;
  BitVec iv, ov;
  int rem_in;
  Type type = Type::Mid;
  std::vector<CipherEdge *> in;
  std::vector<CipherEdge *> out;
};

class CipherEdge : public opa::utils::Initable {
public:
  virtual void init(CipherNode *src, CipherNode *dest, u32 src_pos,
                    u32 dest_pos) {
    opa::utils::Initable::init();
    this->src = src;
    this->src_pos = src_pos;
    this->dest = dest;
    this->dest_pos = dest_pos;
  }

  CipherNode *src, *dest;
  u32 src_pos, dest_pos;
};

class IdBlock : public CipherBlock {
public:
  virtual void init(u32 size) {
    CipherBlock::init(Params()
                        .call(CallInfo(size, 0, true))
                        .input_size(size)
                        .output_size(size)
                        .fast_eval(true));
    desc() = opa::utils::SPrintf("IdBlock sz=%d", size);
  }
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    ov = iv;
  }

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;
};

class InputBlock : public CipherBlock {
public:
  virtual void init(u32 size) {
    CipherBlock::init(Params().call(CallInfo(size, 0, true)).output_size(size));
    desc() = opa::utils::SPrintf("InputBlock sz=%d", size);
  }
  virtual void do_get_relations(Relations &rels) const override {}
  virtual void do_get_relations_diff(Relations &rels) const override {}
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    ov = val;
  }

  virtual void setup_jit(const JitBuilder &builder) const override {}

  mutable BitVec val;
};

class OutputBlock : public CipherBlock {
public:
  virtual void init(u32 size) {
    CipherBlock::init(Params().call(CallInfo(size, 0, true)).input_size(size));
    desc() = opa::utils::SPrintf("OutputBlock sz=%d", size);
  }
  virtual void do_get_relations(Relations &rels) const override {}
  virtual void do_get_relations_diff(Relations &rels) const override {}
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override {
    val = iv;
  }
  virtual void setup_jit(const JitBuilder &builder) const override {}
  mutable BitVec val;
};

typedef InputBlock KeyBlock;

class CipherGraph : public CipherBlock {
public:
  struct PlugDesc {
    PlugDesc(CipherNode *node, int off = 0, int sz = -1) {
      this->node = node;
      this->off = off;
      this->sz = sz;
    }

    int osz() const {
      if (sz != -1)
        return sz;
      return node->block->output_size() - off;
    }

    int isz() const {
      if (sz != -1)
        return sz;
      return node->block->input_size() - off;
    }

    CipherNode *node;
    int off;
    int sz;
  };

  struct RelationRef {
    int node_id;
    Relation rel;
    RelationRef(int node_id, const Relation &rel) {
      this->node_id = node_id;
      this->rel = rel;
    }

    friend std::ostream &operator<<(std::ostream &os, const RelationRef &r) {
      os << RAW_OPA_DISP_VARS(r.node_id, r.rel);
      return os;
    }
  };

  struct InputData {
    CipherNode *node;
    int pos; // for key type, it starts at raw_input_size()
    int par = -1;

    CipherNode::Type type;
    std::set<int> out;
    int prev = -1;
    bool is_input() const { return type == CipherNode::Type::Input; }
    bool is_output() const { return type == CipherNode::Type::Output; }
    bool is_key() const { return type == CipherNode::Type::Key; }
    bool is_mid() const { return type == CipherNode::Type::Mid; }
  };

  struct GraphRel {
    RangeCoverage in, out, mid, key;
    double cost = 0;
    std::vector<RelationRef> rels;

    void add(const InputData &x) {
      if (x.is_input())
        in.add(x.par);
      else if (x.is_output())
        out.add(x.par);
      else if (x.is_key())
        key.add(x.par);
      else if (x.is_mid())
        mid.add(x.par);
    }

    void merge(const GraphRel &peer) {
      cost *= peer.cost;
      mid.do_xor(peer.mid);
      key.do_xor(peer.key);
      in.do_xor(peer.in);
      out.do_xor(peer.out);
    }

    bool operator<(const GraphRel &peer) const {
      OPA_LT_OP(peer, cost, in, out, mid, key);
    }
  };

  struct GraphRelStorer {

    std::set<GraphRel> rels;

    bool check(const GraphRel &rel);
    void add(const GraphRel &rel);
    std::vector<GraphRel> extract(const RangeCoverage &repr);
  };

  CipherNode *add_node() {
    nodes.emplace_back(new CipherNode);
    return nodes.back().get();
  }

  CipherNode *add_node(CipherBlock *block) {
    CipherNode *res = new CipherNode;
    res->init(block);
    nodes.emplace_back(res);
    return res;
  }

  void add_range_edges(CipherNode *src, CipherNode *dest, u32 src_pos,
                       u32 dest_pos, u32 len) {
    REP (i, len)
      add_edge()->init(src, dest, src_pos + i, dest_pos + i);
  }

  CipherEdge *add_edge() {
    edges.emplace_back(new CipherEdge);
    return edges.back().get();
  }

  virtual void init(const Params &params) override;
  virtual void fini() override;
  virtual void build() override;
  virtual void verify() const override;
  virtual void do_evaluate(const BitVec &iv, const BitVec &kv,
                           BitVec &ov) const override;

  void smart_plug(const std::vector<PlugDesc> &ins, PlugDesc out);

  std::string build_graph_desc() const { return "graphDesc"; }

  void init_rels();

  virtual void do_get_relations(Relations &rels) const override;
  virtual void do_get_relations_diff(Relations &rels) const override;

  void do_rel_rec(GraphRelStorer &res, const std::vector<GraphRel> &candidates,
                  const RangeCoverage &targets) const;

  std::pair<RangeCoverage, RangeCoverage> rel_to_io(const GraphRel &rel) const;

  virtual void do_get_pre_fail(const Basis &out, Basis &res) const override;

  virtual void setup_jit(const JitBuilder &builder) const override;

  std::unique_ptr<asmjit::X86GpVar>
  get_var_for_remap(std::map<const CipherNode *, JitBuilder> &mp,
                    JitContext &jc, const JitBuilder &builder, CipherNode *to,
                    int to_pos, CipherNode *from, int from_pos, int len,
                    bool is_key) const;

  struct DiffRelEntry {
    // diff basis, [ B | M ]
    // M diff relation basis, B=M[mask_elems]
    opa::math::common::Matrix<u32> basis;
    // bit pos that are to 1 in mask&sum_mask
    std::vector<u32> mask_elems;
    std::vector<BitVec> diffs;
    const CipherNode *node;
    BitVec mask;
    BitVec sum_mask;
    int basis_size;

    opa::math::common::Matrix<u32> val_to_vec(const BitVec &value) const;
  };

protected:
  struct DiffStatus {
    DiffStatus(int pos, const BitVec &cur, bool has_non_basis)
        : pos(pos), cur(cur), has_non_basis(has_non_basis) {}
    int pos;
    BitVec cur;
    bool has_non_basis;
  };

  std::vector<u32> relation_to_io_space(const CipherNode *node,
                                        const Relation &rel) const;
  u32 node_out_to_io_space(const CipherNode *node, int output_pos) const;
  u32 node_in_to_io_space(const CipherNode *node, int input_pos) const;

  struct DiffData {
    std::vector<DiffRelEntry> rels;
    std::vector<int> proc_order;
    std::vector<BitVec> res_vec;
    Relations *res_rels;
  };

  virtual void init_graph() = 0;
  void compute_collapse();
  void build_graph_rels(std::vector<GraphRel> &rels) const;
  void diff_rec(const DiffStatus &status) const;

  InputBlock input_block;
  IdBlock output_block;
  KeyBlock key_block;
  std::vector<CipherNode *> order;
  std::map<CipherNode *, int> node_order;

  CipherNode *input;
  CipherNode *output;
  CipherNode *key;

  std::vector<UPTR(CipherNode)> nodes;
  std::vector<UPTR(CipherEdge)> edges;
  bool is_built = false;
  bool m_init_rels = false;

  // Store relations for each nodes, either linear or diff
  mutable std::vector<GraphRel> tb_linear;
  std::map<std::pair<const CipherNode *, int>, int> oid;
  std::vector<InputData> data;
  mutable DiffData m_diff_data;
};

OPA_NAMESPACE_END(opa, crypto, la)
