#include <opa_common.h>
#include <opa/crypto/la/graph.h>
#include <opa/crypto/la/relation_driver.h>
#include "jit.h"

using namespace std;
using namespace opa::math::common;

OPA_NAMESPACE(opa, crypto, la)

void IdBlock::setup_jit(const JitBuilder &builder) const {
  auto c = builder.c;
  OPA_CHECK0(builder.input_vars.size() == 1);
  c->mov(*builder.output_var, *builder.input_vars[0]);
}

void IdBlock::do_get_relations_diff(Relations &rels) const {
  REP (i, input_size()) {
    BitVec v(input_size() + output_size());
    v.set(i, 1);
    v.set(input_size() + i, 1);
    rels.add_basis_vec(v);
  }
}

void IdBlock::do_get_relations(Relations &rels) const {
  REP (i, input_size()) {
    Relation r;
    r.out.add(i);
    r.in.add(i);
    r.cost.set_bias(1);
    rels.add_rel(r);
  }
}

bool CipherGraph::GraphRelStorer::check(const GraphRel &rel) {
  if (rel.in.size() > 5)
    return false;
  if (rel.out.size() > 5)
    return false;
  if (rel.cost < 0.1)
    return false;
  return true;
}

void CipherGraph::GraphRelStorer::add(const GraphRel &rel) {
  if (!check(rel))
    return;

  int ms = 0xa0;

  auto it = rels.find(rel);
  if (it != rels.end() && it->cost > rel.cost)
    ;
  else {
    if (it != rels.end())
      rels.erase(it);
    rels.insert(rel);
  }
  if (rels.size() > ms) {
    rels.erase(rels.begin());
  }
}

std::vector<CipherGraph::GraphRel>
CipherGraph::GraphRelStorer::extract(const RangeCoverage &repr) {
  std::vector<GraphRel> res;

  for (auto it = rels.begin(); it != rels.end();) {
    auto cur = it++;
    if (cur->mid.has_one(repr)) {
      res.pb(*cur);
      rels.erase(cur);
    }
  }
  return res;
}

void CipherGraph::fini() {

  order.clear();
  node_order.clear();
  nodes.clear();
  edges.clear();
  tb_linear.clear();
  oid.clear();
  data.clear();
  CipherBlock::fini();
}

void CipherGraph::init(const Params &params) {

  is_built = false;
  m_init_rels = false;
  CipherBlock::init(params);
  input_block.init(raw_input_size());
  output_block.init(output_size());
  key_block.init(key_size());

  input = add_node();
  input->type = CipherNode::Type::Input;
  output = add_node();
  output->type = CipherNode::Type::Output;
  key = add_node();
  key->type = CipherNode::Type::Key;

  input->init(&input_block);
  output->init(&output_block);
  key->init(&key_block);
  desc() = opa::utils::SPrintf("(input=%d,output=%d,key=%d)", raw_input_size(),
                               output_size(), key_size());
  init_graph();
}

void CipherGraph::build() {
  if (is_built)
    return;

  for (auto &node : nodes) {
    node->check_init();
    node->block->build();
  }

  for (auto &e : edges) {
    e->dest->in.push_back(e.get());
    e->src->out.push_back(e.get());
  }

  for (auto &node : nodes) {
    struct Sorter {
      bool operator()(const CipherEdge *a, const CipherEdge *b) const {
        return a->dest_pos < b->dest_pos;
      }
    };
    sort(ALL(node->in), Sorter());
  }

  std::queue<CipherNode *> q;
  for (auto &node : nodes) {
    node->iv.init(node->block->input_size());
    node->ov.init(node->block->output_size());
    node->rem_in = node->in.size();
    if (node->rem_in == 0)
      q.push(node.get());
  }

  while (!q.empty()) {
    CipherNode *node = q.front();
    q.pop();
    node_order[node] = order.size();
    order.push_back(node);

    for (auto &e : node->out) {
      if (!--e->dest->rem_in)
        q.push(e->dest);
    }
  }

  init_rels();
  is_built = true;
}

void CipherGraph::verify() const {
  OPA_CHECK0(is_built);
  for (const auto &node : nodes) {
    node->block->verify();

    bool check_in = true;
    bool check_out = true;

    if (node.get() == input)
      check_in = false;
    if (node.get() == key)
      check_in = false;
    if (node.get() == output)
      check_out = false;

    if (check_in) {
      RangeCoverage range(0, node->block->input_size());
      for (const auto &e : node->in) {
        range.add(e->dest_pos);
      }
      OPA_CHECK0(range.covered() && range.unique());
    }

    if (check_out) {
      RangeCoverage range(0, node->block->output_size());
      for (const auto &e : node->out) {
        range.add(e->src_pos);
      }
      OPA_DISP("FOR NODE >> ", node->block->desc(), range, node.get() == input,
               node.get() == key, node.get() == output, uintptr_t(node.get()),
               uintptr_t(input));
      OPA_CHECK0(range.covered());
    }
  }
}
void CipherGraph::do_evaluate(const BitVec &iv, const BitVec &kv,
                              BitVec &ov) const {
  OPA_CHECK0(is_built);

  input_block.val = iv;
  key_block.val = kv;

  for (auto &node : order) {
    if (context())
      context()->start_evaluate(node);

    node->ov = node->block->evaluate(node->iv);

    if (context()->is_debug()) {
      OPA_DISP("Evaluating node ", node->block->desc(), node->iv.to<u64>(),
               node->ov.to<u64>());
    }

    if (context())
      context()->end_evaluate(node);

    for (auto &e : node->out) {
      e->dest->iv.set(e->dest_pos, node->ov.get(e->src_pos));
    }
  }

  ov = output->ov;
}

void CipherGraph::init_rels() {
  if (m_init_rels)
    return;
  m_init_rels = true;

  for (int i = 0; i < order.size(); ++i) {
    auto node = order[i];
    REP (j, node->block->output_size()) {
      int cur = oid.size();
      oid[MP(node, j)] = cur;
      InputData d;

      d.node = node;
      d.pos = j;
      d.type = d.node->type;
      if (d.node->is_output())
        d.out.insert(cur);
      data.pb(d);
    }
  }

  compute_collapse();

  int bad = 0;
  for (auto &d : data) {
    if (d.prev == -1)
      continue;
    ++bad;
    d.par = d.prev;
  }

  REP (i, data.size()) {
    int cur = i;
    while (data[cur].par != cur && data[cur].par != -1) {
      cur = data[cur].par;
    }

    data[i].par = cur;
    if (data[i].node->is_output())
      data[cur].type = CipherNode::Type::Output,
      data[cur].out.insert(ALL(data[i].out));
  }
}

void CipherGraph::compute_collapse() {
  for (auto &d : data) {
    d.prev = -1;
    if (d.node->is_key() || d.node->is_input())
      continue;
    std::vector<RangeCoverage> cover =
      d.node->block->get_pre(RangeCoverage({ d.pos }));
    if (cover.size() != 1 || cover[0].size() != 1)
      continue;
    int x = cover[0].to_vec()[0];
    CipherEdge *e = d.node->in[x];
    d.prev = oid[MP(e->src, e->src_pos)];
  }
}

u32 CipherGraph::node_in_to_io_space(const CipherNode *node,
                                     int input_pos) const {
  CipherEdge *e = node->in[input_pos];
  return oid.find(MP(e->src, e->src_pos))->second;
}

u32 CipherGraph::node_out_to_io_space(const CipherNode *node,
                                      int output_pos) const {
  return oid.find(MP(node, output_pos))->second;
}

std::vector<u32> CipherGraph::relation_to_io_space(const CipherNode *node,
                                                   const Relation &rel) const {
  vector<u32> res;
  for (auto x : rel.out.all())
    res.pb(node_out_to_io_space(node, x));

  for (auto x : rel.in.all())
    res.pb(node_in_to_io_space(node, x));
  return res;
}

void CipherGraph::build_graph_rels(std::vector<GraphRel> &rels) const {
  for (int i = 0; i < order.size(); ++i) {
    auto node = order[i];
    const auto &lst = order[i]->block->get_relations();

    for (const auto &r : lst.tb) {
      GraphRel a;
      a.cost = r.cost.prob();
      a.rels.emplace_back(i, r);

      a.mid.add_lst(relation_to_io_space(node, r));
      // OPA_DISP("add rel ", in, out);
      rels.pb(a);

      // OPA_CHECK(out.size() > 0, r.out.all());
    }
  }
}

struct CmpRels {
  bool operator()(const CipherGraph::DiffRelEntry &a,
                  const CipherGraph::DiffRelEntry &b) const {
    return get_cost(a) < get_cost(b);
  }

  u64 get_cost(const CipherGraph::DiffRelEntry &a) const {
    u64 cost = 1;
    cost *= max<int>(1, a.diffs.size());
    cost *= 1ull << min<u32>(63, a.basis_size);
    return cost;
  }
};

void CipherGraph::diff_rec(const DiffStatus &status) const {
  if (status.pos == m_diff_data.rels.size()) {
    if (status.cur.isz())
      return;
    m_diff_data.res_vec.pb(status.cur);
    return;
  }
  DiffStatus nstatus = status;
  nstatus.pos += 1;

  const DiffRelEntry &target = m_diff_data.rels[status.pos];
   OPA_DISP("diff rec >> ", status.pos, status.cur, target.mask,
           target.diffs.size(), status.has_non_basis);
  REP (i, target.diffs.size() + 1) {
    DiffStatus next = nstatus;
    BitVec target_bits;

    bool noop = i == target.diffs.size();
    if (noop)
      target_bits.init(data.size());
    else
      target_bits = target.diffs[i];

    if (!status.cur.xorz(target_bits)
           .andz(target.mask)
           .andz(target.sum_mask)
           .isz())
      continue;
    next.cur = status.cur.orz(target_bits);
    next.has_non_basis |= !noop;
    diff_rec(next);
  }

  if (status.has_non_basis && target.basis.getN() > 0) {
    DiffStatus next = nstatus;
    Matrix<u32> vec = target.val_to_vec(status.cur);
    Matrix<u32> base_coord;
    // OPA_DISP0("TRY ", vec);
    // OPA_DISP0("got ", target.basis, target.mask_elems);
    if (target.basis.reduce(vec, target.mask_elems.size(), true, &base_coord)) {
      Matrix<u32> rvec =
        base_coord.mul(target.basis).get_submatrix(0, target.mask_elems.size());
      BitVec v2 = BitVec::FromBinaryVec(rvec.tovec());
      next.cur = next.cur.orz(v2);
      diff_rec(next);
    }
  }
}

opa::math::common::Matrix<u32>
CipherGraph::DiffRelEntry::val_to_vec(const BitVec &value) const {
  opa::math::common::Matrix<u32> res;
  res.initialize(&GF2, 1, basis.getM());
  REP (i, mask_elems.size())
    res(i) = value.get(mask_elems[i]);
  res.set_submatrix(value.to_mat(), 0, mask_elems.size());
  return res;
}

void CipherGraph::do_get_relations_diff(Relations &rels) const {
  if (input_size() <= 16) {
    OPA_DISP("processing ", desc());
    get_relations_walsh_diff(rels);
    return;
  }

  for (auto &node : nodes) {
    OPA_DISP("On node ", node->is_mid(), node->block->desc());
    if (!node->is_mid())
      continue;
    const Relations &node_rels = node->block->get_relations(true);
    m_diff_data.rels.emplace_back();
    DiffRelEntry &entry = m_diff_data.rels.back();
    entry.node = node.get();
    entry.basis_size = node_rels.diff_basis.dim();

    BitVec mask(data.size());
    REP (i, node->block->input_size())
      mask.toggle(node_in_to_io_space(node.get(), i));
    REP (i, node->block->output_size())
      mask.toggle(node_out_to_io_space(node.get(), i));
    entry.mask = mask;

    for (auto &rel : node_rels.tb) {
      bool use_key = false;
      for (auto &v : rel.in.all())
        if (v >= node->block->raw_input_size()) {
          use_key = true;
          break;
        }
      if (use_key)
        continue;
      // TODO: temporary
      if (!rel.cost.is_exact())
        continue;
      BitVec evec =
        BitVec::FromVec(relation_to_io_space(node.get(), rel), data.size());
      OPA_DISP("PUSHING DIFFS >> ", rel.cost.prob(), evec);
      entry.diffs.push_back(evec);
    }
  }

  m_diff_data.res_rels = &rels;
  sort(ALL(m_diff_data.rels), CmpRels());
  BitVec sum_mask(data.size());
  for (auto &e : m_diff_data.rels) {
    e.sum_mask = sum_mask;
    sum_mask = sum_mask.orz(e.mask);
  }

  for (auto &entry : m_diff_data.rels) {
    const Relations &node_rels = entry.node->block->get_relations(true);
    const Matrix<u32> &rel_basis = node_rels.diff_basis.reduced();
    REP (i, data.size())
      if (entry.mask.get(i) && entry.sum_mask.get(i))
        entry.mask_elems.pb(i);
    entry.basis.initialize(&GF2, node_rels.diff_basis.dim(),
                           entry.mask_elems.size() + data.size());
    OPA_DISP0("FUU ", entry.mask_elems);
    int node_is = entry.node->block->input_size();
    int node_os = entry.node->block->output_size();
    REP (i, entry.basis.getN()) {
      BitVec v(data.size());

      OPA_CHECK0(rel_basis.getM() == node_is + node_os);
      REP (j, node_is) {
        if (rel_basis.get(i, j))
          v.toggle(node_in_to_io_space(entry.node, j));
      }
      REP (j, node_os) {
        if (rel_basis.get(i, node_is + j))
          v.toggle(node_out_to_io_space(entry.node, j));
      }
      entry.basis.set_row(i, entry.val_to_vec(v));
    }
    entry.basis.row_echelon(-1, entry.mask_elems.size());
    //OPA_DISP0("pushing entry basis ", entry.basis);
    OPA_DISP0("before ", rel_basis);
  }

  DiffStatus start_status(0, BitVec(data.size()), false);
  diff_rec(start_status);

  for (auto &v : m_diff_data.res_vec) {
    Relation r;
    // TODO: temp
    r.cost.set_bias(1.);

    REP (i, data.size()) {
      if (data[i].is_mid())
        continue;
      // Bit does not belong to differential
      if (!v.get(i))
        continue;

      if (data[i].is_output()) {
        for (auto a : data[i].out) {
          r.out.add(data[a].pos);
        }
      } else {
        r.in.add(data[i].pos);
      }
    }
    m_diff_data.res_rels->add_rel(r);
    OPA_DISP0(v);
  }
  OPA_DISP0(opa::utils::Join("\n", *m_diff_data.res_rels));
}

void CipherGraph::do_get_relations(Relations &rels) const {
  if (input_size() <= 16) {
    get_relations_walsh(rels);
  }
  /*
  puts("HHERE do get graph");

  REP (i, data.size())
    OPA_DISP0(i, int(data[i].type), data[i].pos, data[i].par);
    */

  build_graph_rels(tb_linear);

  GraphRelStorer storer;

  // Store all relations in the graph index space.
  for (auto &r : tb_linear) {
    GraphRel r2;
    std::set<int> seen;
    for (int x : r.mid.all()) {
      int p = data[x].par;
      r2.add(data[p]);
      seen.insert(p);
    }
    if (seen.size() == 1)
      continue;

    r2.rels = r.rels;
    r2.cost = r.cost;

    storer.add(r2);
  }

  for (auto node : order) {
    if (node == input || node == key || node == output)
      continue;

    RangeCoverage repr;

    for (auto &oe : node->out) {
      int x = oid.find(MP(oe->src, oe->src_pos))->ND;
      x = data[x].par;

      if (!data[x].node->is_mid())
        continue;
      repr.add(x);
    }
    if (repr.size() == 0)
      continue;
    auto cur = storer.extract(repr);
    // for (auto &k : cur)
    //  OPA_DISP0(k.mid);

    do_rel_rec(storer, cur, repr);
  }

  std::vector<GraphRel> tmp(ALL(storer.rels));
  REP (i, tmp.size())
    REP (j, i) {
      GraphRel x = tmp[i];
      x.merge(tmp[j]);
      storer.add(x);
    }

  for (auto &x : storer.rels) {
    Relation::RelationCost cost;
    cost.set_bias(x.cost);
    if (!driver()->should_add_rel(cost))
      continue;

    Relation u;
    std::tie(u.in, u.out) = rel_to_io(x);
    if (u.out.size() == 0 || u.in.size() == 0)
      continue;
    // OPA_DISP("GOT REL >> ", u.in, u.out);
    OPA_CHECK0(u.out.size() > 0);
    OPA_CHECK0(u.in.size() > 0);
    u.cost = cost;
    driver()->add_rel(u);
  }
  driver()->set_rels(rels);
}

void CipherGraph::do_get_pre_fail(const Basis &out, Basis &res) const {
  /* the argument is:
     we want to find the input relations that allows us to find every output
     relation (stored in Basis out)
     We maintain at each node in the graph the relations that we must know

     For the next node (in reverse dag order) we do:
     c.get_pre() -> outputs a list of input relations that we must know to
     find out about the output relations
     Propagate these input relations to the next nodes. This is done123

     */
  out.check_init();

  std::map<CipherNode *, Basis> data;

  for (auto node : order) {
    if (node == input || node == key || node == output)
      continue;
    data[node].init(node->block->output_size());
  }

  data[input].init(input_size());
  data[output] = out;
  // OPA_DISP("GET PRE FAIL >> ", out.dim());

  REPV (i, order.size()) {
    CipherNode *node = order[i];
    if (node == input || node == key)
      continue;
    // OPA_DISP("ON NODE ", node->block->desc());
    Basis &cur = data[node];
    cur.check_init();
    Basis res2;
    node->block->do_get_pre(cur, res2);

    std::map<CipherNode *, std::map<int, int> > rmp;
    for (auto &e : node->in) {
      CipherNode *src = e->src;
      int src_pos = e->src_pos;
      if (src == key)
        src = input, src_pos += raw_input_size();
      rmp[src][e->dest_pos] = src_pos;
    }
    // OPA_DISP("AFTER COMPUTING BLOCK ", node->block->desc(), res2.dim(),
    //         res2.mat());
    // OPA_DISP0(rmp);

    REP (i, res2.dim()) {
      for (auto &x : rmp) {
        RangeCoverage tmp;
        for (auto &k : x.ND) {
          if (res2.mat().get(i, k.ST))
            tmp.add(k.ND);
        }
        data[x.ST].add(tmp);
      }
    }
  }
  res = data[input];
  res.reduce();
}

void CipherGraph::setup_jit(const JitBuilder &builder) const {
  std::map<const CipherNode *, JitBuilder> mp;
  JitContext jc;
  jc.c = builder.c;

  for (auto node : order) {
    if (node == input || node == key)
      continue;
    auto x = node->block;

    JitBuilder cur;
    x->create_jit_builder(cur, jc);
    mp[node] = cur;
    std::map<CipherNode *, std::vector<pii> > edge_map;
    for (auto e : node->in)
      edge_map[e->src].pb(MP(e->src_pos, e->dest_pos));
    int input_base_size = x->call().input_base_size;
    for (auto u : cur.input_vars)
      jc.c->xor_(*u, *u);
    if (cur.key_var)
      jc.c->xor_(*cur.key_var, *cur.key_var);

    OPA_CHECK_EQ0(cur.input_vars.size(), x->raw_input_size() / input_base_size);

    for (auto a : edge_map) {
      sort(ALL(a.second));
      REP (i, a.second.size() - 1) {
        OPA_CHECK0(a.second[i] + MP(1, 1) == a.second[i + 1]);
      }

      OPA_CHECK0(a.second.size() > 0);

      int start = a.second[0].ND;
      int end = a.second.back().ND + 1;

      int off = 0;
      int bucket_size = x->call().input_base_size;
      auto other = a.first->block;

      bool is_key = false;
      if (start >= x->raw_input_size()) {
        is_key = true;
        start -= x->raw_input_size();
        end -= x->raw_input_size();
        bucket_size = x->call().key_base_size;
      }

      int diff_tofrom = a.second[0].ST - start;

      for (int target = start / bucket_size; target <= (end - 1) / bucket_size;
           ++target) {

        int start2 = std::max(start, target * bucket_size);
        int end2 = std::min(end, start2 + bucket_size);
        auto tmp =
          get_var_for_remap(mp, jc, builder, node, start2, a.first,
                            start2 + diff_tofrom, end2 - start2, is_key);
        asmjit::X86GpVar *dest;
        if (is_key) {
          if (cur.key_var->getVarType() == asmjit::kVarTypeUIntPtr) {
            OPA_CHECK(0, "not handled");

          } else {
            dest = cur.key_var;
          }
        } else {
          dest = cur.input_vars[target];
        }

        OPA_DISP("GOT INPUT >> ", dest, target, cur.input_vars.size());
        if (dest->getSize() == 1)
          jc.c->add(*dest, tmp->r8());
        else if (dest->getSize() == 2)
          jc.c->add(*dest, tmp->r16());
        else
          jc.c->add(*dest, *tmp);
      }
    }
    puts("");
    x->setup_jit(cur);
  }
  jc.c->mov(*builder.output_var, *mp[output].output_var);
}

std::unique_ptr<asmjit::X86GpVar>
CipherGraph::get_var_for_remap(std::map<const CipherNode *, JitBuilder> &mp,
                               JitContext &jc, const JitBuilder &builder,
                               CipherNode *to, int to_pos, CipherNode *from,
                               int from_pos, int len, bool is_key) const {

  int input_block;
  if (is_key) {
    input_block = to->block->call().key_base_size;
  } else {
    input_block = to->block->call().input_base_size;
  }
  OPA_DISP("get_var_for_remap ", to->block->desc(), to_pos, from->block->desc(),
           from_pos, len);

  int output_block;

  asmjit::X86GpVar *base = nullptr;
  if (from == input) {
    output_block = call().input_base_size;
    OPA_CHECK0(from_pos / call().input_base_size ==
               (from_pos + len - 1) / call().input_base_size);
    base = builder.input_vars[from_pos / output_block];

  } else if (from == key) {
    // TODO: fuck this condition, fuck asmjit
    if (call().key_base_size < key_size()) {
      OPA_CHECK0(from_pos / call().key_base_size ==
                 (from_pos + len - 1) / call().key_base_size);
      jc.vars.emplace_back(new asmjit::X86GpVar(
        *jc.c, context()->jit().sizeToType(call().key_base_size)));
      base = jc.vars.back().get();
      int addr = from_pos / call().key_base_size;
      addr <<= call().key_base_size / 8 - 1;
      asmjit::X86GpVar tmp2(*jc.c, asmjit::kVarTypeUIntPtr);
      jc.c->mov(tmp2, asmjit::Imm(addr));
      jc.c->mov(*base, asmjit::x86::ptr(*builder.key_var, tmp2));
    } else {
      OPA_DISP("HAS KEY INPUT 1", 0);
      base = builder.key_var;
    }

    output_block = call().key_base_size;
  } else {
    base = mp[from].output_var;
    output_block = from->block->output_size();
  }
  std::unique_ptr<asmjit::X86GpVar> res(new asmjit::X86GpVar(
    *jc.c, context()->jit().sizeToType(std::max(output_block, input_block))));

  if (res->getSize() > base->getSize())
    jc.c->movzx(*res, *base);
  else
    jc.c->mov(*res, *base);

  if (len != output_block) { // need masking
    u64 mask = (1ull << (from_pos + len)) - (1 << from_pos);
    jc.c->and_(*res, asmjit::Imm(mask));
    OPA_DISP("MASK >>", mask, output_block, from_pos, to_pos, len);
  }

  int shift = to_pos % input_block - from_pos % output_block;
  OPA_DISP("shift ", from->block->desc(), to->block->desc(), shift);
  if (shift > 0) {
    jc.c->shl(*res, shift);
  } else if (shift != 0) {
    jc.c->shr(*res, -shift);
  }
  return res;
}
void CipherGraph::smart_plug(const std::vector<PlugDesc> &ins, PlugDesc out) {
  OPA_CHECK0(ins.size() > 0);
  int tot = 0;
  for (auto &x : ins)
    tot += x.osz();
  OPA_CHECK_EQ0(tot, out.isz());

  int pos = 0;
  for (auto &x : ins) {
    add_range_edges(x.node, out.node, x.off, out.off + pos, x.osz());
    pos += x.osz();
  }
}

void CipherGraph::do_rel_rec(GraphRelStorer &res,
                             const std::vector<GraphRel> &candidates,
                             const RangeCoverage &targets) const {

  OPA_DISP("Adding do_rec_rel, ", candidates.size());
  // REP (i, candidates.size())
  //  OPA_DISP0(candidates[i].in, candidates[i].out, candidates[i].cost);
  REP (i, candidates.size())
    REP (j, i)
      REP (k, j + 1) {
        GraphRel nw = candidates[i];
        nw.merge(candidates[j]);
        if (k != j)
          nw.merge(candidates[k]);
        if (nw.mid.has_one(targets))
          continue;
        res.add(nw);
      }
}

std::pair<RangeCoverage, RangeCoverage>
CipherGraph::rel_to_io(const GraphRel &rel) const {
  std::pair<RangeCoverage, RangeCoverage> res;
  for (auto x : rel.out.all()) {
    OPA_CHECK0(data[x].out.size() > 0);

    for (auto a : data[x].out) {
      res.ND.add(data[a].pos);
      break;
      // TODO: can generate one relation for each entry
    }
  }

  for (auto x : rel.key.all())
    res.ST.add(data[x].pos + raw_input_size());
  for (auto x : rel.in.all())
    res.ST.add(data[x].pos);

  return res;
}
OPA_NAMESPACE_END(opa, crypto, la)
