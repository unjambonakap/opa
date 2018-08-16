#pragma once

#include <opa_common.h>
#include <yaml-cpp/yaml.h>

OPA_NAMESPACE_DECL2(opa, utils)

class ScopedDebugCtx;
class DebugCtx {
public:
  DebugCtx() { reset(); }

  void reset() {
    node = YAML::Node();
    cur = &node;
  }

  std::string str() const { return OPA_STREAM_STR(node); }
  std::string cur_str() const { return OPA_STREAM_STR(*cur); }
  void set(const char *s, double v) { (*cur)[s] = v; }
  void set(const char *s, int v) { (*cur)[s] = v; }
  void set(const char *s, float v) { (*cur)[s] = v; }

  template <class T> void set(const char *s, T v) {
    (*cur)[s] = OPA_STREAM_STR(v);
  }
  template <class T> void push_back(const char *s, T v) {
    (*cur)[s].push_back(OPA_STREAM_STR(v));
  }
  void push_back(const char *s, double v) { (*cur)[s].push_back(v); }
  void push_back(const char *s, int v) { (*cur)[s].push_back(v); }
  void push_back(const char *s, float v) { (*cur)[s].push_back(v); }

  std::string flush_and_str();
  std::vector<ScopedDebugCtx *> cur_ctx;
  YAML::Node node;
  YAML::Node *cur;
};

class ScopedDebugCtx {
public:
  ScopedDebugCtx(DebugCtx &ctx, const std::string &key = "") : ctx(ctx) {
    this->key = key;
    prev = ctx.cur;
    ctx.cur = &cur;
    depth = ctx.cur_ctx.size();
    ctx.cur_ctx.push_back(this);
  }

  void flush1() {
    // one time usage, can be used on crashes as well
    if (key.empty()) {
      (*prev)["children"].push_back(cur);
    } else {
      (*prev)[key] = cur;
    }
  }

  void flush_all() {
    FORV (i, depth + 1, ctx.cur_ctx.size())
      ctx.cur_ctx[i]->flush1();
  }

  std::string flush_and_str() {
    flush_all();
    return OPA_STREAM_STR(cur);
  }

  ~ScopedDebugCtx() {
    flush1();

    ctx.cur_ctx.pop_back();
    ctx.cur = prev;
  }

  int depth;
  std::string key;
  YAML::Node children;
  YAML::Node cur;
  DebugCtx &ctx;
  YAML::Node *prev;
};



OPA_NAMESPACE_DECL2_END
