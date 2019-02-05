#pragma once

#include <opa/predef.h>
#include <opa/utils/serialize_base.pb.h>
#include <opa/utils/string.h>
#include <opa_common.h>

OPA_NM_UTILS

typedef s32 ClassId;

template <class T> class Singleton {
public:
  T *get() {
    if (obj == nullptr) {
      obj = new T;
    }
    return obj;
  }

private:
  static T *obj;
};
template <class T> T *Singleton<T>::obj = nullptr;

template <class T> class Singleton2 {
public:
  static T *get() {
    if (obj == nullptr) {
      obj = new T;
    }
    return obj;
  }

private:
  static T *obj;
};
#define OPA_SINGLETON_GETTER(T)                                                \
  friend class opa::utils::Singleton2<T>;                                      \
  static T *GetSingleton() { return opa::utils::Singleton2<T>::get(); }

template <class T> T *Singleton2<T>::obj = nullptr;

class ProtobufParams {
public:
  virtual ~ProtobufParams() {}
  virtual void load(const ::opa::utils::AnyMsg &any) = 0;
  virtual void store(::opa::utils::AnyMsg &any) const = 0;
  virtual void after_load() {}
  virtual void before_store() const {}
  void load_entry(const ::opa::utils::AnyMsg &any) {
    this->load(any);
    this->after_load();
  }
  void store_entry(::opa::utils::AnyMsg &any) const {
    this->before_store();
    this->store(any);
  }
};

class BaseStorable : public ProtobufParams {
public:
  ClassId get_class_id() const { return m_id; }
  void set_class_id(ClassId id) { m_id = id; }

  virtual void load(const ::opa::utils::AnyMsg &any) override {}
  virtual void store(::opa::utils::AnyMsg &any) const override {}

private:
  ClassId m_id = -1;
};

template <class T> class Storable : public BaseStorable {
public:
  Storable(const T &a) { cl = a; }
  T &get() { return cl; };
  T &operator*() { return cl; }

private:
  T cl;
  ClassId id;
};

template <class In, class Out> class MapperFunc : public BaseStorable {
public:
  virtual ~MapperFunc() {}
  virtual void operator()(const In &in, Out &out) const = 0;
  std::shared_ptr<MapperFunc<In, Out> > non_owned_sptr() {
    return std::shared_ptr<MapperFunc<In, Out> >(this,
                                                 [](MapperFunc<In, Out> *) {});
  }
};
OPA_DECL_SPTR(MapperFunc<opa::stolen::StringRef COMMA std::string>,
              StrStrMapperFuncSptr);

class ClassStore {
public:
  template <class T> using ClassInstanciator = std::function<T *()>;
  typedef ClassInstanciator<BaseStorable> DefaultInstanciator;

  template <class T>
  ClassId reg(const ClassInstanciator<T> &func, const std::string &key) {
    return reg_base(DefaultInstanciator(func), key);
  }

  ClassId reg_base(const DefaultInstanciator &func, const std::string &key);
  template <class T> ClassId reg2(const ClassInstanciator<T> &func) {
    return reg(func, typeid(T).name());
  }

  bool check_id(ClassId id) const;

  template <class T> T *get_by_key(const std::string &key) const {
    OPA_CHECK(m_bykey.count(key), "key >> ", key.c_str());
    // OPA_DISP("GET KEY >> ", key);
    return get<T>(m_bykey.find(key)->ND);
  }

  template <class T> T *get(ClassId id) const {
    // OPA_DISP("GET by ID >> ", id, uintptr_t(this), m_mp.size());
    OPA_CHECK0(m_mp.count(id));
    T *res = (T *)m_mp.find(id)->ND();
    res->set_class_id(id);
    return res;
  }

  ClassId get_id(const std::string &key) const {
    OPA_CHECK(m_bykey.count(key), key);
    return m_bykey.find(key)->ND;
  }

  template <class T> T *get2() const {
    return this->get_by_key<T>(typeid(T).name());
  }

  bool has_key(const std::string &key) const { return m_bykey.count(key); }
  template <class T> bool has() const { return has_key(typeid(T).name()); }

  static Singleton<ClassStore> instance;

#define OPACS_TYPE(T) ::opa::utils::Storable<T>
#define OPACS_BY_CL_KEY(cl) cl
#define OPACS_GET(key, T)                                                      \
  ::opa::utils::ClassStore::instance.get()->get_by_key<T>(OPA_XSTR(key))

#define OPACS_GETTER_BY_CL(cl)                                                 \
  static SPTR(cl) getCl() {                                                    \
    return SPTR(cl)(OPACS_GET(OPACS_BY_CL_KEY(cl), cl));                       \
  }
#define OPACS_GETTER_BY_KEY(cl, key)                                           \
  static SPTR(cl) getCl() { return SPTR(cl)(OPACS_GET(key, cl)); }

#define OPACS_HAS2(T)                                                          \
  ::opa::utils::ClassStore::instance.get()->get_by_key<T>(OPA_XSTR(key))

#define _OPA_CLASS_STORE_REGISTER_BY_KEY(key, code)                            \
  _OPA_CLASS_STORE_REGISTER_BY_KEY2(key, code)
#define _OPA_CLASS_STORE_REGISTER_BY_KEY2(key, code)                           \
  void opa_classstore_init_##key() {                                           \
    ::opa::utils::ClassStore::instance.get()->reg(                             \
      std::function< ::opa::utils::BaseStorable *()>(code), #key);             \
  }                                                                            \
  OPA_REGISTER_INIT(opa_classstore_init_##key, opa_classstore_init_##key)

#define _OPA_CLASS_STORE_REGISTER_BY_T(T, code)                                \
  namespace {                                                                  \
  void opa_classstore_init_##T() {                                             \
    ::opa::utils::ClassStore::instance.get()->reg2<T>(                         \
      std::function<T *()>(code));                                             \
  }                                                                            \
  OPA_REGISTER_INIT(opa_classstore_init_##T, opa_classstore_init_##T)          \
  }

#define OPA_CLASS_STORE_REGISTER_BY_KEY(cl, key)                               \
  _OPA_CLASS_STORE_REGISTER_BY_KEY(key, []() { return new cl; })

#define OPA_CLASS_STORE_REGISTER_BY_CL(cl)                                     \
  _OPA_CLASS_STORE_REGISTER_BY_KEY(OPACS_BY_CL_KEY(cl), []() { return new cl; })

#define OPA_CLASS_STORE_REGISTER_BY_KEY_MAKE(cl, key, typ)                     \
  _OPA_CLASS_STORE_REGISTER_BY_KEY(                                            \
    key, []() { return new ::opa::utils::Storable<typ>(cl); })

private:
  static ClassStore *g_instance;
  std::map<ClassId, DefaultInstanciator> m_mp;
  std::map<std::string, ClassId> m_bykey;
};

#define OPA_LOADER_ONE(field, id) ::opa::utils::tgen_load(field, any.any(id));
#define OPA_STORER_ONE(field, id)                                              \
  ::opa::utils::tgen_store(field, *any.add_any());

#define OPA_TGEN_IMPL(...)                                                     \
  virtual void load(const ::opa::utils::AnyMsg &any) override {                \
    OPA_TGEN_LOADER(__VA_ARGS__);                                              \
  }                                                                            \
  virtual void store(::opa::utils::AnyMsg &any) const override {               \
    OPA_TGEN_STORER(__VA_ARGS__);                                              \
  }

#define OPA_TGEN_IMPL_INIT2(...)                                               \
  virtual void load(const ::opa::utils::AnyMsg &any) override {                \
    this->init();                                                              \
    OPA_TGEN_LOADER(__VA_ARGS__);                                              \
  }                                                                            \
  virtual void store(::opa::utils::AnyMsg &any) const override {               \
    OPA_TGEN_STORER(__VA_ARGS__);                                              \
  }

#define OPA_TGEN_IMPL_INIT(...)                                                \
  OPA_TGEN_BASE_LOADER(typ, prototyp, OPA_TGEN_LOADER(__VA_ARGS__);            \
                       init(__VA_ARGS__));                                     \
  OPA_TGEN_BASE_STORER(typ, prototyp, OPA_TGEN_STORER(__VA_ARGS__));

#define OPA_TGEN_CODE(loadcode, storecode)                                     \
  virtual void load(const ::opa::utils::AnyMsg &any) override { loadcode; }    \
  virtual void store(::opa::utils::AnyMsg &any) const override { storecode; }

#define OPA_TGEN_LOADER(...)                                                   \
  OPA_DISPATCH(_OPA_TGEN_LOADER, __VA_ARGS__)(__VA_ARGS__)

#define _OPA_TGEN_LOADER0()                                                    \
  {}
#define _OPA_TGEN_LOADER1(a1) OPA_LOADER_ONE(a1, 0)
#define _OPA_TGEN_LOADER2(a1, a2)                                              \
  _OPA_TGEN_LOADER1(a1);                                                       \
  OPA_LOADER_ONE(a2, 1);

#define _OPA_TGEN_LOADER3(a1, a2, a3)                                          \
  _OPA_TGEN_LOADER2(a1, a2);                                                   \
  OPA_LOADER_ONE(a3, 2);

#define _OPA_TGEN_LOADER4(a1, a2, a3, a4)                                      \
  _OPA_TGEN_LOADER3(a1, a2, a3);                                               \
  OPA_LOADER_ONE(a4, 3);

#define _OPA_TGEN_LOADER5(a1, a2, a3, a4, a5)                                  \
  _OPA_TGEN_LOADER4(a1, a2, a3, a4);                                           \
  OPA_LOADER_ONE(a5, 4);

#define _OPA_TGEN_LOADER6(a1, a2, a3, a4, a5, a6)                              \
  _OPA_TGEN_LOADER5(a1, a2, a3, a4, a5);                                       \
  OPA_LOADER_ONE(a6, 5);

#define _OPA_TGEN_LOADER7(a1, a2, a3, a4, a5, a6, a7)                          \
  _OPA_TGEN_LOADER6(a1, a2, a3, a4, a5, a6);                                   \
  OPA_LOADER_ONE(a7, 6);

#define _OPA_TGEN_LOADER8(a1, a2, a3, a4, a5, a6, a7, a8)                      \
  _OPA_TGEN_LOADER7(a1, a2, a3, a4, a5, a6, a7);                               \
  OPA_LOADER_ONE(a8, 7);

#define _OPA_TGEN_LOADER9(a1, a2, a3, a4, a5, a6, a7, a8, a9)                  \
  _OPA_TGEN_LOADER8(a1, a2, a3, a4, a5, a6, a7, a8);                           \
  OPA_LOADER_ONE(a9, 8);

#define _OPA_TGEN_LOADER10(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)            \
  _OPA_TGEN_LOADER9(a1, a2, a3, a4, a5, a6, a7, a8, a9);                       \
  OPA_LOADER_ONE(a10, 9);

#define OPA_TGEN_STORER(...)                                                   \
  OPA_DISPATCH(_OPA_TGEN_STORER, __VA_ARGS__)(__VA_ARGS__)

#define _OPA_TGEN_STORER1(a1) OPA_STORER_ONE(a1, 0)
#define _OPA_TGEN_STORER2(a1, a2)                                              \
  _OPA_TGEN_STORER1(a1);                                                       \
  OPA_STORER_ONE(a2, 1);

#define _OPA_TGEN_STORER3(a1, a2, a3)                                          \
  _OPA_TGEN_STORER2(a1, a2);                                                   \
  OPA_STORER_ONE(a3, 2);

#define _OPA_TGEN_STORER4(a1, a2, a3, a4)                                      \
  _OPA_TGEN_STORER3(a1, a2, a3);                                               \
  OPA_STORER_ONE(a4, 3);

#define _OPA_TGEN_STORER5(a1, a2, a3, a4, a5)                                  \
  _OPA_TGEN_STORER4(a1, a2, a3, a4);                                           \
  OPA_STORER_ONE(a5, 4);

#define _OPA_TGEN_STORER6(a1, a2, a3, a4, a5, a6)                              \
  _OPA_TGEN_STORER5(a1, a2, a3, a4, a5);                                       \
  OPA_STORER_ONE(a6, 5);

#define _OPA_TGEN_STORER7(a1, a2, a3, a4, a5, a6, a7)                          \
  _OPA_TGEN_STORER6(a1, a2, a3, a4, a5, a6);                                   \
  OPA_STORER_ONE(a7, 6);

#define _OPA_TGEN_STORER8(a1, a2, a3, a4, a5, a6, a7, a8)                      \
  _OPA_TGEN_STORER7(a1, a2, a3, a4, a5, a6, a7);                               \
  OPA_STORER_ONE(a8, 7);

#define _OPA_TGEN_STORER9(a1, a2, a3, a4, a5, a6, a7, a8, a9)                  \
  _OPA_TGEN_STORER8(a1, a2, a3, a4, a5, a6, a7, a8);                           \
  OPA_STORER_ONE(a9, 8);

#define _OPA_TGEN_STORER10(a1, a2, a3, a4, a5, a6, a7, a8, a9, a10)            \
  _OPA_TGEN_STORER9(a1, a2, a3, a4, a5, a6, a7, a8, a9);                       \
  OPA_STORER_ONE(a10, 9);

#define OPA_TGEN_BASE_LOADER(typ, prototyp, code)                              \
  static void tgen_load(typ &x, const google::protobuf::Any &any) {            \
    prototyp y;                                                                \
    any.UnpackTo(&y);                                                          \
    code;                                                                      \
  }
#define OPA_TGEN_BASE_STORER(typ, prototyp, code)                              \
  static void tgen_store(const typ &x, google::protobuf::Any &any) {           \
    prototyp y;                                                                \
    code;                                                                      \
    any.PackFrom(y);                                                           \
  }

#define OPA_TGEN_BASE_LOADER_SPE(typ, spe, prototyp, code)                     \
  template <> void tgen_load<spe>(typ & x, const google::protobuf::Any &any) { \
    prototyp y;                                                                \
    any.UnpackTo(&y);                                                          \
    code;                                                                      \
  }
#define OPA_TGEN_BASE_STORER_SPE(typ, spe, prototyp, code)                     \
  template <> void tgen_store<spe>(const typ &x, google::protobuf::Any &any) { \
    prototyp y;                                                                \
    code;                                                                      \
    any.PackFrom(y);                                                           \
  }

#define OPA_TGEN_BASE(typ, prototyp)                                           \
  OPA_TGEN_BASE_LOADER(typ, prototyp, x = y.x());                              \
  OPA_TGEN_BASE_STORER(typ, prototyp, y.set_x(x));

#define OPA_TGEN_BASE2(typ, prototyp, codeload, codestore)                     \
  OPA_TGEN_BASE_LOADER(typ, prototyp, codeload);                               \
  OPA_TGEN_BASE_STORER(typ, prototyp, codestore);

#define OPA_TGEN_BASE2_SPE(typ, spe, prototyp, codeload, codestore)            \
  OPA_TGEN_BASE_LOADER_SPE(typ, spe, prototyp, codeload);                      \
  OPA_TGEN_BASE_STORER_SPE(typ, spe, prototyp, codestore);

#define ARGS_TO_TEMPLATE(...)                                                  \
  OPA_DISPATCH(_ARGS_TO_TEMPLATE, __VA_ARGS__)(__VA_ARGS__)
#define _ARGS_TO_TEMPLATE1(a) class a
#define _ARGS_TO_TEMPLATE2(a, b) class a, class b
#define _ARGS_TO_TEMPLATE3(a, b, c) class a, class b, class c

#define ARGS_TO_TUPLE(...)                                                     \
  OPA_DISPATCH(_ARGS_TO_TUPLE, __VA_ARGS__)(__VA_ARGS__)
#define _ARGS_TO_TUPLE1(a) a
#define _ARGS_TO_TUPLE2(a, b) a, b
#define _ARGS_TO_TUPLE3(a, b, c) a, b, c

#define OPA_TGEN_BASE_TEMPLATE_RAW(typ, prototyp, codeload, codestore, ...)    \
  template <ARGS_TO_TEMPLATE(__VA_ARGS__)>                                     \
  OPA_TGEN_BASE_LOADER(typ<ARGS_TO_TUPLE(__VA_ARGS__)>, prototyp, codeload);   \
  template <ARGS_TO_TEMPLATE(__VA_ARGS__)>                                     \
  OPA_TGEN_BASE_STORER(typ<ARGS_TO_TUPLE(__VA_ARGS__)>, prototyp, codestore);

#define OPA_TGEN_BASE_TEMPLATE(typ, prototyp, codeload, codestore)             \
  template <class T> OPA_TGEN_BASE_LOADER(typ, prototyp, codeload);            \
  template <class T> OPA_TGEN_BASE_STORER(typ, prototyp, codestore);

OPA_TGEN_BASE(std::string, ::opa::utils::BaseString);
OPA_TGEN_BASE(float, ::opa::utils::BaseFloat);
OPA_TGEN_BASE(double, ::opa::utils::BaseDouble);
OPA_TGEN_BASE(char, ::opa::utils::BaseChar);
OPA_TGEN_BASE(u32, ::opa::utils::BaseInt32);
OPA_TGEN_BASE(s32, ::opa::utils::BaseInt32);
OPA_TGEN_BASE(u64, ::opa::utils::BaseInt64);
OPA_TGEN_BASE(ull, ::opa::utils::BaseInt64);
OPA_TGEN_BASE(ll, ::opa::utils::BaseInt64);
OPA_TGEN_BASE(s64, ::opa::utils::BaseInt64);
OPA_TGEN_BASE(u16, ::opa::utils::BaseInt32);
OPA_TGEN_BASE(s16, ::opa::utils::BaseInt32);

template <class T>
static void tgen_store(const std::vector<T> &x, google::protobuf::Any &any);

template <class T>
static void tgen_load(std::vector<T> &x, const google::protobuf::Any &any);

OPA_TGEN_BASE(bool, ::opa::utils::BaseBool);
OPA_TGEN_BASE2(ProtobufParams, ::opa::utils::AnyMsg, { x.load_entry(y); },
               { x.store_entry(y); });
OPA_TGEN_BASE_TEMPLATE(::opa::utils::Storable<std::function<T> >,
                       ::opa::utils::AnyMsg, {}, {})

OPA_TGEN_BASE_TEMPLATE(
  std::shared_ptr<T>, ::opa::utils::AnyMsg,
  {
    int id;
    ::opa::utils::tgen_load(id, y.any(0));
    OPA_CHECK(::opa::utils::ClassStore::instance.get()->check_id(id), id);
    x.reset(::opa::utils::ClassStore::instance.get()->get<T>(id));
    ::opa::utils::tgen_load(*x, y.any(1));
  },
  {
    y.add_any();
    y.add_any();
    OPA_CHECK(
      ::opa::utils::ClassStore::instance.get()->check_id(x->get_class_id()),
      x->get_class_id());
    ::opa::utils::tgen_store(x->get_class_id(), *y.mutable_any(0));
    ::opa::utils::tgen_store(*x, *y.mutable_any(1));
  });

OPA_TGEN_BASE_TEMPLATE(
  std::shared_ptr<const T>, ::opa::utils::AnyMsg,
  {
    int id;
    ::opa::utils::tgen_load(id, y.any(0));
    OPA_CHECK(::opa::utils::ClassStore::instance.get()->check_id(id), id);
    T *tmp = ::opa::utils::ClassStore::instance.get()->get<T>(id);
    ::opa::utils::tgen_load(*tmp, y.any(1));
    x.reset(tmp);
  },
  {
    y.add_any();
    y.add_any();
    OPA_CHECK(
      ::opa::utils::ClassStore::instance.get()->check_id(x->get_class_id()),
      x->get_class_id());
    ::opa::utils::tgen_store(x->get_class_id(), *y.mutable_any(0));
    ::opa::utils::tgen_store(*x, *y.mutable_any(1));
  });

OPA_TGEN_BASE_TEMPLATE_RAW(std::pair, ::opa::utils::AnyMsg,
                           {
                             ::opa::utils::tgen_load(x.ST, y.any(0));
                             ::opa::utils::tgen_load(x.ND, y.any(1));
                           },
                           {
                             y.add_any();
                             y.add_any();
                             ::opa::utils::tgen_store(x.ST, *y.mutable_any(0));
                             ::opa::utils::tgen_store(x.ND, *y.mutable_any(1));
                           },
                           A1, A2);

OPA_TGEN_BASE_TEMPLATE_RAW(
  std::tuple, ::opa::utils::AnyMsg,
  {
    ::opa::utils::tgen_load(std::get<0>(x), y.any(0));
    ::opa::utils::tgen_load(std::get<1>(x), y.any(1));
    ::opa::utils::tgen_load(std::get<2>(x), y.any(2));
  },
  {
    y.add_any();
    y.add_any();
    y.add_any();
    ::opa::utils::tgen_store(std::get<0>(x), *y.mutable_any(0));
    ::opa::utils::tgen_store(std::get<1>(x), *y.mutable_any(1));
    ::opa::utils::tgen_store(std::get<2>(x), *y.mutable_any(2));
  },
  A1, A2, A3);

OPA_TGEN_BASE_TEMPLATE_RAW(std::tuple, ::opa::utils::AnyMsg,
                           {
                             ::opa::utils::tgen_load(std::get<0>(x), y.any(0));
                             ::opa::utils::tgen_load(std::get<1>(x), y.any(1));
                           },
                           {
                             y.add_any();
                             y.add_any();
                             ::opa::utils::tgen_store(std::get<0>(x),
                                                      *y.mutable_any(0));
                             ::opa::utils::tgen_store(std::get<1>(x),
                                                      *y.mutable_any(1));
                           },
                           A1, A2);

OPA_TGEN_BASE_TEMPLATE(std::vector<T>, ::opa::utils::AnyMsg,
                       {
                         x.resize(y.any_size());
                         REP (i, x.size())
                           tgen_load(x[i], y.any(i));
                       },
                       {
                         REP (i, x.size())
                           tgen_store(x[i], *y.add_any());
                       });

OPA_TGEN_BASE_TEMPLATE(std::set<T>, ::opa::utils::AnyMsg,
                       {
                         REP (i, y.any_size()){
                         T tmp;
                           tgen_load(tmp, y.any(i));
                           x.insert(tmp);
}
}
, {
  for (auto &k : x) tgen_store(k, *y.add_any());
});

template <class T>
static void anymsg_load(T &a, const ::opa::utils::AnyMsg &any);
template <class T>
static void anymsg_store(const T &a, ::opa::utils::AnyMsg &any);

template <>
void anymsg_load(opa::utils::ProtobufParams &a,
                 const ::opa::utils::AnyMsg &any) {
  a.load_entry(any);
}

template <>
void anymsg_store(const opa::utils::ProtobufParams &a,
                  ::opa::utils::AnyMsg &any) {
  a.store_entry(any);
}

template <class T> void anymsg_load(T &a, const ::opa::utils::AnyMsg &any) {
  OPA_LOADER_ONE(a, 0);
}

template <class T> void anymsg_store(const T &a, ::opa::utils::AnyMsg &any) {
  OPA_STORER_ONE(a, 0);
}

struct DoubleRes : public ProtobufParams {

  DoubleRes(double v) { this->v = v; }
  DoubleRes() {}
  double v;
  OPA_TGEN_IMPL(v);
};

OPA_NM_UTILS_END
