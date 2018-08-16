#pragma once
#include <opa_common.h>

OPA_NAMESPACE(opa, utils)

//class RegisterObj {
//public:
//  virtual ~RegisterObj() {}
//};
//
//class Register {
//public:
//  void reg(const std::string &key, UPTR(RegisterObj) obj) {
//    OPA_CHECK(m_objs.count(key) == 0, key);
//    m_objs[key] = std::move(obj);
//  }
//
//  template <class T> T *get(const std::string &key) const {
//    OPA_CHECK(m_objs.count(key), key);
//    return (T *)m_objs.find(key)->ND.get();
//  }
//
//  static Singleton<Register> singleton;
//
//private:
//  std::map<std::string, UPTR(RegisterObj)> m_objs;
//};

OPA_NAMESPACE_END(opa, utils)
