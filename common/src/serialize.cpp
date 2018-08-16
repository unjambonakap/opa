#include <opa/utils/serialize.h>

OPA_NAMESPACE(opa, utils)
Singleton<ClassStore> ClassStore::instance;

ClassId ClassStore::reg_base(const DefaultInstanciator &func,
                           const std::string &key) {
  ClassId id = m_mp.size();
  m_mp[id] = func;
  OPA_DISPERR("REGISTER ", key, id, uintptr_t(this), m_mp.size());
  if (key.size() > 0) {
    OPA_CHECK0(!m_bykey.count(key));
    m_bykey[key] = id;
  }

  return id;
}
bool ClassStore::check_id(ClassId id) const { return m_mp.count(id); }

OPA_NAMESPACE_END(opa, utils)
