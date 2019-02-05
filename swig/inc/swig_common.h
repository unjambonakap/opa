
struct swig_tsf_data {
  std::vector<std::pair<PyObject *, bool> > objs;
  std::string str;
  opa::stolen::StringRef ref;
  glib::StringPiece piece;
  bool ok = true;
  int n;
  u8 *str_buf;

  glib::StringPiece get_piece(PyObject *init) { return *get_piece_ptr(init); }
  glib::StringPiece *get_piece_ptr(PyObject *init) {
    get_u8(init);
    piece = glib::StringPiece((const char *)str_buf, n);
    return &piece;
  }
  opa::stolen::StringRef get_ref(PyObject *init) { return *get_ref_ptr(init); }

  opa::stolen::StringRef *get_ref_ptr(PyObject *init) {
    get_u8(init);
    ref = opa::stolen::StringRef((const char *)str_buf, n);
    return &ref;
  }

  void add_obj(PyObject *obj, bool created) { objs.emplace_back(obj, created); }

  bool check(PyObject *obj) const {
    return obj != nullptr && (PyBytes_Check(obj) || PyUnicode_Check(obj));
  }

  u8 *tsf_str(PyObject *init, bool owned = false) {
    add_obj(init, owned);
    if (PyUnicode_Check(init)) add_obj(PyUnicode_AsUTF8String(init), true);
    if (PyByteArray_Check(init)) add_obj(PyBytes_FromObject(init), true);
    if (!PyBytes_Check(objs.back().first)) {
      PyErr_SetString(PyExc_ValueError, "Expected a bytes shit");
      ok = false;
      return nullptr;
    }

    n = PyBytes_Size(objs.back().first);
    return (u8 *)PyBytes_AsString(objs.back().first);
  }

  std::string *get_str(PyObject *init, bool owned = false) {
    if (!get_u8(init, owned)) return nullptr;
    str = std::string((const char *)str_buf, n);
    return &str;
  }

  u8 *get_u8(PyObject *init, bool owned = false) {
    str_buf = tsf_str(init, owned);
    return str_buf;
  }

  ~swig_tsf_data() {
    for (auto &e : objs) {
      if (e.second) Py_XDECREF(e.first);
    }
  }

  std::string get_obj_str(PyObject *obj) {
    PyObject *str_obj = PyObject_Str(obj);
    if (str_obj == nullptr) return "NULLOBJ";
    return *get_str(obj, true);
  }
};

namespace swig_helper {

PyObject *create_list() { return PyList_New(0); }
template <typename T> PyObject *convert(const T &x) { OPA_CHECK0(false); }
PyObject *convert(PyObject *o) { return o; }
PyObject *convert(float v) { return PyFloat_FromDouble(v); }
PyObject *convert(bool v) {
  if (v) Py_RETURN_TRUE;
  Py_RETURN_FALSE;
}
PyObject *convert(double v) { return PyFloat_FromDouble(v); }
PyObject *convert(int v) { return PyInt_FromLong(v); }
PyObject *convert(const std::string &v) {
  return PyBytes_FromStringAndSize(v.data(), v.size());
}

template <typename T> void add_to_list(PyObject *list, const T &obj) {
  PyList_Append(list, convert(obj));
}

template <typename T> PyObject *convert(const std::vector<T> &v) {
  PyObject *res = create_list();
  for (auto &e : v) add_to_list(res, e);
  return res;
}
template <typename T, std::size_t N>
PyObject *convert(const std::array<T, N> &v) {
  std::vector<T> vec;
  for (int i = 0; i < N; ++i) vec.push_back(v[i]);
  return convert(vec);
}

template <typename T> PyObject *convert(const T *v) { return convert(*v); }
template <typename U, typename V> PyObject *convert(const std::pair<U, V> &x) {
  PyObject *res = create_list();
  add_to_list(res, x.first);
  add_to_list(res, x.second);
  return res;
}
template <typename T> PyObject *convert(const SwigValueWrapper<T> &v) {
  return convert((T) * (&v));
}

#define DEF_LOAD_FUNC(typ, func_check, func_conv)                              \
  bool load(PyObject *o, typ &v) {                                             \
    if (func_check(o)) {                                                       \
      v = (typ)func_conv(o);                                                   \
      return true;                                                             \
    }                                                                          \
    return false;                                                              \
  }

DEF_LOAD_FUNC(float, PyNumber_Check, PyFloat_AsDouble);
DEF_LOAD_FUNC(double, PyNumber_Check, PyFloat_AsDouble);
DEF_LOAD_FUNC(int, PyInt_Check, PyInt_AsLong);

template <typename T>
bool load_vec(PyObject *obj, std::vector<T> &res, int expected_length = -1) {
  if (!PySequence_Check(obj)) {
    return false;
  }

  int length = PySequence_Length(obj);
  if (expected_length != -1 && expected_length != length) {
    return false;
  }

  res.resize(length);
  REP (i, length) {
    PyObject *no = PySequence_GetItem(obj, i);
    res.emplace_back();
    if (!load(no, res[i])) return false;
  }
  return true;
}

bool load(PyObject *o, bool &res) {
  if (!PyBool_Check(o)) return false;
  return o == Py_True;
}

bool load(PyObject *o, std::string &res) {
  swig_tsf_data data;
  res = *data.get_str(o);
  return data.ok;
}

template <typename T> bool load(PyObject *o, std::vector<T> &res) {
  return load_vec(o, res);
}
template <typename T, int N> bool load(PyObject *o, std::array<T, N> &res) {
  std::vector<T> tb;
  if (!load_vec(o, tb, N)) return false;
  REP (i, N)
    res[i] = tb[i];
  return true;
}

} // namespace swig_helper

std::string obj_str(PyObject *obj) { return swig_tsf_data().get_obj_str(obj); }
