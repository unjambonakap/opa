
struct swig_tsf_data {
  std::vector<std::pair<PyObject *, bool> > objs;
  std::string str;
  opa::stolen::StringRef ref;
  glib::StringPiece piece;
  bool ok = true;
  int n;
  u8 *str_buf;
  glib::StringPiece get_piece(PyObject *init) {
    get_u8(init);
    piece = glib::StringPiece((const char *)str_buf, n);
    return piece;
  }
  glib::StringPiece *get_piece_ptr(PyObject *init) {
    get_u8(init);
    piece = glib::StringPiece((const char *)str_buf, n);
    return &piece;
  }
  opa::stolen::StringRef *get_ref(PyObject *init) {
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

  std::string get_obj_str(PyObject *obj){
    PyObject *str_obj = PyObject_Str(obj);
    if (str_obj==nullptr) return "NULLOBJ";
    return *get_str(obj, true);
  }
};

struct SwigHelper {
  static PyObject *FromString(const std::string &s) {
    return PyBytes_FromStringAndSize(s.data(), s.size());
  }
};
