#include <opa_common.h>

// trnasition to swig::from?

template <int N, typename BaseType>
std::vector<BaseType> glm_to_vec(const glm::vec<N, BaseType> &e) {
  std::vector<BaseType> res(N);
  REP (i, N)
    res[i] = e[i];
  return res;
}

template <int N, typename BaseType>
glm::vec<N, BaseType> vec_to_glm(const std::vector<BaseType> &e) {
  OPA_CHECK0(e.size() == N);
  glm::vec<N, BaseType> res;
  REP (i, N)
    res[i] = e[i];
  return res;
}

template <int N, typename BaseType>
std::vector<BaseType> glm_to_mat(const glm::mat<N, N, BaseType> &e) {
  std::vector<BaseType> res(N * N);
  REP (i, N)
    REP (j, N)
      res[(i * N + j)] = e[i][j];
  return res;
}

template <int N, typename BaseType>
glm::mat<N, N, BaseType> mat_to_glm(const std::vector<BaseType> &e) {
  OPA_CHECK0(e.size() == N * N);
  glm::mat<N, N, BaseType> res;
  REP (i, N)
    REP (j, N)
      res[i][j] = e[(i * N + j)];
  return res;
}

struct SwigVectorHelper {
  PyObject *obj;
  SwigVectorHelper() { create_list(); }
  SwigVectorHelper &create_list() {
    obj = PyList_New(0);
    return *this;
  }

  SwigVectorHelper &add(PyObject *o) {
    PyList_Append(obj, o);
    return *this;
  }

  SwigVectorHelper &add(float v) {
    PyList_Append(obj, PyFloat_FromDouble(v));
    return *this;
  }

  SwigVectorHelper &add(double v) {
    PyList_Append(obj, PyFloat_FromDouble(v));
    return *this;
  }

  SwigVectorHelper &add(int v) {
    PyList_Append(obj, PyInt_FromLong(v));
    return *this;
  }

  template <class T> SwigVectorHelper &add_vec(const std::vector<T> &v) {
    for (auto &e : v) {
      add(e);
    }
    return *this;
  }

  SwigVectorHelper &add(const glm::quat &v) {
    std::vector<double> vec;
    for (int i = 0; i < 4; ++i) vec.push_back(v[i]);
    return add_vec(vec);
  }

  template <int N, typename BaseType>
  SwigVectorHelper &add(const glm::vec<N, BaseType> &v) {
    return add_vec(glm_to_vec(v));
  }

  template <int N, typename BaseType>
  SwigVectorHelper &add(const glm::mat<N, N, BaseType> &v) {
    return add_vec(glm_to_mat(v));
  }

  template <int N, typename BaseType>
  SwigVectorHelper &add(const glm::vec<N, BaseType> *v) {
    return add(*v);
  }

  template <int N, typename BaseType>
  SwigVectorHelper &add(const glm::mat<N, N, BaseType> *v) {
    return add(*v);
  }
  SwigVectorHelper &add(const glm::quat *v) { return add(*v); }

  template <typename T, std::size_t N>
  SwigVectorHelper &add(const std::array<T, N> &v) {
    std::vector<T> vec;
    for (int i = 0; i < N; ++i) vec.push_back(v[i]);
    return add_vec(vec);
  }

  template <typename T, std::size_t N>
  SwigVectorHelper &add(const std::array<T, N> *v) {
    return add(*v);
  }

  template <class T> SwigVectorHelper &add(const SwigValueWrapper<T> &v) {
    return this->add(*(&v));
  }
};


struct SwigSeqLoader {
  bool ok = false;
  bool check_type;

  SwigSeqLoader(bool check_type = false) : check_type(check_type) {}

  void set_err(const std::string &str) const {
    if (!check_type) PyErr_SetString(PyExc_ValueError, str.c_str());
  }

  template <class T>
  std::vector<T> load(PyObject *obj, int expected_length = -1) {
    std::vector<T> res;
    ok = this->convert_vec(obj, res, expected_length);
    return res;
  }

  template <typename T> T load_obj(PyObject *o) const {
    T res;
    OPA_CHECK0(this->convert(o, res));
    return res;
  }

  template <typename T, std::size_t N>
  bool convert(PyObject *o, std::array<T, N> &dest) const {
    SwigSeqLoader load2;
    auto tmp = load2.load<T>(o, N);
    if (!load2.ok) return false;
    REP (i, N)
      dest[i] = tmp[i];
    return true;
  }

  template <typename BaseType>
  bool convert_vec(PyObject *o, std::vector<BaseType> &dest, int n = -1) {

    if (!PySequence_Check(o)) {
      set_err(OPA_TRACE_STR("Expected a sequence, got ", obj_str(o)));
      return false;
    }

    int length = PySequence_Length(o);
    if (n != -1 && n != length) {
      set_err(OPA_TRACE_STR("Bad length ", n, length));
      return false;
    }

    REP (i, length) {
      PyObject *no = PySequence_GetItem(o, i);
      dest.emplace_back();
      if (!this->convert(no, dest.back())) return false;
    }
    return true;
  }

  template <int N, typename BaseType>
  bool convert(PyObject *o, glm::vec<N, BaseType> &dest) {
    std::vector<BaseType> tmp;
    if (!convert_vec(o, tmp, N)) return false;
    dest = vec_to_glm<N, BaseType>(tmp);
    return true;
  }

  template <int N, typename BaseType>
  bool convert(PyObject *o, glm::mat<N, N, BaseType> &dest) {
    std::vector<BaseType> tmp;
    if (!convert_vec(o, tmp, N)) return false;
    dest = mat_to_glm<N, BaseType>(tmp);
    return true;
  }

  bool convert(PyObject *o, float &dest) const {
    if (PyNumber_Check(o)) {
      dest = (double)PyFloat_AsDouble(o);
      return true;
    }
    set_err(
      OPA_TRACE_STR("Sequence elements must be numbers, got: ", obj_str(o)));
    return false;
  }

  bool convert(PyObject *o, double &dest) const {
    if (PyNumber_Check(o)) {
      dest = (double)PyFloat_AsDouble(o);
      return true;
    }
    set_err(
      OPA_TRACE_STR("Sequence elements must be numbers, got: ", obj_str(o)));
    return false;
  }

  bool convert(PyObject *o, int &dest) const {
    if (PyInt_Check(o)) {
      dest = (double)PyInt_AsLong(o);
      return true;
    }
    set_err(
      OPA_TRACE_STR("Sequence elements must be numbers, got: ", obj_str(o)));
    return false;
  }
  template <class T> bool convert(PyObject *o, std::vector<T> &dest) {
    return this->convert_vec(o, dest);
  }

  template <class T> bool can_convert(PyObject *o) {
    check_type = true;
    T e;
    return this->convert(o, e);
  }
};
