namespace opa{
namespace math{
  namespace common{
    static PyObject *PyFromBignum(const bignum &s) {
      std::string sv = s.str(16);
      return PyLong_FromString(sv.data(), nullptr, 16);
    }
  }
}
}
