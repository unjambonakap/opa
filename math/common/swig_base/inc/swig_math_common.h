namespace opa{
namespace math{
  namespace common{
    static PyObject *PyFromBignum(const bignum &s) {
      std::string sv = s.str(16);
      return PyLong_FromString(sv.data(), nullptr, 16);
    }

    bool PyToBignum(PyObject *o, bignum *output) {
      swig_tsf_data tmp;
      std::string  *x =tmp.get_repr(o);
      if (!x) goto bad;
      *output = bignum::fromstr(*x, 10);
      if (output->bad) goto bad;
      return true;
bad:;
        PyErr_SetString(PyExc_ValueError, "Bad integer provided");
        return false;
    }
  }
}
}
