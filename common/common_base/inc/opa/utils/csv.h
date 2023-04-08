#pragma once

#include <opa/utils/stacktrace.h>
#include <opa/utils/string_base.h>
#include <opa_common_base.h>

OPA_NAMESPACE(opa, utils)

#if OPA_PIN == 0
#define DEFINE_BASE_SPTR_FUNC_1(cl, T)                                         \
  std::shared_ptr<cl<T> > get_sptr() { return make_sptr(this); }               \
  std::shared_ptr<cl<T> > get_dummy_sptr() { return make_dummy_sptr(this); }
#define DEFINE_BASE_SPTR_FUNC_2(cl, T, U)                                      \
  std::shared_ptr<cl<T, U> > get_sptr() { return make_sptr(this); }            \
  std::shared_ptr<cl<T, U> > get_dummy_sptr() { return make_dummy_sptr(this); }

template <typename T, typename Stream> class CsvRecordWriter {
public:
  virtual ~CsvRecordWriter() {}
  virtual std::vector<std::string> get_headers() const = 0;
  virtual void dump(Stream &s, const T &a) const = 0;
  DEFINE_BASE_SPTR_FUNC_2(CsvRecordWriter, T, Stream);
};
OPA_DECL_SPTR_TMPL2(CsvRecordWriter, CsvRecordWriterSptr);

template <typename T, typename Stream>
class AggregateCsvRecordWriter
    : public CsvRecordWriter<std::vector<T>, Stream> {
public:
  void init(const std::vector<CsvRecordWriterSptr<T, Stream> > &tb) {
    m_tb = tb;
  }

  void init(CsvRecordWriterSptr<T, Stream> writer, int count) {
    REP (i, count)
      m_tb.pb(writer);
  }

  virtual std::vector<std::string> get_headers() const override {
    std::vector<std::string> res;
    REP (i, m_tb.size()) {
      auto tmp = m_tb[i]->get_headers();
      for (auto &s : tmp) s += utils::stdsprintf("_%02d", i);
      res.insert(res.end(), ALL(tmp));
    }
    return res;
  }

  virtual void dump(Stream &s, const std::vector<T> &a) const override {
    OPA_CHECK0(a.size() == m_tb.size());
    REP (i, a.size()) {
      if (i) s << ",";
      m_tb[i]->dump(s, a[i]);
    }
  }

  std::vector<CsvRecordWriterSptr<T, Stream> > m_tb;
};

template <typename T, typename Stream>
class CsvRecordWriterFromFunc : public CsvRecordWriter<T, Stream> {
public:
  typedef std::function<void(Stream &s, const T &a)> DumpFunc;
  CsvRecordWriterFromFunc(const std::vector<std::string> &headers,
                          DumpFunc func) {
    m_func = func;
    m_headers = headers;
  }
  virtual std::vector<std::string> get_headers() const override {
    return m_headers;
  }

  virtual void dump(Stream &s, const T &a) const override { m_func(s, a); }

private:
  std::vector<std::string> m_headers;
  DumpFunc m_func;
};

#if OPA_CPP14 == 1
template <typename Stream, typename... Args>
class CsvTupleRecordWriter
    : public CsvRecordWriter<std::tuple<Args...>, Stream> {
public:
  typedef std::tuple<Args...> CurType;
  CsvTupleRecordWriter(const std::vector<std::string> &headers) {
    m_headers = headers;
  }
  virtual std::vector<std::string> get_headers() const override {
    return m_headers;
  }

  virtual void dump(Stream &s, const CurType &a) const override {
    auto vec_str = opa::utils::tuple_to_string_vec(a);
    REP (i, vec_str.size()) {
      if (i) s << ",";
      s << vec_str[i];
    }
  }

private:
  std::vector<std::string> m_headers;
};
#endif

template <typename T, typename Stream> class CsvWriter {
public:
  CsvWriter() {}
  CsvWriter(CsvRecordWriterSptr<T, Stream> writer, Stream *stream) {
    init(writer, stream);
  }
  void init(CsvRecordWriterSptr<T, Stream> writer, Stream *stream) {
    m_stream = stream;
    m_writer = writer;
    *m_stream << opa::utils::Join(",", m_writer->get_headers()) << "\n";
  }

  void add(const T &a) {
    m_writer->dump(*m_stream, a);
    *m_stream << "\n";
  }

private:
  Stream *m_stream;
  CsvRecordWriterSptr<T, Stream> m_writer;
};
OPA_DECL_SPTR_TMPL2(CsvWriter, CsvWriterSptr);

template <typename Stream>
class AutoCsvRecordWriter : public CsvRecordWriter<std::string, Stream> {
public:
  std::vector<std::string> headers;
  AutoCsvRecordWriter(const std::vector<std::string> &headers) {
    this->headers = headers;
  }

  virtual std::vector<std::string> get_headers() const override {
    return headers;
  }
  virtual void dump(Stream &s, const std::string &a) const override { s << a; }
};
OPA_DECL_SPTR_TMPL(AutoCsvRecordWriter, AutoCsvRecordWriterSptr);

template <typename Stream>
class AutoCsvWriter : public CsvWriter<std::string, Stream> {
public:
  AutoCsvWriter() {}
  AutoCsvWriter(Stream *stream, const std::vector<std::string> &headers) {
    m_ustream.reset(stream);
    CsvWriter<std::string, Stream>::init(
      (new AutoCsvRecordWriter<Stream>(headers))->get_sptr(), stream);
  }

  template <typename... Args> void push(const Args &... args) {
    std::string res = opa::utils::Join(",", args...);
    this->add(res);
  }

private:
  UPTR(Stream) m_ustream;
};
OPA_DECL_SPTR_TMPL(AutoCsvWriter, AutoCsvWriterSptr);

AutoCsvWriterSptr<std::ofstream>
SimpleFileWriter(StringRef filename, const std::vector<std::string> &headers);

class CsvReaderBase {
public:
  virtual ~CsvReaderBase() {}
  virtual const std::vector<opa::stolen::StringRef> get() = 0;
  virtual bool is_end() const = 0;
  virtual bool has_more() const = 0;
};

template <typename Stream> class CsvReader : public CsvReaderBase {
public:
  CsvReader() {}
  CsvReader(Stream *stream, bool has_headers = true, char sep=',') {
    init(stream, has_headers, sep);
  }
  void init(Stream *stream, bool has_headers = true, char sep=',') {
    m_stream = stream;
    m_has_headers = has_headers;
    m_sep = sep;

    set_next_line();
    if (m_has_headers) {
      get(); // do one more round to retrieve correct data
    }
  }

  const std::vector<opa::stolen::StringRef> get() {
    m_fields.clear();
    if (is_end()) return m_fields;

    int pos = 0;
    std::string &cur = m_cur_line[m_p ^ 1];
    while (pos < cur.size()) {
      int nxt = cur.find(m_sep, pos);
      if (nxt == std::string::npos)
        nxt = cur.size();
      else
        cur[nxt] = 0;
      m_fields.pb(opa::stolen::StringRef(cur.data() + pos, nxt - pos));
      pos = nxt + 1;
    }
    if (m_num_fields == -1) m_num_fields = m_fields.size();
    OPA_CHECK0(m_num_fields == -1 || m_fields.size() == 0 ||
               m_fields.size() == m_num_fields);

    set_next_line();

    return m_fields;
  }

  virtual bool is_end() const { return m_is_end; }
  virtual bool has_more() const { return !is_end(); }

private:
  void set_next_line() {
    m_cur_line[m_p].clear();

    while (true) {
      if (m_stream->eof()) {
        m_is_end = true;
        break;
      }
      getline(*m_stream, m_cur_line[m_p]);
      if (m_cur_line[m_p].size() > 0) break;
    }
    m_p ^= 1;
  }

  Stream *m_stream;
  bool m_has_headers;
  int m_num_fields = -1;
  // TODO: optimize, go no copy
  int m_p = 0;
  bool m_is_end = false;
  std::string m_cur_line[2];
  std::vector<opa::stolen::StringRef> m_fields;
  char m_sep;
};

template <typename... Args> class CsvFieldReader {
public:
  typedef std::vector<opa::stolen::StringRef> FieldArray;
  typedef std::tuple<Args...> OType;
  CsvFieldReader() {}
  CsvFieldReader(CsvReaderBase *reader) { init(reader); }
  void init(CsvReaderBase *reader) {
    m_is_end = false;
    m_reader = reader;
  }

  OType next() {
    const FieldArray &fields = m_reader->get();
    OPA_CHECK0(fields.size() == sizeof...(Args));
    return opa::utils::string_vec_to_tuple<Args...>(fields);
  }

  bool has() const { return m_reader->has_more(); }

  std::vector<OType> get_all() {
    std::vector<OType> res;
    while (has()) res.pb(next());
    return res;
  }

private:
  bool m_is_end;
  CsvReaderBase *m_reader;
};
#endif

OPA_NAMESPACE_END(opa, utils)
