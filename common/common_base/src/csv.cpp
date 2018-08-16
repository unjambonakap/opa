#include <opa/utils/csv.h>

OPA_NAMESPACE(opa, utils)
#if OPA_PIN == 0
AutoCsvWriterSptr<std::ofstream>
SimpleFileWriter(StringRef filename, const std::vector<std::string> &headers) {
  return std::make_shared<AutoCsvWriter<std::ofstream> >(
    new std::ofstream(filename.data(), std::ofstream::binary), headers);
}
#endif

OPA_NAMESPACE_END(opa, utils)
