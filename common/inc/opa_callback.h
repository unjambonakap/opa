#pragma once

namespace opa {
class OpaCallback {
public:
  virtual ~OpaCallback() {}
  virtual void run() {}
};

}
