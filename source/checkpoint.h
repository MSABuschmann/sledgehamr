#ifndef SLEDGEHAMR_CHECKPOINT_H_
#define SLEDGEHAMR_CHECKPOINT_H_

#include "sledgehamr.h"

namespace sledgehamr {

class Checkpoint {
  public:
    Checkpoint(Sledgehamr* owner) {
        sim = owner;
    };

    void Write(std::string prefix);
    void Read(std::string prefix, int id);

  private:
    static void GotoNextLine(std::istream& is);

    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_CHECKPOINT_H_
