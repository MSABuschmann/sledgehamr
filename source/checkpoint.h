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
    void Read(std::string folder);
    void UpdateOutputModules(std::string folder);
    void UpdateOutputModules(std::string prefix, int id);

  private:
    static void GotoNextLine(std::istream& is);
    void ChangeNGhost(int new_nghost);
    void RegridCoarse();
    void UpdateLevels(std::string filename);
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_CHECKPOINT_H_
