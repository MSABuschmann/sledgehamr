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
    bool ReadHeader(std::string folder);
    void UpdateOutputModules(std::string folder);
    void UpdateOutputModules(std::string prefix, int id);

    double GetTime() const {
        return time;
    }

  private:
    static void GotoNextLine(std::istream& is);
    void ChangeNGhost(int new_nghost);
    void RegridCoarse();
    void UpdateLevels(std::string folder);

    std::string GetHeaderName(std::string folder) const {
        return folder + "/Meta.hdf5";
    }

    Sledgehamr* sim;

    double time;
    int MPIranks;
    int finest_level;
    int dim0;
    int nghost;
    int nscalars;
    int noutput;
    int npredefoutput;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_CHECKPOINT_H_
