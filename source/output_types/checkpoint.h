#ifndef SLEDGEHAMR_CHECKPOINT_H_
#define SLEDGEHAMR_CHECKPOINT_H_

#include "sledgehamr.h"

namespace sledgehamr {

class Checkpoint {
  public:
    Checkpoint(Sledgehamr* owner, std::string chk_folder)
        : sim(owner), folder(chk_folder) {};

    void Write();
    void Read();
    bool ReadHeader();
    void UpdateOutputModules();
    void Delete();

    double GetTime() const {
        return time;
    }

  private:
    static void GotoNextLine(std::istream& is);
    void UpdateLevels();

    std::string GetHeaderName() const {
        return folder + "/Meta.hdf5";
    }

    std::string GetBoxArrayName() const {
        return folder + "/BoxArrays";
    }

    std::string GetLevelDirName(const int lev) const {
        return folder + "/Level_" + std::to_string(lev);
    }

    Sledgehamr* sim;
    std::string folder;

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
