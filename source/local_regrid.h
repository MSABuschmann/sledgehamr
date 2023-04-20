#ifndef SLEDGEHAMR_LOCAL_REGRID_H_
#define SLEDGEHAMR_LOCAL_REGRID_H_

#include "sledgehamr.h"
#include "unique_layout.h"

namespace sledgehamr {

class Sledgehamr;
class UniqueLayout;

class LocalRegrid {
  public:
    LocalRegrid(Sledgehamr* owner);

    bool AttemptRegrid(const int lev);

    std::vector< std::vector<int> > comm_matrix;

  private:
    bool DoAttemptRegrid(const int lev);
    void ParseInput();
    void CreateCommMatrix();
    int DetermineNewBoxArray(const int lev);

    double volume_threshold_strong = 1.1;
    double volume_threshold_weak = 1.05;
    int veto_level = -1;
    bool force_global_regrid_at_restart = 0;
    int n_error_buf = 0;

    std::vector<long long> numPts;
    std::vector<bool> no_local_regrid;

    std::vector< std::vector< std::unique_ptr<UniqueLayout> > > layouts;

    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LOCAL_REGRID_H_
