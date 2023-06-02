#ifndef SLEDGEHAMR_PROJECTION_H_
#define SLEDGEHAMR_PROJECTION_H_

#include <AMReX_AmrCore.H>

#include "sledgehamr.h"

namespace sledgehamr {

class Sledgehamr;

typedef std::function<double(amrex::Array4<amrex::Real const> const&, const int,
                             const int, const int, const int, const double,
                             const double, const double)> projection_fct;

/** @brief TODO
 */
class Projection {
  public:
    Projection(projection_fct function, std::string identification,
               int projection_mode = 0)
        : fct{function}, ident{identification}, mode{projection_mode} {};

    void Compute(const int id, const hid_t file_id, Sledgehamr* sim);

    projection_fct fct;
    std::string ident = "None";
    int mode = 0;

  private:
    void Add(const int i, const int j, const int k, double* projection,
             int* n_projection, const int ratio, const int dimN,
             const double val);
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PROJECTION_H_
