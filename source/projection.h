#ifndef SLEDGEHAMR_PROJECTION_H_
#define SLEDGEHAMR_PROJECTION_H_

#include <AMReX_AmrCore.H>

namespace sledgehamr {

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

    projection_fct fct;
    std::string ident = "None";
    int mode = 0;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PROJECTION_H_
