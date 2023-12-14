#ifndef SLEDGEHAMR_PROJECTION_H_
#define SLEDGEHAMR_PROJECTION_H_

#include <hdf5.h>

#include <AMReX_AmrCore.H>

#include "sledgehamr.h"

namespace sledgehamr {

class Sledgehamr;

/** @brief Function typdef that computes the quantity to project.
 */
typedef std::function<double(amrex::Array4<amrex::Real const> const&, const int,
        const int, const int, const int, const double, const double,
        const double, const std::vector<double>&)> projection_fct;

/** @brief Computes a line-of-sight projection for an arbitrary quantity and
 *         saves it to disk. Will up-sample the 2D projection to a finer level.
 */
class Projection {
  public:
    /** @brief Gather metadata.
     * @param   function        Computes quantity to project.
     * @param   identification  Unique name of projection.
     * @param   projection_mode Projection mode. See mode for more details.
     */
    Projection(projection_fct function, std::string identification,
               int projection_mode = 0)
        : fct{function}, ident{identification}, mode{projection_mode} {};

    void Compute(const int id, const hid_t file_id, Sledgehamr* sim);

    /** @brief Function pointer to projection function.
     */
    projection_fct fct;

    /** @brief Unique name of projection.
     */
    std::string ident = "None";

    /** @brief Projection mode:
                0: Line-of-sight sum.
                1: Maximum along the line-of-sight.
     */
    int mode = 0;

  private:
    void Add(const int i, const int j, const int k, double* projection,
             int* n_projection, const int ratio, const int dimN,
             const double val);
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_PROJECTION_H_
