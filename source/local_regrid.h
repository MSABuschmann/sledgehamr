#ifndef SLEDGEHAMR_LOCAL_REGRID_H_
#define SLEDGEHAMR_LOCAL_REGRID_H_

#include "sledgehamr.h"
#include "unique_layout.h"

namespace sledgehamr {

class Sledgehamr;
class UniqueLayout;

/** @brief Class to perform a local regrid (if possible).
 */
class LocalRegrid {
  public:
    LocalRegrid(Sledgehamr* owner);

    /** @brief Perform a local regrid on level lev after checking various
     *         criteria.
     */
    bool AttemptRegrid(const int lev);

    /** @brief To be called if a global regrid on level lev has been performed.
     *         Needed in case this LocalRegrid module invoked the global regrid
     *         such that we can reset the relevant flags.
     */
    void DidGlobalRegrid(const int lev);

    /** @brief Flag that will be checked the TimeStepper module to force a
     *         global regrid.
     */
    std::vector<bool> do_global_regrid;

    /** @brief  N<->M MPI communication matrix. Look-up table used by the
     *          UniqueLayout class.
     */
    std::vector< std::vector<int> > comm_matrix;

  private:
    bool DoAttemptRegrid(const int lev);
    void ParseInput();

    /** @brief Creates comm_matrix loop-up table.
     */
    void CreateCommMatrix();

    /** @brief Pre-computes wrapped index loop-up table for periodic boundary
     *         conditions.
     */
    void WrapIndices(const int lev);

    /** @brief Determines the amrex::BoxArray that will need to be added at a
     *         given level.
     */
    double DetermineNewBoxArray(const int lev);

    /** @brief Ensures that the new BoxArray at level lev is properly nested
     *         with the lower levels by expanding the lower levels if needed.
     */
    void FixNesting(const int lev);

    /** @brief Wrapps an amrex::BoxArray across periodic boundary conditions.
     */
    amrex::BoxArray WrapBoxArray(amrex::BoxArray& ba, int N);

    /** @brief Adds an amrex::BoxArray to the existing BoxArray at level lev.
     *         Data will be allocated and filled.
     */
    void AddBoxes(const int lev, amrex::BoxArray& ba);

    /** @brief Strong and weak threshold that decide whether we want to do a
     *         local or global regrid.
     */
    double volume_threshold_strong = 1.1;
    double volume_threshold_weak = 1.05;

    /** @brief Lowest level that would rather do a global regrid than a local.
     */
    int veto_level = -1;

    /** @brief If flag set to true we will skip the local regrid once.
     */
    bool force_global_regrid_at_restart = 0;

    /** @brief How often we are allowed to do a local regrid in a row and 
     *         corresponding counter.
     */
    int max_local_regrids = 10;
    int nregrids = 0;

    /** @brief Size of error buffer.
     */
    int n_error_buf = 0;

    /** @brief Number of cells contained in each level after the last global
     *         regrid.
     */
    std::vector<long long> last_numPts;

    /** @brief If we want to skip the local regrid entirely at a given level.
     */
    std::vector<bool> no_local_regrid;

    /** @brief Look-up table of wrapped indices for periodic boundary
     *         conditions.
     */
    std::vector< std::vector<int> > wrapped_index;

    /** @brief Vector of UniqueLayouts for each level and each core.
     */
    std::vector< std::vector< std::unique_ptr<UniqueLayout> > > layouts;

    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LOCAL_REGRID_H_
