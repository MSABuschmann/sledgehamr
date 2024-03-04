#ifndef SLEDGEHAMR_LOCAL_REGRID_H_
#define SLEDGEHAMR_LOCAL_REGRID_H_

#include "boost/multi_array.hpp"

#include "sledgehamr.h"
#include "unique_layout.h"
#include "location.h"

namespace sledgehamr {

class Sledgehamr;
class UniqueLayout;
class Location;

/** @brief Class to perform a local regrid (if possible).
 */
class LocalRegrid {
  public:
    LocalRegrid(Sledgehamr* owner);

    bool AttemptRegrid(const int lev);
    void DidGlobalRegrid(const int lev);

    void InitializeLayout(const int max_lev);
    void ClearLayout();
    void FinalizeLayout(const int lev);
    void WrapIndices(const int lev);

    void JoinBoxArrays(const int lev, amrex::BoxArray& ba);
    void AddBoxes(const int lev, amrex::BoxArray& ba);
    void FixNesting(const int lev);
    void AddToLayout(const int lev, const int thread, const int i, const int j,
                     const int k);

    /** @brief Flag that will be checked the TimeStepper module to force a
     *         global regrid.
     */
    std::vector<bool> do_global_regrid;

    /** @brief MPI communication matrix. Look-up table used by the
     *         UniqueLayout class.
     */
    std::vector< std::vector<int> > comm_matrix;

  private:
    /* @brief Enum with possible course of actions after a veto.
     */
    enum VetoResult {
        DoGlobalRegrid = 0,
        DoNoRegrid = 1,
        DoLocalRegrid = 2,
    };

    bool Prechecks(const int lev);
    void InitializeLocalRegrid();
    void ParseInput();
    void CreateCommMatrix();

    bool DoAttemptRegrid(const int lev);
    void DetermineAllBoxArrays(const int lev);
    double DetermineNewBoxArray(const int lev);
    void FixAllNesting();
    void JoinAllBoxArrays(std::vector<amrex::BoxArray>& box_arrays);
    void AddAllBoxes(std::vector<amrex::BoxArray>& box_arrays);
    amrex::BoxArray WrapBoxArray(amrex::BoxArray& ba, int N);

    inline int GetBoxCoarseFineBorders(
        const amrex::Box& tilebox, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int lev,
        boost::multi_array<bool, 3>& border);

    inline void TagAndMeasure(
        const amrex::Dim3& lo, const amrex::Dim3& hi, int remaining,
        const amrex::Array4<char>& tag_arr, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int lev, const int ibff,
        const double bff, boost::multi_array<bool, 3>& border,
        std::vector<Location>& closest_locations, const double threshold,
        const int omp_thread_num);

    inline void CheckBorders(
        const amrex::IntVect& ci, const amrex::IntVect& c0,
        const amrex::IntVect& c1, const int ibff, const double bff,
        int remaining, const int lev, boost::multi_array<bool, 3>& border,
        std::vector<Location>& closest_locations, const double threshold,
        const int omp_thread_num);

    bool CheckThresholds(const int lev, amrex::BoxArray& box_array);
    void ComputeLatestPossibleRegridTime(const int l, const int lev);
    bool CheckForVeto(const int lev,
                      std::vector<amrex::BoxArray>& box_arrays);
    VetoResult DealWithVeto(const int lev);

    /** @brief Strong threshold that decide whether we want to do a local or
     *         global regrid.
     */
    double volume_threshold_accumulated = 1.1;

    /** @brief Weak threshold that decide whether we want to do a local or
     *         global regrid.
     */
    double volume_threshold_single      = 1.05;

    /** @brief Lowest level that would rather do a global regrid than a local.
     */
    int veto_level = -1;

    /** @brief If flag set to true we will skip the local regrid once.
     */
    bool force_global_regrid_at_restart = 0;

    /** @brief How often we are allowed to do a local regrid in a row.
     */
    int max_local_regrids = 10;

    /** @brief Counts number of consecutive local regrids.
     */
    int nregrids = 0;

    /** @brief Size of error buffer.
     */
    int n_error_buf = 1;

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

    /** @brief Contains the time by which we should have done a regrid.
     */
    std::vector<double> latest_possible_regrid_time;

    /** @brief Minimum distance of a tagged cell to a coarse/fine boundary for
     *         each OpenMP thread.
     */
    std::vector<double> min_distance;

    /** @brief Vector of UniqueLayouts for each level and each core.
     */
    std::vector< std::vector< std::unique_ptr<UniqueLayout> > > layouts;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LOCAL_REGRID_H_
