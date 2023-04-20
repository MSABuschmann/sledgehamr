#ifndef SLEDGEHAMR_LEVEL_DATA_H
#define SLEDGEHAMR_LEVEL_DATA_H

#include <AMReX_MultiFab.H>

namespace sledgehamr {

/** @brief Class that holds the MultiFab data while also keeping track of time
 *         and step numbers.
 */
class LevelData : public amrex::MultiFab {
  public:
    LevelData() : amrex::MultiFab{} {};

    /** @brief Construct data using the given grid layout.
     * @param   ba      BoxArray.
     * @param   dm      DistributionMapping.
     * @param   ncomp   Total number of scalar fields components.
     * @param   nghost  Number of ghost cells.
     * @param   time    Time.
     */
    LevelData(amrex::BoxArray ba, amrex::DistributionMapping dm, int ncomp,
              int nghost, double time=0)
        : amrex::MultiFab{ba, dm, ncomp, nghost}, t{time} {};

    using amrex::MultiFab::MultiFab;

    /** @brief Overload amrex::MultiFab::define function to include time.
     * @param   ba      BoxArray.
     * @param   dm      DistributionMapping.
     * @param   ncomp   Total number of scalar components.
     * @param   nghost  Number of ghost cells.
     * @param   time    Time.
     */
    void define(amrex::BoxArray ba, amrex::DistributionMapping dm, int ncomp,
                int nghost, double time) {
        define(ba, dm, ncomp, nghost);
        t = time;
    }

    using amrex::MultiFab::define;

    /** @brief Static method that returns vector of times from a given LevelData
     *         vector.
     * @param   mfs Vector of pointers to MultiFab objects. MultiFab needs to be
     *              castable to LevelData.
     * @return  Vector of times. Needs to be amrex::Vector to be usable with
     *          e.g. amrex::FillPatchTwoLevels.
     */
    static amrex::Vector<double> getTimes(std::vector<amrex::MultiFab*>& mfs) {
        amrex::Vector<double> times;

        for (amrex::MultiFab* mf : mfs) {
            times.push_back( static_cast<LevelData*>(mf)->t );
        }

        return times;
    };

    /** @brief Time corresponding to amrex::MultiFab data.
     */
    double t = 0.;

    /** @brief How often this level has been advanced.
     */
    int istep = 0;

    /** @brief Flag as to whether the MutliFab has truncation errors encoded
     *         into it.
     */
    bool contains_truncation_errors = false;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LEVEL_DATA_H_
