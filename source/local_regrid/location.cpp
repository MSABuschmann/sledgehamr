#include "location.h"

namespace sledgehamr {

/** @brief Overrides current location if provided new location is closer.
 * @param   i_new           New i-th location.
 * @param   j_new           New j-th location.
 * @param   k_new           New k-th location.
 * @param   distance_sq_new New distance square. 
 */
void Location::SelectClosest(const int i_new, const int j_new, const int k_new,
                             const int distance_sq_new) {
    if (distance_sq_new < distance_sq) {
        i = i_new;
        j = j_new;
        k = k_new;
        distance_sq = distance_sq_new;
    }
}

/** @brief Overrides current location if provided new location is closer.
 * @param   location    New location.
 */
void Location::SelectClosest(Location location) {
    SelectClosest(location.i, location.j, location.k, location.distance_sq);
}

/** @brief Finds the location with the globally shortest distance after doing
 *         OpenMP and MPI reduction.
 * @param   locations   Location at each OpenMP thread.
 * @return  Closest location.
 */
Location Location::FindClosestGlobally(
        const std::vector<Location>& locations) {
    Location closest;

    // OpenMP reduction
    for (Location location : locations)
        closest.SelectClosest(location);

    // MPI reducation
    std::vector<int> all_i = AllGather(closest.i);
    std::vector<int> all_j = AllGather(closest.j);
    std::vector<int> all_k = AllGather(closest.k);
    std::vector<int> all_d = AllGather(closest.distance_sq);

    for (int l = 0; l < all_d.size(); ++l) {
        closest.SelectClosest(all_i[l], all_j[l], all_k[l], all_d[l]);
    }

    return closest;
}

/** @brief MPI Gather.
 * @param   val Value to gather.
 * @return Gathered values.
 */
std::vector<int> Location::Gather(const int val) {
    std::vector<int> all(amrex::ParallelDescriptor::NProcs());
    amrex::ParallelDescriptor::Gather(&val, 1, &all[0], 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
    return all;
}

/** @brief MPI AllGather.
 * @param   val Value to gather.
 * @return Gathered values.
 */
std::vector<int> Location::AllGather(const int val) {
    std::vector<int> all(amrex::ParallelDescriptor::NProcs());
    amrex::ParallelAllGather::AllGather(&val, 1, &all[0], MPI_COMM_WORLD);
    return all;
}

}; // namespace sledgehamr

