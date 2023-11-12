#include "location.h"

namespace sledgehamr {

void Location::SelectClosest(const int i_new, const int j_new, const int k_new,
                             const int distance_sq_new) {
    if (distance_sq_new < distance_sq) {
        i = i_new;
        j = j_new;
        k = k_new;
        distance_sq = distance_sq_new;
    }
}

void Location::SelectClosest(Location location) {
    SelectClosest(location.i, location.j, location.k, location.distance_sq);
}

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

std::vector<int> Location::Gather(const int val) {
    std::vector<int> all(amrex::ParallelDescriptor::NProcs());
    amrex::ParallelDescriptor::Gather(&val, 1, &all[0], 1,
            amrex::ParallelDescriptor::IOProcessorNumber());
    return all;
}

std::vector<int> Location::AllGather(const int val) {
    std::vector<int> all(amrex::ParallelDescriptor::NProcs());
    amrex::ParallelAllGather::AllGather(&val, 1, &all[0], MPI_COMM_WORLD);
    return all;
}


}; // namespace sledgehamr

