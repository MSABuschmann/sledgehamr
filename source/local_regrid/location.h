#ifndef SLEDGEHAMR_LOCATION_H_
#define SLEDGEHAMR_LOCATION_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Keeps track of a location and a distance and can do MPI and OpenMP
 *         reductions to determine the one with the smallest distance.
 */
class Location {
    public:
        Location() {};
        Location(const int i_new, const int j_new, const int k_new,
                 const int distance_sq_new)
            : i(i_new), j(j_new), k(k_new), distance_sq(distance_sq_new) {};

        void SelectClosest(const int i_new, const int j_new, const int k_new,
                           const int distance_sq_new);
        void SelectClosest(const Location location);

        static Location FindClosestGlobally(
                const std::vector<Location>& locations);

        /** @brief Locations i-th index.
         */
        int i = -1;

        /** @brief Locations j-th index.
         */
        int j = -1;

        /** @brief Locations k-th index.
         */
        int k = -1;

        /** @brief Distance square of location.
         */
        int distance_sq = INT_MAX;

    private:
        static std::vector<int> Gather(const int val);
        static std::vector<int> AllGather(const int val);
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LOCATION_H_
