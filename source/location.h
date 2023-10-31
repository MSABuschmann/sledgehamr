#ifndef SLEDGEHAMR_LOCATION_H_
#define SLEDGEHAMR_LOCATION_H_

#include "sledgehamr.h"

namespace sledgehamr {

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

        int i = -1;
        int j = -1;
        int k = -1;
        int distance_sq = INT_MAX;

    private:
        static std::vector<int> Gather(const int val);
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_LOCATION_H_
