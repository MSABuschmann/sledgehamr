#ifndef SLEDGEHAMR_SLEDGEHAMR_UTILS_H_
#define SLEDGEHAMR_SLEDGEHAMR_UTILS_H_

#include <chrono>

namespace sledgehamr{
    namespace utils{

/* brief TODO
 */
template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f) {
    if constexpr (Start < End) {
        f(std::integral_constant<decltype(Start), Start>());
        constexpr_for<Start + Inc, End, Inc>(f);
    }
}

typedef std::chrono::steady_clock::time_point sctp;

/* brief TODO
 */
static sctp StartTimer() {
    amrex::ParallelDescriptor::Barrier();
    return std::chrono::steady_clock::now();
}

/* brief TODO
 */
static double DurationSeconds(sctp start) {
    amrex::ParallelDescriptor::Barrier();
    sctp stop = std::chrono::steady_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(
            stop - start);
    return static_cast<double>(duration.count())/1e3;
}

/** @brief Calculates the Laplacian at a given order. Only 0th, 1st and 2nd
 *         order are currently implemented.
 * @param   state_fab   Data from which Laplacian is to be calculated.
 * @param   i           i-th center cell.
 * @param   j           j-th center cell.
 * @param   k           k-th center cell.
 * @param   c           Scalar component.
 * @param   dx2         Squared grid spacing.
 */
template<int> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian(amrex::Array4<amrex::Real const> const& state_fab, const int i,
                 const int j, const int k, const int c, const double dx2);

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<0>(amrex::Array4<amrex::Real const> const& state_fab,
                    const int i, const int j, const int k, const int c,
                    const double dx2) {
    return 0;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<1>(amrex::Array4<amrex::Real const> const& state_fab,
                    const int i, const int j, const int k, const int c,
                    const double dx2) {
    return (state_fab(i + 1, j, k, c) + state_fab(i - 1, j, k, c) +
            state_fab(i, j + 1, k, c) + state_fab(i, j - 1, k, c) +
            state_fab(i, j, k + 1, c) + state_fab(i, j, k - 1, c) -
            6.*state_fab(i,j,k,c)) / dx2;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<2>(amrex::Array4<amrex::Real const> const& state_fab,
                    const int i, const int j, const int k, const int c,
                    const double dx2) {
    return ( - (state_fab(i + 2, j, k, c) + state_fab(i - 2, j, k, c) +
                state_fab(i, j + 2, k, c) + state_fab(i, j - 2, k, c) +
                state_fab(i, j, k + 2, c) + state_fab(i, j, k - 2, c) )
         + 16.*(state_fab(i + 1, j, k, c) + state_fab(i - 1, j, k, c) +
                state_fab(i, j + 1, k, c) + state_fab(i, j - 1, k, c) +
                state_fab(i, j, k + 1, c) + state_fab(i, j, k - 1, c) )
         - 90.* state_fab(i,j,k,c) ) / (12.*dx2);
};

/** @brief Returns string of ordinal number suffix for small positive numbers.
 * @param   num Number.
 */
static std::string OrdinalNumberSuffix(int num) {
    if (num == 1) return "st";
    if (num == 2) return "nd";
    if (num == 3) return "rd";
    return "th";
}

/** @brief Returns string containing the name of the level.
 * @param   lev             Level.
 * @param   shadow_hierachy Whether a shadow hierarchy has been used.
 */
static std::string LevelName(int lev, bool shadow_hierarchy) {
    if (lev==0 && shadow_hierarchy)
        return "shadow level";

    if ((lev==0 && !shadow_hierarchy) ||
        (lev==1 &&  shadow_hierarchy))
        return "coarse level";

    int refnum = lev - shadow_hierarchy;
    return std::to_string(refnum) + OrdinalNumberSuffix(refnum) + " refinement";
};

}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_UTILS
