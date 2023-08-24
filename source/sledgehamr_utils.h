#ifndef SLEDGEHAMR_SLEDGEHAMR_UTILS_H_
#define SLEDGEHAMR_SLEDGEHAMR_UTILS_H_

#include <chrono>

namespace sledgehamr{
namespace utils{

/* @brief for constexpr approximation.
 */
template <auto Start, auto End, auto Inc, class F>
constexpr void constexpr_for(F&& f) {
    if constexpr (Start < End) {
        f(std::integral_constant<decltype(Start), Start>());
        constexpr_for<Start + Inc, End, Inc>(f);
    }
}

typedef std::chrono::steady_clock::time_point sctp;

/* @brief Starts a timer.
 * @return Timer.
 */
static sctp StartTimer() {
    amrex::ParallelDescriptor::Barrier();
    return std::chrono::steady_clock::now();
}

/* @brief Computes elapsed time since start of a timer in seconds.
 * @param   start   Timer.
 * @return Elapsed time.
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
double Laplacian(amrex::Array4<amrex::Real const> const& state, const int i,
                 const int j, const int k, const int c, const double dx2);

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<0>(amrex::Array4<amrex::Real const> const& state, const int i,
                    const int j, const int k, const int c, const double dx2) {
    return 0;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<1>(amrex::Array4<amrex::Real const> const& state, const int i,
                    const int j, const int k, const int c, const double dx2) {
    return (state(i+1,j,  k,  c) + state(i-1,j,  k,  c)
          + state(i,  j+1,k,  c) + state(i,  j-1,k,  c)
          + state(i,  j,  k+1,c) + state(i,  j,  k-1,c)
          - 6.*state(i,j,k,c)) / dx2;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<2>(amrex::Array4<amrex::Real const> const& state, const int i,
                    const int j, const int k, const int c, const double dx2) {
    return ( - (state(i+2,j,  k,  c) + state(i-2,j,  k,  c) +
                state(i,  j+2,k,  c) + state(i,  j-2,k,  c) +
                state(i,  j,  k+2,c) + state(i,  j,  k-2,c) )
         + 16.*(state(i+1,j,  k,  c) + state(i-1,j,  k,  c) +
                state(i,  j+1,k,  c) + state(i,  j-1,k,  c) +
                state(i,  j,  k+1,c) + state(i,  j,  k-1,c) )
         - 90.* state(i,  j,  k,  c) ) / (12.*dx2);
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Laplacian<3>(amrex::Array4<amrex::Real const> const& state, const int i,
                    const int j, const int k, const int c, const double dx2) {
    return (2.*(state(i+3,j,  k,  c) + state(i-3,j,  k,  c) +
                state(i,  j+3,k,  c) + state(i,  j-3,k,  c) +
                state(i,  j,  k+3,c) + state(i,  j,  k-3,c) )
          -27.*(state(i+2,j,  k,  c) + state(i-2,j,  k,  c) +
                state(i,  j+2,k,  c) + state(i,  j-2,k,  c) +
                state(i,  j,  k+2,c) + state(i,  j,  k-2,c) )
         +270.*(state(i+1,j,  k,  c) + state(i-1,j,  k,  c) +
                state(i,  j+1,k,  c) + state(i,  j-1,k,  c) +
                state(i,  j,  k+1,c) + state(i,  j,  k-1,c) )
         -1470.*state(i,  j,  k,  c) ) / (180.*dx2);
};

/** @brief TODO
 */
template<int> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Gradient(amrex::Array4<amrex::Real const> const& state, const int i,
                const int j, const int k, const int c, const double dx,
                const char axis);

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Gradient<0>(amrex::Array4<amrex::Real const> const& state, const int i,
                   const int j, const int k, const int c, const double dx,
                   const char axis) {
    return 0;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Gradient<1>(amrex::Array4<amrex::Real const> const& state, const int i,
                   const int j, const int k, const int c, const double dx,
                   const char axis) {
    double result = 0;

    switch (axis) {
        case 'x':
            result = (state(i+1,j,  k,  c) - state(i-1,j,  k,  c)) / (2.*dx);
            break;
        case 'y':
            result = (state(i,  j+1,k,  c) - state(i,  j-1,k,  c)) / (2.*dx);
            break;
        case 'z':
            result = (state(i,  j,  k+1,c) - state(i,  j,  k-1,c)) / (2.*dx);
            break;
    }

    return result;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Gradient<2>(amrex::Array4<amrex::Real const> const& state, const int i,
                   const int j, const int k, const int c, const double dx,
                   const char axis) {
    double result = 0;

    switch (axis) {
        case 'x':
            result = (-state(i+2,j,  k,  c) + 8.*state(i+1,j,  k,  c)
                      +state(i-2,j,  k,  c) - 8.*state(i-1,j,  k,  c))/(12.*dx);
            break;
        case 'y':
            result = (-state(i,  j+2,k,  c) + 8.*state(i,  j+1,k,  c)
                      +state(i,  j-2,k,  c) - 8.*state(i,  j-1,k,  c))/(12.*dx);
            break;
        case 'z':
            result = (-state(i,  j,  k+2,c) + 8.*state(i,  j,  k+1,c)
                      +state(i,  j,  k-2,c) - 8.*state(i,  j,  k-1,c))/(12.*dx);
            break;
    }

    return result;
};

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double Gradient<3>(amrex::Array4<amrex::Real const> const& state, const int i,
                   const int j, const int k, const int c, const double dx,
                   const char axis) {
    double result = 0;

    switch (axis) {
        case 'x':
            result = (+    state(i+3,j,k,c) -    state(i-3,j,k,c)
                      - 9.*state(i+2,j,k,c) + 9.*state(i-2,j,k,c)
                      +45.*state(i+1,j,k,c) -45.*state(i-1,j,k,c)) / (60.*dx);
            break;
        case 'y':
            result = (+    state(i,j+3,k,c) -    state(i,j-3,k,c)
                      - 9.*state(i,j+2,k,c) + 9.*state(i,j-2,k,c)
                      +45.*state(i,j+1,k,c) -45.*state(i,j-1,k,c)) / (60.*dx);
            break;
        case 'z':
            result = (+    state(i,j,k+3,c) -    state(i,j,k-3,c)
                      - 9.*state(i,j,k+2,c) + 9.*state(i,j,k-2,c)
                      +45.*state(i,j,k+1,c) -45.*state(i,j,k-1,c)) / (60.*dx);
            break;
    }

    return result;
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
static std::string LevelName(int lev) {
    switch (lev) {
        case -1:
            return "shadow level";
        case 0:
            return "coarse level";
        default:
            return std::to_string(lev) + OrdinalNumberSuffix(lev)
                    + " refinement";
    }
}

static bool IsPowerOfTwo(int val) {
    return (val > 0) && ((val & (val - 1)) == 0);
}

}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_UTILS
