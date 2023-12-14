#ifndef SLEDGEHAMR_SLEDGEHAMR_UTILS_H_
#define SLEDGEHAMR_SLEDGEHAMR_UTILS_H_

#include <chrono>

namespace sledgehamr{
namespace utils{

/* @brief for loop constexpr approximation.
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
 * @param   state   Data from which Laplacian is to be calculated.
 * @param   i       i-th center cell.
 * @param   j       j-th center cell.
 * @param   k       k-th center cell.
 * @param   c       Scalar component.
 * @param   dx2     squared grid spacing.
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

/** @brief Calculates the finite-difference gradient for a field.
 * @param   state   Data from which gradient is to be calculated.
 * @param   i       i-th center cell.
 * @param   j       j-th center cell.
 * @param   k       k-th center cell.
 * @param   c       Scalar component.
 * @param   dx      Grid spacing.
 * @param   axis    Axis along which the gradient is to be calculated. Valid
 *                  values are 'x', 'y', and 'z'.
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

/** @brief Checks if a value is a power of two.
 */
static bool IsPowerOfTwo(int val) {
    return (val > 0) && ((val & (val - 1)) == 0);
}

/** @brief Checks if two values are approximately equal within a given accuracy.
 */
static bool ApproxEqual(double a, double b, double eps = 1e-8) {
    return (fabs(a - b) < a*eps);
}

enum ErrorState {
    ERROR = 0,
    OK = 1,
    WARNING = 2
};

/** @brief Prints the value of a parameter to the display.
 */
template <typename T>
static void PrintParamState(std::string param_name, T val, std::string state) {
    amrex::Print() << param_name << " = " << val << " : " << state << std::endl;
}

/** @brief Prints whether the value of a parameter is okay.
 */
template <typename T>
static void AssessParam(ErrorState validity, std::string param_name, T val,
                std::string errror_msg, std::string warning_msg, int& nerrors,
                bool do_thorough_checks) {
    switch (validity) {
        case OK:
            if (do_thorough_checks) {
                PrintParamState(param_name, val, "OK");
            }
            break;
        case WARNING:
            PrintParamState(param_name, val, "WARNING: " + warning_msg);
            break;
        case ERROR:
            PrintParamState(param_name, val, "ERROR: " + errror_msg);
            nerrors++;
            break;
    }
}

template <typename T>
static void AssessParamOK(std::string param_name, T val,
                          bool do_thorough_checks) {
    int tmp = 0;
    AssessParam(ErrorState::OK, param_name, val, "", "", tmp,
                do_thorough_checks);
}

}; // namespace utils
}; // namespace sledgehamr

#endif // SLEDGEHAMR_SLEDGEHAMR_UTILS
