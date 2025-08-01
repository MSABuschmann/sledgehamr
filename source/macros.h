#ifndef SLEDGEHAMR_MACROS_H_
#define SLEDGEHAMR_MACROS_H_

#include "kernels.h"
#include "sledgehamr_utils.h"

/** Yes, I know this file is disgusting. This is what we get for not
 *  wanting everyone to write a bunch of complicated boilerplate every single
 *  time but also want to avoid all this function overhead for performance ...
 */

namespace sledgehamr {

/* @brief Force boost to define variadics.
 */
#define BOOST_PP_VARIADICS 1
#include <boost/preprocessor.hpp>

/** @brief Element-wise vector reduction for OpenMP.
 */
#pragma omp declare reduction(                                                 \
        vec_int_plus : std::vector<int> : std::transform(                      \
                omp_out.begin(), omp_out.end(), omp_in.begin(),                \
                    omp_out.begin(), std::plus<int>()))                        \
    initializer(omp_priv = omp_orig)

#pragma omp declare reduction(                                                 \
        vec_long_plus : std::vector<long> : std::transform(                    \
                omp_out.begin(), omp_out.end(), omp_in.begin(),                \
                    omp_out.begin(), std::plus<long>()))                       \
    initializer(omp_priv = omp_orig)

/** @brief Expand OMP #pragma statement.
 */
#define SLEDGEHAMR_DO_PRAGMA(x) _Pragma(#x)

/** @brief Declare and initialize single scalar field within project namespace.
 */
#define SLEDGEHAMR_ADD_SCALAR(field)                                           \
    static sledgehamr::ScalarField BOOST_PP_CAT(_s_, field) = {                \
        #field, l_scalar_fields, false};                                       \
    namespace Scalar {                                                         \
    constexpr int field = _Scalar::field;                                      \
    };

/** @brief Expand single scalar field.
 */
#define SLEDGEHAMR_EXPAND_SCALARS(r, data, field) SLEDGEHAMR_ADD_SCALAR(field)

/** @brief Declare and initialize single momentum field within project
 *         namespace.
 */
#define SLEDGEHAMR_ADD_MOMENTUM(field)                                         \
    static sledgehamr::ScalarField BOOST_PP_CAT(_s_, field) = {                \
        #field, l_scalar_fields, true};                                        \
    namespace Scalar {                                                         \
    constexpr int field = _Momentum::field + _Scalar::NScalarFields;           \
    };

/** @brief Expand single momnetum field.
 */
#define SLEDGEHAMR_EXPAND_MOMENTA(r, data, field) SLEDGEHAMR_ADD_MOMENTUM(field)

/** @brief Expand enum element.
 */
#define SLEDGEHAMR_SCALAR_ENUM_VALUE(r, data, elem) elem,

/** @brief Macros to create enum of scalar fields for fast and convinient
 *         component access within the project namespace.
 */
#define SLEDGEHAMR_SCALAR_ENUM(...)                                            \
    enum _Scalar {                                                             \
        BOOST_PP_SEQ_FOR_EACH(SLEDGEHAMR_SCALAR_ENUM_VALUE, _,                 \
                              BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))           \
            NScalarFields                                                      \
    };

/** @brief Same as SLEDGEHAMR_SCALAR_ENUM but for momentum fields.
 */
#define SLEDGEHAMR_MOMENTUM_ENUM(...)                                          \
    enum _Momentum {                                                           \
        BOOST_PP_SEQ_FOR_EACH(SLEDGEHAMR_SCALAR_ENUM_VALUE, _,                 \
                              BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))           \
            NMomentumFields                                                    \
    };

/** @brief Same as SLEDGEHAMR_SCALAR_ENUM but for gravitational wave tensor
 *         components.
 */
#define SLEDGEHAMR_GW_ENUM                                                     \
    enum Gw {                                                                  \
        u_xx = Scalar::NScalars,                                               \
        u_yy,                                                                  \
        u_zz,                                                                  \
        u_xy,                                                                  \
        u_xz,                                                                  \
        u_yz,                                                                  \
        du_xx,                                                                 \
        du_yy,                                                                 \
        du_zz,                                                                 \
        du_xy,                                                                 \
        du_xz,                                                                 \
        du_yz,                                                                 \
        NGwScalars                                                             \
    };

/* brief Default implementation for f(\tau) > \tau_{crit} criteria:
 *       f(\tau) = \tau. This template can be specialized for each scalar field
 *       component by the project.
 */
#define SLEDGEHAMR_TRUNCATION_MODIFIER                                         \
    template <int>                                                             \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE double TruncationModifier(        \
        const amrex::Array4<const double> &state, const int i, const int j,    \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, const double truncation_error,                        \
        const double *params) {                                                \
        return truncation_error;                                               \
    };

/* brief User-defined tagging criteria, switched off by default. Template <true>
 *       will be used by tagger and this specialisation can be implemented
 *       within the project.
 */
#define SLEDGEHAMR_TAG_CELL_FOR_REFINEMENT                                     \
    template <bool>                                                            \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool TagCellForRefinement(        \
        const amrex::Array4<const double> &state, const int i, const int j,    \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, const double *params) {                               \
        return false;                                                          \
    };

/* @brief Template declaration for procedure that computes the Rhs of the
 *        gravitational wave tensor components.
 */
#define SLEDGEHAMR_GRAVITATIONAL_WAVES_RHS                                     \
    template <bool>                                                            \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void GravitationalWavesRhs(       \
        const amrex::Array4<double> &rhs,                                      \
        const amrex::Array4<const double> &state, const int i, const int j,    \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, const double *params) {};

/* @brief Template declaration for computing the backreaction of the
 *        gravitational tensor onto our scalar fields.
 */
#define SLEDGEHAMR_GRAVITATIONAL_WAVES_BACKREACTION                            \
    template <bool>                                                            \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE void                              \
    GravitationalWavesBackreaction(                                            \
        const amrex::Array4<double> &rhs,                                      \
        const amrex::Array4<const double> &state, const int i, const int j,    \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, const double *params_scalars,                         \
        const double *params_gw) {};

/** @brief Macro to add multiple scalar fields and default template functions to
 *         the project namespace.
 */
#define SLEDGEHAMR_ADD_SCALARS(...)                                            \
    SLEDGEHAMR_SCALAR_ENUM(__VA_ARGS__)                                        \
    static std::vector<sledgehamr::ScalarField *> l_scalar_fields;             \
    BOOST_PP_SEQ_FOR_EACH(SLEDGEHAMR_EXPAND_SCALARS, _,                        \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

/** @brief Same as SLEDGEHAMR_ADD_SCALARS but for momentum fields.
 */
#define SLEDGEHAMR_ADD_CONJUGATE_MOMENTA(...)                                  \
    SLEDGEHAMR_MOMENTUM_ENUM(__VA_ARGS__)                                      \
    BOOST_PP_SEQ_FOR_EACH(SLEDGEHAMR_EXPAND_MOMENTA, _,                        \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))               \
    namespace Scalar {                                                         \
    constexpr int NScalars =                                                   \
        _Scalar::NScalarFields + _Momentum::NMomentumFields;                   \
    };                                                                         \
    SLEDGEHAMR_GW_ENUM                                                         \
    SLEDGEHAMR_TRUNCATION_MODIFIER                                             \
    SLEDGEHAMR_TAG_CELL_FOR_REFINEMENT                                         \
    SLEDGEHAMR_GRAVITATIONAL_WAVES_RHS                                         \
    SLEDGEHAMR_GRAVITATIONAL_WAVES_BACKREACTION

/** @brief Identifies cells that violate the truncation error threshold. To be
 *         run on CPU code.
 * @param   state       Current grid.
 * @param   te          Grid containing truncation errors.
 * @param   i           i-th cell index.
 * @param   j           j-th cell index.
 * @param   k           k-th cell index.
 * @param   lev         Current level.
 * @param   time        Current time.
 * @param   dt          Time step size.
 * @param   dx          Grid spacing.
 * @param   te_crit     Array containing truncation error thresholds.
 * @param   ntags_trunc Counts number of tags.
 * @param   params      User-defined parameter.
 * @return  If truncation error threshold has been exceeded.
 */
#define SLEDGEHAMR_TRUNCATION_ERROR_TAG_CPU                                    \
    template <int NScalars>                                                    \
    AMREX_FORCE_INLINE bool TruncationErrorTagCpu(                             \
        const amrex::Array4<const double> &state,                              \
        const amrex::Array4<const double> &te, const int i, const int j,       \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, std::vector<double> &te_crit, long *ntags_trunc,      \
        const std::vector<double> &params) {                                   \
        if (i % 2 != 0 || j % 2 != 0 || k % 2 != 0)                            \
            return false;                                                      \
        bool res = false;                                                      \
        sledgehamr::utils::constexpr_for<0, NScalars, 1>([&](auto n) {         \
            double mte =                                                       \
                TruncationModifier<n>(state, i, j, k, lev, time, dt, dx,       \
                                      te(i, j, k, n), params.data());          \
            if (mte >= te_crit[n]) {                                           \
                res = true;                                                    \
                ntags_trunc[n] += 8;                                           \
            }                                                                  \
        });                                                                    \
        return res;                                                            \
    };

/** @brief Same as TRUNCATION_ERROR_TAG_CPU but to be used for GPU code.
 */
#define SLEDGEHAMR_TRUNCATION_ERROR_TAG_GPU                                    \
    template <int NScalars>                                                    \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE bool TruncationErrorTagGpu(       \
        const amrex::Array4<double const> &state,                              \
        const amrex::Array4<double const> &te, const int i, const int j,       \
        const int k, const int lev, const double time, const double dt,        \
        const double dx, double *te_crit, const double *params) {              \
        if (i % 2 != 0 || j % 2 != 0 || k % 2 != 0)                            \
            return false;                                                      \
        bool res = false;                                                      \
        sledgehamr::utils::constexpr_for<0, NScalars, 1>([&](auto n) {         \
            double mte = TruncationModifier<n>(state, i, j, k, lev, time, dt,  \
                                               dx, te(i, j, k, n), params);    \
            if (mte >= te_crit[n]) {                                           \
                res = true;                                                    \
            }                                                                  \
        });                                                                    \
        return res;                                                            \
    };

/** @brief Add functions to project namespace that depend on custom template
 *         specialisations.
 */
#define SLEDGEHAMR_FINISH_SETUP                                                \
    SLEDGEHAMR_TRUNCATION_ERROR_TAG_CPU                                        \
    SLEDGEHAMR_TRUNCATION_ERROR_TAG_GPU

/** @brief Constructor of project class to initialize scalar fields
 *         automatically.
 */
#define SLEDGEHAMR_PRJ_CONSTRUCTOR(prj)                                        \
    prj() {                                                                    \
        scalar_fields = l_scalar_fields;                                       \
        amrex::Print() << "Starting " << #prj << " project..." << std::endl;   \
        amrex::Print() << "Number of field components: "                       \
                       << scalar_fields.size() << std::endl;                   \
        amrex::Print() << std::endl;                                           \
    };

/** @brief Initializes async arrays for Kreiss-Oliger meta data.
 */
#define SLEDGEHAMR_KO_LOCAL_SETUP                                              \
    amrex::Gpu::AsyncArray<double> async_dissipation_strength(                 \
        dissipation_strength.data(), dissipation_strength.size());             \
    double *l_dissipation_strength = async_dissipation_strength.data();        \
    const int l_dissipation_order = dissipation_order;                         \
    const bool l_with_dissipation = with_dissipation;

/** @brief Initializes async arrays for meta data needed by Rhs computation.
 */
#define SLEDGEHAMR_RHS_PARAMS_LOCAL_SETUP                                      \
    std::vector<double> params_rhs;                                            \
    SetParamsRhs(params_rhs, time, lev);                                       \
    amrex::Gpu::AsyncArray<double> async_params_rhs(params_rhs.data(),         \
                                                    params_rhs.size());        \
    double *l_params_rhs = async_params_rhs.data();

/** @brief Initializes async arrays for meta data needed by Rhs computation of
 *         gravitational waves.
 */
#define SLEDGEHAMR_RHS_GW_PARAMS_LOCAL_SETUP                                   \
    std::vector<double> params_gw_rhs;                                         \
    if (with_gravitational_waves) {                                            \
        SetParamsGravitationalWaveRhs(params_gw_rhs, time, lev);               \
    }                                                                          \
    amrex::Gpu::AsyncArray<double> async_params_gw_rhs(params_gw_rhs.data(),   \
                                                       params_gw_rhs.size());  \
    double *l_params_gw_rhs = async_params_gw_rhs.data();

/** @brief Computes Rhs for the entire level.
 * @param   rhs_mf      Container to fill the Rhs with.
 * @param   state_mf    Current grid.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   dt          Time step size.
 * @param   dx          Grid spacing.
 */
#define SLEDGEHAMR_PRJ_FILL_RHS                                                \
    virtual void FillRhs(amrex::MultiFab &rhs_mf,                              \
                         const amrex::MultiFab &state_mf, const double time,   \
                         const int lev, const double dt, const double dx)      \
        override {                                                             \
        performance_monitor->Start(performance_monitor->idx_rhs, lev);         \
        SLEDGEHAMR_KO_LOCAL_SETUP                                              \
        SLEDGEHAMR_RHS_PARAMS_LOCAL_SETUP                                      \
        SLEDGEHAMR_RHS_GW_PARAMS_LOCAL_SETUP                                   \
        SLEDGEHAMR_DO_PRAGMA(                                                  \
            omp parallel if (amrex::Gpu::notInLaunchRegion()))                 \
        for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU());               \
             mfi.isValid(); ++mfi) {                                           \
            const amrex::Box &bx = mfi.tilebox();                              \
            const amrex::Array4<double> &rhs_fab = rhs_mf.array(mfi);          \
            const amrex::Array4<double const> &state_fab =                     \
                state_mf.array(mfi);                                           \
            if (with_gravitational_waves) {                                    \
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,      \
                                                            int k) noexcept {  \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs);                                         \
                    GravitationalWavesRhs<true>(rhs_fab, state_fab, i, j, k,   \
                                                lev, time, dt, dx,             \
                                                l_params_gw_rhs);              \
                    GravitationalWavesBackreaction<true>(                      \
                        rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs, l_params_gw_rhs);                        \
                    switch (l_with_dissipation * l_dissipation_order) {        \
                    case 2:                                                    \
                        sledgehamr::utils::constexpr_for<0, Gw::NGwScalars,    \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    2>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    case 3:                                                    \
                        sledgehamr::utils::constexpr_for<0, Gw::NGwScalars,    \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    3>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    }                                                          \
                });                                                            \
            } else {                                                           \
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,      \
                                                            int k) noexcept {  \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs);                                         \
                    switch (l_with_dissipation * l_dissipation_order) {        \
                    case 2:                                                    \
                        sledgehamr::utils::constexpr_for<0, Scalar::NScalars,  \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    2>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    case 3:                                                    \
                        sledgehamr::utils::constexpr_for<0, Scalar::NScalars,  \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    3>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    }                                                          \
                });                                                            \
            }                                                                  \
        }                                                                      \
        performance_monitor->Stop(performance_monitor->idx_rhs, lev);          \
    };

/** @brief Computes Rhs for the entire level and adds it weighted to the
 *         existing values stored in the rhs_mf, rhs_mf += weight*rhs, instead
 *         of overwriting them as SLEDGEHAMR_PRJ_FILL_RHS does.
 * @param   rhs_mf      Container to fill the Rhs with.
 * @param   state_mf    Current grid.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   dt          Time step size.
 * @param   dx          Grid spacing.
 * @param   weight      Relative weight between values of rhs_mf and rhs.
 */
#define SLEDGEHAMR_PRJ_FILL_ADD_RHS                                            \
    virtual void FillAddRhs(amrex::MultiFab &rhs_mf,                           \
                            const amrex::MultiFab &state_mf,                   \
                            const double time, const int lev, const double dt, \
                            const double dx, const double weight) override {   \
        performance_monitor->Start(performance_monitor->idx_rhs, lev);         \
        const int ncomp = rhs_mf.nComp();                                      \
        SLEDGEHAMR_KO_LOCAL_SETUP                                              \
        SLEDGEHAMR_RHS_PARAMS_LOCAL_SETUP                                      \
        SLEDGEHAMR_RHS_GW_PARAMS_LOCAL_SETUP                                   \
        SLEDGEHAMR_DO_PRAGMA(                                                  \
            omp parallel if (amrex::Gpu::notInLaunchRegion()))                 \
        for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU());               \
             mfi.isValid(); ++mfi) {                                           \
            const amrex::Box &bx = mfi.tilebox();                              \
            const amrex::Array4<double> &rhs_fab = rhs_mf.array(mfi);          \
            const amrex::Array4<double const> &state_fab =                     \
                state_mf.array(mfi);                                           \
            if (with_gravitational_waves) {                                    \
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,      \
                                                            int k) noexcept {  \
                    double tmp_rhs[Gw::NGwScalars];                            \
                    sledgehamr::utils::constexpr_for<0, Gw::NGwScalars, 1>(    \
                        [&](auto n) { tmp_rhs[n] = rhs_fab(i, j, k, n); });    \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs);                                         \
                    GravitationalWavesRhs<true>(rhs_fab, state_fab, i, j, k,   \
                                                lev, time, dt, dx,             \
                                                l_params_gw_rhs);              \
                    GravitationalWavesBackreaction<true>(                      \
                        rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs, l_params_gw_rhs);                        \
                    switch (l_with_dissipation * l_dissipation_order) {        \
                    case 2:                                                    \
                        sledgehamr::utils::constexpr_for<0, Gw::NGwScalars,    \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    2>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    case 3:                                                    \
                        sledgehamr::utils::constexpr_for<0, Gw::NGwScalars,    \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    3>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    }                                                          \
                    sledgehamr::utils::constexpr_for<0, Gw::NGwScalars, 1>(    \
                        [&](auto n) {                                          \
                            rhs_fab(i, j, k, n) += weight * tmp_rhs[n];        \
                        });                                                    \
                });                                                            \
            } else {                                                           \
                amrex::ParallelFor(bx, [=] AMREX_GPU_DEVICE(int i, int j,      \
                                                            int k) noexcept {  \
                    double tmp_rhs[Scalar::NScalars];                          \
                    sledgehamr::utils::constexpr_for<0, Scalar::NScalars, 1>(  \
                        [&](auto n) { tmp_rhs[n] = rhs_fab(i, j, k, n); });    \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx,        \
                        l_params_rhs);                                         \
                    switch (l_with_dissipation * l_dissipation_order) {        \
                    case 2:                                                    \
                        sledgehamr::utils::constexpr_for<0, Scalar::NScalars,  \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    2>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    case 3:                                                    \
                        sledgehamr::utils::constexpr_for<0, Scalar::NScalars,  \
                                                         1>([&](auto n) {      \
                            rhs_fab(i, j, k, n) +=                             \
                                sledgehamr::kernels::KreissOligerDissipation<  \
                                    3>(state_fab, i, j, k, n, dx,              \
                                       l_dissipation_strength[n]);             \
                        });                                                    \
                        break;                                                 \
                    }                                                          \
                    sledgehamr::utils::constexpr_for<0, Scalar::NScalars, 1>(  \
                        [&](auto n) {                                          \
                            rhs_fab(i, j, k, n) += weight * tmp_rhs[n];        \
                        });                                                    \
                });                                                            \
            }                                                                  \
        }                                                                      \
        performance_monitor->Stop(performance_monitor->idx_rhs, lev);          \
    };

/** @brief Overrides function in project class. Does tagging on CPUs.
 */
#define SLEDGEHAMR_PRJ_TAG_WITH_TRUNCATION_CPU                                 \
    virtual void TagWithTruncationCpu(                                         \
        const amrex::Array4<double const> &state_fab,                          \
        const amrex::Array4<double const> &state_fab_te,                       \
        const amrex::Array4<char> &tagarr, const amrex::Box &tilebox,          \
        double time, int lev, long *ntags_total, long *ntags_user,             \
        long *ntags_trunc, const std::vector<double> &params_tag,              \
        const std::vector<double> &params_mod) override {                      \
        const amrex::Dim3 lo = amrex::lbound(tilebox);                         \
        const amrex::Dim3 hi = amrex::ubound(tilebox);                         \
        if (with_gravitational_waves) {                                        \
            for (int k = lo.z; k <= hi.z; ++k) {                               \
                for (int j = lo.y; j <= hi.y; ++j) {                           \
                    AMREX_PRAGMA_SIMD                                          \
                    for (int i = lo.x; i <= hi.x; ++i) {                       \
                        tagarr(i, j, k) = amrex::TagBox::CLEAR;                \
                        bool res = false;                                      \
                        res = TagCellForRefinement<true>(                      \
                            state_fab, i, j, k, lev, time, dt[lev], dx[lev],   \
                            params_tag.data());                                \
                        if (res) {                                             \
                            tagarr(i, j, k) = amrex::TagBox::SET;              \
                            (*ntags_user)++;                                   \
                            (*ntags_total)++;                                  \
                        }                                                      \
                        bool te_res = TruncationErrorTagCpu<Gw::NGwScalars>(   \
                            state_fab, state_fab_te, i, j, k, lev, time,       \
                            dt[lev], dx[lev], te_crit, ntags_trunc,            \
                            params_mod);                                       \
                        if (te_res) {                                          \
                            tagarr(i, j, k) = amrex::TagBox::SET;              \
                            tagarr(i + 1, j, k) = amrex::TagBox::SET;          \
                            tagarr(i, j + 1, k) = amrex::TagBox::SET;          \
                            tagarr(i, j, k + 1) = amrex::TagBox::SET;          \
                            tagarr(i + 1, j + 1, k) = amrex::TagBox::SET;      \
                            tagarr(i, j + 1, k + 1) = amrex::TagBox::SET;      \
                            tagarr(i + 1, j, k + 1) = amrex::TagBox::SET;      \
                            tagarr(i + 1, j + 1, k + 1) = amrex::TagBox::SET;  \
                            (*ntags_total) += 8 - (int)res;                    \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        } else {                                                               \
            for (int k = lo.z; k <= hi.z; ++k) {                               \
                for (int j = lo.y; j <= hi.y; ++j) {                           \
                    AMREX_PRAGMA_SIMD                                          \
                    for (int i = lo.x; i <= hi.x; ++i) {                       \
                        tagarr(i, j, k) = amrex::TagBox::CLEAR;                \
                        bool res = false;                                      \
                        res = TagCellForRefinement<true>(                      \
                            state_fab, i, j, k, lev, time, dt[lev], dx[lev],   \
                            params_tag.data());                                \
                        if (res) {                                             \
                            tagarr(i, j, k) = amrex::TagBox::SET;              \
                            (*ntags_user)++;                                   \
                            (*ntags_total)++;                                  \
                        }                                                      \
                        bool te_res = TruncationErrorTagCpu<Scalar::NScalars>( \
                            state_fab, state_fab_te, i, j, k, lev, time,       \
                            dt[lev], dx[lev], te_crit, ntags_trunc,            \
                            params_mod);                                       \
                        if (te_res) {                                          \
                            tagarr(i, j, k) = amrex::TagBox::SET;              \
                            tagarr(i + 1, j, k) = amrex::TagBox::SET;          \
                            tagarr(i, j + 1, k) = amrex::TagBox::SET;          \
                            tagarr(i, j, k + 1) = amrex::TagBox::SET;          \
                            tagarr(i + 1, j + 1, k) = amrex::TagBox::SET;      \
                            tagarr(i, j + 1, k + 1) = amrex::TagBox::SET;      \
                            tagarr(i + 1, j, k + 1) = amrex::TagBox::SET;      \
                            tagarr(i + 1, j + 1, k + 1) = amrex::TagBox::SET;  \
                            (*ntags_total) += 8 - (int)res;                    \
                        }                                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    };

/** @brief Overrides function in project class. Does tagging on GPUs.
 */
#define SLEDGEHAMR_PRJ_TAG_WITH_TRUNCATION_GPU                                 \
    virtual void TagWithTruncationGpu(                                         \
        const amrex::Array4<double const> &state_fab,                          \
        const amrex::Array4<double const> &state_fab_te,                       \
        const amrex::Array4<char> &tagarr, const amrex::Box &tilebox,          \
        double time, int lev, const std::vector<double> &params_tag,           \
        const std::vector<double> &params_mod) override {                      \
        amrex::Gpu::AsyncArray<double> async_te_crit(te_crit.data(),           \
                                                     te_crit.size());          \
        double *l_te_crit = async_te_crit.data();                              \
        amrex::Gpu::AsyncArray<double> async_params_mod(params_mod.data(),     \
                                                        params_mod.size());    \
        double *l_params_mod = async_params_mod.data();                        \
        amrex::Gpu::AsyncArray<double> async_params_tag(params_tag.data(),     \
                                                        params_tag.size());    \
        double *l_params_tag = async_params_tag.data();                        \
        double l_dt = dt[lev];                                                 \
        double l_dx = dx[lev];                                                 \
        if (with_gravitational_waves) {                                        \
            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j,     \
                                                             int k) noexcept { \
                tagarr(i, j, k) = amrex::TagBox::CLEAR;                        \
                bool res = TagCellForRefinement<true>(                         \
                    state_fab, i, j, k, lev, time, l_dt, l_dx, l_params_tag);  \
                if (res) {                                                     \
                    tagarr(i, j, k) = amrex::TagBox::SET;                      \
                }                                                              \
                bool te_res = TruncationErrorTagGpu<Gw::NGwScalars>(           \
                    state_fab, state_fab_te, i, j, k, lev, time, l_dt, l_dx,   \
                    l_te_crit, l_params_mod);                                  \
                if (te_res) {                                                  \
                    tagarr(i, j, k) = amrex::TagBox::SET;                      \
                    tagarr(i + 1, j, k) = amrex::TagBox::SET;                  \
                    tagarr(i, j + 1, k) = amrex::TagBox::SET;                  \
                    tagarr(i, j, k + 1) = amrex::TagBox::SET;                  \
                    tagarr(i + 1, j + 1, k) = amrex::TagBox::SET;              \
                    tagarr(i, j + 1, k + 1) = amrex::TagBox::SET;              \
                    tagarr(i + 1, j, k + 1) = amrex::TagBox::SET;              \
                    tagarr(i + 1, j + 1, k + 1) = amrex::TagBox::SET;          \
                }                                                              \
            });                                                                \
        } else {                                                               \
            amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j,     \
                                                             int k) noexcept { \
                tagarr(i, j, k) = amrex::TagBox::CLEAR;                        \
                bool res = TagCellForRefinement<true>(                         \
                    state_fab, i, j, k, lev, time, l_dt, l_dx, l_params_tag);  \
                if (res) {                                                     \
                    tagarr(i, j, k) = amrex::TagBox::SET;                      \
                }                                                              \
                bool te_res = TruncationErrorTagGpu<Scalar::NScalars>(         \
                    state_fab, state_fab_te, i, j, k, lev, time, l_dt, l_dx,   \
                    l_te_crit, l_params_mod);                                  \
                if (te_res) {                                                  \
                    tagarr(i, j, k) = amrex::TagBox::SET;                      \
                    tagarr(i + 1, j, k) = amrex::TagBox::SET;                  \
                    tagarr(i, j + 1, k) = amrex::TagBox::SET;                  \
                    tagarr(i, j, k + 1) = amrex::TagBox::SET;                  \
                    tagarr(i + 1, j + 1, k) = amrex::TagBox::SET;              \
                    tagarr(i, j + 1, k + 1) = amrex::TagBox::SET;              \
                    tagarr(i + 1, j, k + 1) = amrex::TagBox::SET;              \
                    tagarr(i + 1, j + 1, k + 1) = amrex::TagBox::SET;          \
                }                                                              \
            });                                                                \
        }                                                                      \
    };

/** @brief Overrides function in project class. Does tagging on CPUs but without
 *         using truncation error estimates.
 */
#define SLEDGEHAMR_PRJ_TAG_WITHOUT_TRUNCATION_CPU                              \
    virtual void TagWithoutTruncationCpu(                                      \
        const amrex::Array4<double const> &state_fab,                          \
        const amrex::Array4<char> &tagarr, const amrex::Box &tilebox,          \
        double time, int lev, long *ntags_total,                               \
        const std::vector<double> &params) override {                          \
        const amrex::Dim3 lo = amrex::lbound(tilebox);                         \
        const amrex::Dim3 hi = amrex::ubound(tilebox);                         \
        for (int k = lo.z; k <= hi.z; ++k) {                                   \
            for (int j = lo.y; j <= hi.y; ++j) {                               \
                AMREX_PRAGMA_SIMD                                              \
                for (int i = lo.x; i <= hi.x; ++i) {                           \
                    tagarr(i, j, k) = amrex::TagBox::CLEAR;                    \
                    bool res = TagCellForRefinement<true>(                     \
                        state_fab, i, j, k, lev, time, dt[lev], dx[lev],       \
                        params.data());                                        \
                    if (res) {                                                 \
                        tagarr(i, j, k) = amrex::TagBox::SET;                  \
                        (*ntags_total)++;                                      \
                    }                                                          \
                }                                                              \
            }                                                                  \
        }                                                                      \
    };

/** @brief Overrides function in project class. Does tagging on GPUs but without
 *         using truncation error estimates.
 */
#define SLEDGEHAMR_PRJ_TAG_WITHOUT_TRUNCATION_GPU                              \
    virtual void TagWithoutTruncationGpu(                                      \
        const amrex::Array4<double const> &state_fab,                          \
        const amrex::Array4<char> &tagarr, const amrex::Box &tilebox,          \
        double time, int lev, const std::vector<double> &params) override {    \
        amrex::Gpu::AsyncArray async_params(params.data(), params.size());     \
        double *l_params = async_params.data();                                \
        double l_dt = dt[lev];                                                 \
        double l_dx = dx[lev];                                                 \
        amrex::ParallelFor(tilebox, [=] AMREX_GPU_DEVICE(int i, int j,         \
                                                         int k) noexcept {     \
            tagarr(i, j, k) = amrex::TagBox::CLEAR;                            \
            bool res = TagCellForRefinement<true>(state_fab, i, j, k, lev,     \
                                                  time, l_dt, l_dx, l_params); \
            if (res) {                                                         \
                tagarr(i, j, k) = amrex::TagBox::SET;                          \
            }                                                                  \
        });                                                                    \
    };

/** @brief Add boilerplate functions to project class.
 */
#define SLEDGEHAMR_INITIALIZE_PROJECT(prj)                                     \
    SLEDGEHAMR_PRJ_CONSTRUCTOR(prj)                                            \
    SLEDGEHAMR_PRJ_FILL_RHS                                                    \
    SLEDGEHAMR_PRJ_FILL_ADD_RHS                                                \
    SLEDGEHAMR_PRJ_TAG_WITH_TRUNCATION_CPU                                     \
    SLEDGEHAMR_PRJ_TAG_WITH_TRUNCATION_GPU                                     \
    SLEDGEHAMR_PRJ_TAG_WITHOUT_TRUNCATION_CPU                                  \
    SLEDGEHAMR_PRJ_TAG_WITHOUT_TRUNCATION_GPU

}; // namespace sledgehamr

#endif // SLEDGEHAMR_MACROS_H_
