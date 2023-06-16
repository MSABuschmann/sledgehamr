#ifndef SLEDGEHAMR_MACROS_H_
#define SLEDGEHAMR_MACROS_H_

#include "kernels.h"

/** Yes, I know this file is disgusting. This is what we get for not 
 *  wanting everyone to write a bunch of complicated boilerplate every single
 *  time but also want to avoid all this function overhead for performance ...
 */

namespace sledgehamr {

// Force boost to define variadics.
#define BOOST_PP_VARIADICS 1
#include <boost/preprocessor.hpp>

/** @brief Element-wise vector reduction for OpenMP.
 */
#pragma omp declare reduction(vec_int_plus : std::vector<int> : \
        std::transform(omp_out.begin(), omp_out.end(), omp_in.begin(), \
                       omp_out.begin(), std::plus<int>())) \
        initializer(omp_priv = omp_orig)

#define DO_PRAGMA(x) _Pragma(#x)

/** @brief Macros to declare and initialize scalar fields within project
 *         namespace.
 */
#define ADD_SCALAR(field) \
    static sledgehamr::ScalarField BOOST_PP_CAT(_s_, field)  = \
            {#field, l_scalar_fields, false}; \
    namespace Scalar { \
        constexpr int field = _Scalar::field; \
    };

#define EXPAND_SCALARS(r, data, field) ADD_SCALAR(field)

#define ADD_MOMENTUM(field) \
    static sledgehamr::ScalarField BOOST_PP_CAT(_s_, field)  = \
            {#field, l_scalar_fields, true}; \
    namespace Scalar { \
        constexpr int field = _Momentum::field + _Scalar::NScalarFields; \
    };

#define EXPAND_MOMENTA(r, data, field) ADD_MOMENTUM(field)

/** @brief Macros to create enum of fields for fast and convinient component
 *         access within the project namespace.
 */
#define SCALAR_ENUM_VALUE(r, data, elem) elem,

#define SCALAR_ENUM(...) \
    enum _Scalar { BOOST_PP_SEQ_FOR_EACH(SCALAR_ENUM_VALUE, _, \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) NScalarFields \
    };

#define MOMENTUM_ENUM(...) \
    enum _Momentum { BOOST_PP_SEQ_FOR_EACH(SCALAR_ENUM_VALUE, _, \
            BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) NMomentumFields \
    };

#define GW_ENUM enum Gw { u_xx = Scalar::NScalars, u_yy, u_zz, u_xy, u_xz, \
                          u_yz, du_xx, du_yy, du_zz, du_xy, du_xz, du_yz, \
                          NGwScalars \
    };

/* brief Default implementation for f(\tau) > \tau_{crit} criteria:
 *       f(\tau) = \tau. This template can be specialized for each scalar field
 *       component by the project.
 */
#define TRUNCATION_MODIFIER template<int> \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
    double TruncationModifier(const amrex::Array4<const double>& state, \
                              const int i, const int j, const int k, \
                              const int lev, const double time, \
                              const double dt, const double dx, \
                              const double truncation_error) { \
        return truncation_error; \
    };

/* brief User-defined tagging criteria, switched off by default. Template <true>
 *       will be used by tagger and this specialisation can be implemented
 *       within the project.
 */
#define TAG_CELL_FOR_REFINEMENT template<bool> \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
    bool TagCellForRefinement(const amrex::Array4<const double>& state, \
                              const int i, const int j, const int k, \
                              const int lev, const double time, \
                              const double dt, const double dx) { \
        return false; \
    };

/* brief TODO
 */
#define GRAVITATIONAL_WAVES_RHS template<bool> \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
    void GravitationalWavesRhs(const amrex::Array4<double>& rhs, \
            const amrex::Array4<const double>& state, \
            const int i, const int j, const int k, const int lev, \
            const double time, const double dt, const double dx) { \
    };

/** @brief Macro to add multiple scalar fields and default template functions to
 *         the project namespace.
 */
#define ADD_SCALARS(...) \
    SCALAR_ENUM(__VA_ARGS__) \
    static std::vector<sledgehamr::ScalarField*> l_scalar_fields; \
    BOOST_PP_SEQ_FOR_EACH(EXPAND_SCALARS, _, \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__))

#define ADD_CONJUGATE_MOMENTA(...) \
    MOMENTUM_ENUM(__VA_ARGS__) \
    BOOST_PP_SEQ_FOR_EACH(EXPAND_MOMENTA, _, \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    namespace Scalar { \
        constexpr int NScalars = _Scalar::NScalarFields \
                               + _Momentum::NMomentumFields; \
    }; \
    GW_ENUM \
    TRUNCATION_MODIFIER \
    TAG_CELL_FOR_REFINEMENT \
    GRAVITATIONAL_WAVES_RHS


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
 * @return  If truncation error threshold has been exceeded.
 */
#define TRUNCATION_ERROR_TAG_CPU template<int NScalars> \
    AMREX_FORCE_INLINE \
    bool TruncationErrorTagCpu(const amrex::Array4<const double>& state, \
                               const amrex::Array4<const double>& te, \
                               const int i, const int j, const int k, \
                               const int lev, const double time, \
                               const double dt, const double dx, \
                               std::vector<double>& te_crit, \
                               int* ntags_trunc) { \
        if (i%2 != 0 || j%2 != 0 || k%2 != 0) \
            return false; \
        bool res = false; \
        sledgehamr::utils::constexpr_for<0, NScalars, 1>([&](auto n) { \
            double mte = TruncationModifier<n>(state, i, j, k, lev, time, dt, \
                                               dx, te(i,j,k,n)); \
            if (mte >= te_crit[n]) { \
                res = true; \
                ntags_trunc[n] += 8; \
            } \
        }); \
        return res; \
    };

/** @brief Same as TRUNCATION_ERROR_TAG_CPU but to be used for GPU code.
 */
#define TRUNCATION_ERROR_TAG_GPU template<int NScalars> \
    AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE \
    bool TruncationErrorTagGpu(const amrex::Array4<double const>& state, \
                               const amrex::Array4<double const>& te, \
                               const int i, const int j, const int k, \
                               const int lev, const double time, \
                               const double dt, const double dx, \
                               double* te_crit) { \
        if (i%2 != 0 || j%2 != 0 || k%2 != 0) \
            return false; \
        bool res = false; \
        sledgehamr::utils::constexpr_for<0, NScalars, 1>([&](auto n) { \
            double mte = TruncationModifier<n>(state, i, j, k, lev, time, dt, \
                                               dx, te(i,j,k,n)); \
            if (mte >= te_crit[n]) { \
                res = true; \
            } \
        }); \
        return res; \
    };

/** @brief Add functions to project namespace that depend on custom template
 *         specialisations.
 */
#define FINISH_SLEDGEHAMR_SETUP TRUNCATION_ERROR_TAG_CPU \
    TRUNCATION_ERROR_TAG_GPU

/** @brief Constructor of project class to initialize scalar fields
 *         automatically.
 */
#define PRJ_CONSTRUCTOR(prj) \
    prj (){ \
        scalar_fields = l_scalar_fields; \
        amrex::Print() << "Starting "  << #prj << " project..." << std::endl; \
        amrex::Print() << "Number of field components: " \
                       << scalar_fields.size() << std::endl; \
        amrex::Print() << std::endl; \
    };

/** @brief Computes Rhs for the entire grid.
 * @param   rhs_mf      Container to fill the Rhs with.
 * @param   state_mf    Current grid.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   dt          Time step size.
 * @param   dx          Grid spacing.
 */
#define PRJ_FILL_RHS virtual void FillRhs(amrex::MultiFab& rhs_mf, \
                const amrex::MultiFab& state_mf, const double time, \
                const int lev, const double dt, const double dx) override { \
        double* l_dissipation_strength = dissipation_strength.data(); \
        const int l_dissipation_order = dissipation_order; \
        DO_PRAGMA(omp parallel if (amrex::Gpu::notInLaunchRegion())) \
        for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); \
             mfi.isValid(); ++mfi) { \
            const amrex::Box& bx = mfi.tilebox(); \
            const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi); \
            const amrex::Array4<double const>& state_fab = \
                    state_mf.array(mfi); \
            if (with_gravitational_waves) { \
                amrex::ParallelFor(bx, \
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
                { \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx); \
                    GravitationalWavesRhs<true>(rhs_fab, state_fab, i, j, k, \
                                                lev, time, dt, dx); \
                    switch (with_dissipation * dissipation_order) { \
                        case 2: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Gw::NGwScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<2>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                        case 3: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Gw::NGwScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<3>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                    } \
                }); \
            } else { \
                amrex::ParallelFor(bx, \
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
                { \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx); \
                    switch (with_dissipation * dissipation_order) { \
                        case 2: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Scalar::NScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<2>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                        case 3: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Scalar::NScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<3>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                    } \
               }); \
            } \
        } \
    };

/** @brief
 */
#define PRJ_FILL_ADD_RHS virtual void FillAddRhs(amrex::MultiFab& rhs_mf, \
                const amrex::MultiFab& state_mf, const double time, \
                const int lev, const double dt, const double dx, \
                const double weight) override { \
        const int ncomp = rhs_mf.nComp(); \
        double* l_dissipation_strength = dissipation_strength.data(); \
        const int l_dissipation_order = dissipation_order; \
        DO_PRAGMA(omp parallel if (amrex::Gpu::notInLaunchRegion())) \
        for (amrex::MFIter mfi(rhs_mf, amrex::TilingIfNotGPU()); \
             mfi.isValid(); ++mfi) { \
            const amrex::Box& bx = mfi.tilebox(); \
            const amrex::Array4<double>& rhs_fab = rhs_mf.array(mfi); \
            const amrex::Array4<double const>& state_fab = \
                    state_mf.array(mfi); \
            if (with_gravitational_waves) { \
                amrex::ParallelFor(bx, \
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
                { \
                    double tmp_rhs[Gw::NGwScalars]; \
                    sledgehamr::utils::constexpr_for<0, Gw::NGwScalars, 1> \
                            ([&](auto n) { \
                        tmp_rhs[n] = rhs_fab(i, j, k, n); \
                    }); \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx); \
                    GravitationalWavesRhs<true>(rhs_fab, state_fab, i, j, k, \
                                                lev, time, dt, dx); \
                    switch (with_dissipation * dissipation_order) { \
                        case 2: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Gw::NGwScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<2>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                        case 3: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Gw::NGwScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<3>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                    } \
                    sledgehamr::utils::constexpr_for<0, Gw::NGwScalars, 1> \
                            ([&](auto n) { \
                        rhs_fab(i, j, k, n) += weight * tmp_rhs[n]; \
                    }); \
                }); \
            } else { \
                amrex::ParallelFor(bx, \
                [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept \
                { \
                    double tmp_rhs[Scalar::NScalars]; \
                    sledgehamr::utils::constexpr_for<0, Scalar::NScalars, 1> \
                            ([&](auto n) { \
                        tmp_rhs[n] = rhs_fab(i, j, k, n); \
                    }); \
                    Rhs(rhs_fab, state_fab, i, j, k, lev, time, dt, dx); \
                    switch (with_dissipation * dissipation_order) { \
                        case 2: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Scalar::NScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<2>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                        case 3: \
                            sledgehamr::utils::constexpr_for \
                                    <0, Scalar::NScalars, 1> ([&](auto n) { \
                                rhs_fab(i,j,k,n) += \
                            sledgehamr::kernels::KreissOligerDissipation<3>( \
                                        state_fab, i, j, j, n, dx, \
                                        l_dissipation_strength[n]); \
                            }); \
                            break; \
                    } \
                    sledgehamr::utils::constexpr_for<0, Scalar::NScalars, 1> \
                            ([&](auto n) { \
                        rhs_fab(i, j, k, n) += weight * tmp_rhs[n]; \
                    }); \
                }); \
            } \
        } \
    };

/** @brief Overrides function in project class.
 */
#define PRJ_TAG_WITH_TRUNCATION_CPU virtual void TagWithTruncationCpu( \
            const amrex::Array4<double const>& state_fab, \
            const amrex::Array4<double const>& state_fab_te, \
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox, \
            double time, int lev, int* ntags_total, int* ntags_user, \
            int* ntags_trunc) override { \
        const amrex::Dim3 lo = amrex::lbound(tilebox); \
        const amrex::Dim3 hi = amrex::ubound(tilebox); \
        if (with_gravitational_waves) { \
            for (int k = lo.z; k <= hi.z; ++k) { \
                for (int j = lo.y; j <= hi.y; ++j) { \
                    AMREX_PRAGMA_SIMD \
                    for (int i = lo.x; i <= hi.x; ++i) { \
                        tagarr(i,j,k) = amrex::TagBox::CLEAR; \
                        bool res = false; \
                        res = TagCellForRefinement<true>( \
                                state_fab, i, j, k, lev, time, dt[lev], \
                                dx[lev]); \
                        if( res ){ \
                            tagarr(i,j,k) = amrex::TagBox::SET; \
                            (*ntags_user)++; \
                            (*ntags_total)++; \
                        } \
                        bool te_res = TruncationErrorTagCpu<Gw::NGwScalars>( \
                                state_fab, state_fab_te, i, j, k, lev, time, \
                                dt[lev], dx[lev], te_crit, ntags_trunc); \
                        if (te_res) { \
                            tagarr(i  ,j  ,k  ) = amrex::TagBox::SET; \
                            tagarr(i+1,j  ,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j+1,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j  ,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j+1,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j+1,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j  ,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j+1,k+1) = amrex::TagBox::SET; \
                            (*ntags_total) += 8 - (int)res; \
                        } \
                    } \
                } \
            } \
        } else { \
            for (int k = lo.z; k <= hi.z; ++k) { \
                for (int j = lo.y; j <= hi.y; ++j) { \
                    AMREX_PRAGMA_SIMD \
                    for (int i = lo.x; i <= hi.x; ++i) { \
                        tagarr(i,j,k) = amrex::TagBox::CLEAR; \
                        bool res = false; \
                        res = TagCellForRefinement<true>( \
                                state_fab, i, j, k, lev, time, dt[lev], \
                                dx[lev]); \
                        if( res ){ \
                            tagarr(i,j,k) = amrex::TagBox::SET; \
                            (*ntags_user)++; \
                            (*ntags_total)++; \
                        } \
                        bool te_res = TruncationErrorTagCpu<Scalar::NScalars>( \
                                state_fab, state_fab_te, i, j, k, lev, time, \
                                dt[lev], dx[lev], te_crit, ntags_trunc); \
                        if (te_res) { \
                            tagarr(i  ,j  ,k  ) = amrex::TagBox::SET; \
                            tagarr(i+1,j  ,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j+1,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j  ,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j+1,k  ) = amrex::TagBox::SET; \
                            tagarr(i  ,j+1,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j  ,k+1) = amrex::TagBox::SET; \
                            tagarr(i+1,j+1,k+1) = amrex::TagBox::SET; \
                            (*ntags_total) += 8 - (int)res; \
                        } \
                    } \
                } \
            } \
        } \
    };

/** @brief Overrides function in project class.
 */
#define PRJ_TAG_WITH_TRUNCATION_GPU virtual void TagWithTruncationGpu( \
            const amrex::Array4<double const>& state_fab, \
            const amrex::Array4<double const>& state_fab_te, \
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox, \
            double time, int lev) override { \
        amrex::Gpu::AsyncArray<double> l_te_crit_arr(&te_crit[0], \
                                                     te_crit.size()); \
        double l_dt = dt[lev]; \
        double l_dx = dx[lev]; \
        double* l_te_crit = l_te_crit_arr.data(); \
        if (with_gravitational_waves) { \
            amrex::ParallelFor(tilebox, \
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { \
                tagarr(i,j,k) = amrex::TagBox::CLEAR; \
                bool res = TagCellForRefinement<true>(state_fab, i, j, k, lev, \
                                                      time, l_dt, l_dx); \
                if (res) { \
                    tagarr(i,j,k) = amrex::TagBox::SET; \
                } \
                bool te_res = TruncationErrorTagGpu<Gw::NGwScalars>( \
                        state_fab, state_fab_te, i, j, k, lev, time, l_dt, \
                        l_dx, l_te_crit); \
                if (te_res) { \
                    tagarr(i  ,j  ,k  ) = amrex::TagBox::SET; \
                    tagarr(i+1,j  ,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j+1,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j  ,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j+1,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j+1,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j  ,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j+1,k+1) = amrex::TagBox::SET; \
                } \
            }); \
        } else { \
            amrex::ParallelFor(tilebox, \
            [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { \
                tagarr(i,j,k) = amrex::TagBox::CLEAR; \
                bool res = TagCellForRefinement<true>(state_fab, i, j, k, lev, \
                                                      time, l_dt, l_dx); \
                if (res) { \
                    tagarr(i,j,k) = amrex::TagBox::SET; \
                } \
                bool te_res = TruncationErrorTagGpu<Scalar::NScalars>( \
                        state_fab, state_fab_te, i, j, k, lev, time, l_dt, \
                        l_dx, l_te_crit); \
                if (te_res) { \
                    tagarr(i  ,j  ,k  ) = amrex::TagBox::SET; \
                    tagarr(i+1,j  ,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j+1,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j  ,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j+1,k  ) = amrex::TagBox::SET; \
                    tagarr(i  ,j+1,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j  ,k+1) = amrex::TagBox::SET; \
                    tagarr(i+1,j+1,k+1) = amrex::TagBox::SET; \
                } \
            }); \
        } \
    };

/** @brief Overrides function in project class.
 */
#define PRJ_TAG_WITHOUT_TRUNCATION_CPU virtual void TagWithoutTruncationCpu( \
            const amrex::Array4<double const>& state_fab, \
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox, \
            double time, int lev, int* ntags_total) override { \
        const amrex::Dim3 lo = amrex::lbound(tilebox); \
        const amrex::Dim3 hi = amrex::ubound(tilebox); \
        for (int k = lo.z; k <= hi.z; ++k) { \
            for (int j = lo.y; j <= hi.y; ++j) { \
                AMREX_PRAGMA_SIMD \
                for (int i = lo.x; i <= hi.x; ++i) { \
                    tagarr(i,j,k) = amrex::TagBox::CLEAR; \
                    bool res = TagCellForRefinement<true>(state_fab, i, j, k, \
                            lev, time, dt[lev], dx[lev]); \
                    if (res) { \
                        tagarr(i,j,k) = amrex::TagBox::SET; \
                        (*ntags_total)++; \
                    } \
                } \
            } \
        } \
    };

/** @brief Overrides function in project class.
 */
#define PRJ_TAG_WITHOUT_TRUNCATION_GPU virtual void TagWithoutTruncationGpu( \
            const amrex::Array4<double const>& state_fab, \
            const amrex::Array4<char>& tagarr, const amrex::Box& tilebox, \
            double time, int lev) override { \
        double l_dt = dt[lev]; \
        double l_dx = dx[lev]; \
        amrex::ParallelFor(tilebox, \
        [=] AMREX_GPU_DEVICE (int i, int j, int k) noexcept { \
            tagarr(i,j,k) = amrex::TagBox::CLEAR; \
            bool res = TagCellForRefinement<true>(state_fab, i, j, k, lev, \
                                                  time, l_dt, l_dx); \
            if (res) { \
                tagarr(i,j,k) = amrex::TagBox::SET; \
            } \
        }); \
    };

/** @brief Add boilerplate functions to project class.
 */
#define START_PROJECT(prj) \
    PRJ_CONSTRUCTOR(prj) \
    PRJ_FILL_RHS \
    PRJ_FILL_ADD_RHS \
    PRJ_TAG_WITH_TRUNCATION_CPU \
    PRJ_TAG_WITH_TRUNCATION_GPU \
    PRJ_TAG_WITHOUT_TRUNCATION_CPU \
    PRJ_TAG_WITHOUT_TRUNCATION_GPU

}; // namespace sledgehamr

#endif // SLEDGEHAMR_MACROS_H_
