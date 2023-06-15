#ifndef SLEDGEHAMR_KERNELS_H_
#define SLEDGEHAMR_KERNELS_H_

namespace sledgehamr {
namespace kernels {

/** @brief Averages down a set of fine level cells onto a coarse level cell.
 *         Computes truncation errors at the same time.
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
 * @param   ncomp   Number of scalar fields.
 * @param   crse    Coarse level data.
 * @param   fine    Fine level data.
 * @param   te      Container to save truncation errors. Same dimensions as
 *                  fine.
 */
AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
void AverageDownWithTruncationError(int i, int j, int k, const int ncomp,
                                    amrex::Array4<double> const& crse,
                                    amrex::Array4<double const> const& fine,
                                    amrex::Array4<double> const& te) {
    const int ratio = 2;
    const double volfrac = 1.0/(double)(ratio*ratio*ratio);
    const int ii = i*ratio;
    const int jj = j*ratio;
    const int kk = k*ratio;

    for (int n = 0; n < ncomp; ++n ) {
        double c = fine(ii, jj, kk, n);
        c += fine(ii+1, jj  , kk  , n);
        c += fine(ii  , jj+1, kk  , n);
        c += fine(ii  , jj  , kk+1, n);
        c += fine(ii+1, jj+1, kk  , n);
        c += fine(ii+1, jj  , kk+1, n);
        c += fine(ii  , jj+1, kk+1, n);
        c += fine(ii+1, jj+1, kk+1, n);

        te(ii,jj,kk,n) = fabs(crse(i,j,k,n) - volfrac * c);
        crse(i,j,k,n) = volfrac * c;
    }
}

template<int> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double KreisOligerDissipation(const amrex::Array4<double const>& state,
                              const int i, const int j, const int k,
                              const int c, const double dx,
                              const double dissipation_strength);

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double KreisOligerDissipation<2>(const amrex::Array4<double const>& state,
                                 const int i, const int j, const int k,
                                 const int c, const double dx,
                                 const double dissipation_strength) {
    double dx4 = state(i+2,j,k,c) - 4.*state(i+1,j,k,c)
               + state(i-2,j,k,c) - 4.*state(i-1,j,k,c) 
               + 6.*state(i,j,k,c);

    double dy4 = state(i,j+2,k,c) - 4.*state(i,j+1,k,c)
               + state(i,j-2,k,c) - 4.*state(i,j-1,k,c) 
               + 6.*state(i,j,k,c);

    double dz4 = state(i,j,k+2,c) - 4.*state(i,j,k+1,c)
               + state(i,j,k-2,c) - 4.*state(i,j,k-1,c)
               + 6.*state(i,j,k,c);

    return - dissipation_strength * (dx4+dy4+dz4) / 16. / dx;
}

template<> AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
double KreisOligerDissipation<3>(const amrex::Array4<double const>& state,
                                 const int i, const int j, const int k,
                                 const int c, const double dx,
                                 const double dissipation_strength) {
    double dx6 = state(i+3,j,k,c) - 6.*state(i+2,j,k,c) + 15.*state(i+1,j,k,c)
               + state(i-3,j,k,c) - 6.*state(i-2,j,k,c) + 15.*state(i-1,j,k,c)
               - 20.*state(i,j,k,c);

    double dy6 = state(i,j+3,k,c) - 6.*state(i,j+2,k,c) + 15.*state(i,j+1,k,c)
               + state(i,j-3,k,c) - 6.*state(i,j-2,k,c) + 15.*state(i,j-1,k,c)
               - 20.*state(i,j,k,c);

    double dz6 = state(i,j,k+3,c) - 6.*state(i,j,k+2,c) + 15.*state(i,j,k+1,c)
               + state(i,j,k-3,c) - 6.*state(i,j,k-2,c) + 15.*state(i,j,k-1,c)
               - 20.*state(i,j,k,c);

    return dissipation_strength * (dx6+dy6+dz6) / 64. / dx;
}

}; // namespace kernels
}; // namespace sledgehamr

#endif // SLEDGEHAMR_KERNELS_H_
