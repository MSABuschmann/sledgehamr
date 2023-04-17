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

    // TODO: Do custom functions

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

/** @brief Identifies cells that violate the truncation error threshold.
 * @param   i           i-th cell index.
 * @param   j           j-th cell index.
 * @param   k           k-th cell index.
 * @param   time        Current time.
 * @param   lev         Current level.
 * @param   te          Truncation errors.
 * @param   ntags_trunc Counts number of tags.
 * @return  If truncation error threshold has been exceeded.
 */
/*AMREX_FORCE_INLINE
bool TruncationErrorTagCpu(int i, int j, int k, double time, int lev,
                           const amrex::Array4<double const>& te,
                           std::vector<double>& te_crit, int* ntags_trunc) {
    // truncation errors for scalar field components saved in even indices.
    if (i%2 != 0 || j%2 != 0 || k%2 != 0)
        return false;

    bool res = false;
    for (int n=0; n<te.nComp(); ++n) {
        if (te(i,j,k,n) >= te_crit[n] ){
            res = true;
            ntags_trunc[n] += 8;
        }
    }

    return res;
};*/

/* @brief Identifies cells that violate the truncation error threshold.
 * @param   i       i-th cell index.
 * @param   j       j-th cell index.
 * @param   k       k-th cell index.
 * @param   time    Current time.
 * @param   lev     Current level.
 * @param   te      Truncation errors.
 * @return  If truncation error threshold has been exceeded.
 */
/*AMREX_GPU_HOST_DEVICE AMREX_FORCE_INLINE
bool TruncationErrorTagGpu(int i, int j, int k, double time, int lev,
                           const amrex::Array4<double const>& te,
                           double* te_crit) {
    // truncation errors for scalar field components saved in even indices.
    if (i%2 != 0 || j%2 != 0 || k%2 != 0)
        return false;

    bool res = false;
    for (int n=0; n<te.nComp(); ++n) {
        if (te(i,j,k,n) >= te_crit[n]) {
            res = true;
        }
    }

    return res;
};*/

}; // namespace kernels
}; // namespace sledgehamr

#endif // SLEDGEHAMR_KERNELS_H_
