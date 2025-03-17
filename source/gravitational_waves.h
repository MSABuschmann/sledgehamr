#ifndef SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
#define SLEDGEHAMR_GRAVITATIONAL_WAVES_H_

#include <hdf5.h>

#include "sledgehamr.h"

namespace sledgehamr {

struct GravitationalWavesSpectrumModifier;
class Sledgehamr;

/** @brief This will setup everything needed to simulate gravitational waves and
 *         compute the gravitational wave spectrum.
 */
class GravitationalWaves {
  public:
    GravitationalWaves(Sledgehamr *owner);

    void
    ComputeSpectrum(hid_t file_id,
                    GravitationalWavesSpectrumModifier *modifier = nullptr);

    int GetProjectionType() const { return projection_type; }
    int GetZeroPadding() const { return zero_padding; }

    /** @brief enum with the tensor components.
     */
    enum Gw {
        u_xx = 0,
        u_yy,
        u_zz,
        u_xy,
        u_xz,
        u_yz,
        du_xx,
        du_yy,
        du_zz,
        du_xy,
        du_xz,
        du_yz,
        NGwScalars
    };

  private:
    /** @brief Pointer to the simulation.
     */
    Sledgehamr *sim;

    /** @brief Offset corresponding the number of user-defined scalar fields.
     *         We need this as all gravitional wave fields are appended to the
     *         end of that list.
     */
    int idx_offset;

    /** @brief Number of gravitational wave components
     */
    static constexpr int NScalars = 12;

    /** @brief Projection type when computing the spectrum. This should be set
     *         to the same order as the gradient computation in the EOM.
     *         However, since this defined in a kernel function we need to keep
     *         this is a free parameter here and rely on the user to select
     *         the appropiate projection type through the input file.
     *         Will be set by 'output.gw_spectra.projection_type'.
     */
    int projection_type = 2;

    /** @brief Defines the amount of zero padding during the FFT when computing
     *         the spectrum. Must be a power of 2. 1 means no zero padding.
     */
    int zero_padding = 1;

    bool unbinned = true;

    std::unique_ptr<GravitationalWavesSpectrumModifier> default_modifier;
};

struct GravitationalWavesSpectrumModifier {
    virtual void SelectComponents(int components[6]) {
        components[0] = GravitationalWaves::Gw::du_xx;
        components[1] = GravitationalWaves::Gw::du_xy;
        components[2] = GravitationalWaves::Gw::du_xz;
        components[3] = GravitationalWaves::Gw::du_yy;
        components[4] = GravitationalWaves::Gw::du_yz;
        components[5] = GravitationalWaves::Gw::du_zz;
    };

    virtual void FourierSpaceModifications(amrex::MultiFab du_real[6],
                                           amrex::MultiFab du_imag[6],
                                           const double dk, const int dimN) {};
};

/** @brief Projects all indicies.
 * @param   i   i-index.
 * @param   j   j-index.
 * @param   abc a-, b-, and c- index.
 * @param   N   Total index length.
 * @return Projected value.
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE double gw_GetProjection(int i, int j,
                                                            double abc[3]) {
    double norm = abc[0] * abc[0] + abc[1] * abc[1] + abc[2] * abc[2];
    int l = amrex::min(i, j);
    int m = amrex::max(i, j);

    double proj = abc[l] * abc[m];

    return static_cast<double>(l == m) - proj / norm;
}

/** @brief Computes lambda from projections.
 * @param   i   i-index.
 * @param   j   j-index.
 * @param   l   l-index.
 * @param   m   m-index.
 * @param   abc a-, b-, and c-index.
 * @param   N   Total index length.
 * @return Lambda
 */
AMREX_GPU_DEVICE AMREX_FORCE_INLINE double
gw_GetLambda(int i, int j, int l, int m, int abc[3], const double *index_to_k) {
    if (abc[0] == 0 && abc[1] == 0 && abc[2] == 0)
        return 0.;

    double abc_d[3];
    abc_d[0] = index_to_k[abc[0]];
    abc_d[1] = index_to_k[abc[1]];
    abc_d[2] = index_to_k[abc[2]];

    return gw_GetProjection(i, l, abc_d) * gw_GetProjection(j, m, abc_d) -
           gw_GetProjection(i, j, abc_d) * gw_GetProjection(l, m, abc_d) / 2.;
}

}; // namespace sledgehamr

#endif // SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
