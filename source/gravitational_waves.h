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
    double IndexToK(int a, int N);
    double GetProjection(int i, int j, int abc[3], int N);
    double GetLambda(int i, int j, int l, int m, int abc[3], int N);

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
                                           const double dk, const int dimN){};
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
