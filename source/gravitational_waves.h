#ifndef SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
#define SLEDGEHAMR_GRAVITATIONAL_WAVES_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief This will setup everything needed to simulate gravitational waves and
 *         compute the gravitational wave spectrum.
 */
class GravitationalWaves {
  public:
    GravitationalWaves(Sledgehamr* owner);
    void ComputeSpectrum(hid_t file_id);

  private:
    double IndexToK(int a, int N);
    double GetProjection(int i, int j, int abc[3], int N);
    double GetLambda(int i, int j, int l, int m, int abc[3], int N);

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;

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

    /** @brief enum with the tensor components.
     */
    enum Gw {
        u_xx = 0, u_yy, u_zz, u_xy, u_xz, u_yz, du_xx, du_yy, du_zz, du_xy,
        du_xz, du_yz, NGwScalars
    };
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_GRAVITATIONAL_WAVES_H_
