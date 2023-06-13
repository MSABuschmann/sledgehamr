#include "gravitational_waves.h"

namespace sledgehamr {

GravitationalWaves::GravitationalWaves(Sledgehamr* owner) {
    sim = owner;

    ScalarField u_xx("u_xx", sim->scalar_fields);
    ScalarField u_yy("u_yy", sim->scalar_fields);
    ScalarField u_zz("u_zz", sim->scalar_fields);
    ScalarField u_xy("u_xy", sim->scalar_fields);
    ScalarField u_xz("u_xz", sim->scalar_fields);
    ScalarField u_yz("u_yz", sim->scalar_fields);
    ScalarField du_xx("du_xx", sim->scalar_fields);
    ScalarField du_yy("du_yy", sim->scalar_fields);
    ScalarField du_zz("du_zz", sim->scalar_fields);
    ScalarField du_xy("du_xy", sim->scalar_fields);
    ScalarField du_xz("du_xz", sim->scalar_fields);
    ScalarField du_yz("du_yz", sim->scalar_fields);
}

void GravitationalWaves::ComputeSpectrum(hid_t file_id) {

}

}; // namespace sledgehamr
