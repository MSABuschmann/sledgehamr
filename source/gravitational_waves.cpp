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

inline double GravitationalWaves::IndexToK(int a, int N) {
    double n_tilde = a-N <= -N/2-1 ? a : a-N;
    double two_pi_n_tilde = 2.*M_PI/static_cast<double>(N)*n_tilde;

    if (sim->nghost == 1) {
        return sin(two_pi_n_tilde);
    } else if (sim->nghost > 1) {
        return (8.*sin(two_pi_n_tilde) - sin(2.*two_pi_n_tilde)) / 6.;
    }

    return 0.;
}

inline double GravitationalWaves::GetProjection(int i, int j, int abc[3],
                                                int N) {
    if( abc[0] == 0 && abc[1] == 0 && abc[2] == 0)
        return 0.;

    double abc_d[3];
    abc_d[0] = IndexToK(abc[0], N);
    abc_d[1] = IndexToK(abc[1], N);
    abc_d[2] = IndexToK(abc[2], N);

    double norm = abc_d[0]*abc_d[0] + abc_d[1]*abc_d[1] + abc_d[2]*abc_d[2];
    int l = amrex::min(i, j);
    int m = amrex::max(i, j);

    double proj = abc_d[l] * abc_d[m];

    return static_cast<double>(l==m) - proj/norm;
}

inline double GravitationalWaves::GetLambda(int i, int j, int l, int m,
                                            int abc[3], int N) {
    return GetProjection(i, l, abc, N) * GetProjection(j, m, abc, N)
         - GetProjection(i, j, abc, N) * GetProjection(l, m, abc, N) / 2.;
}

}; // namespace sledgehamr
