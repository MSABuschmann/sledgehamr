#ifndef SLEDGEHAMR_SPECTRUM_H_
#define SLEDGEHAMR_SPECTRUM_H_

#include <AMReX_AmrCore.H>

#include "sledgehamr.h"

namespace sledgehamr {

class Sledgehamr;

typedef std::function<double(amrex::Array4<amrex::Real const> const&, const int,
        const int, const int, const int, const double, const double,
        const double, const std::vector<double>&)> spectrum_fct;

/** @brief Computes the spectrum given a quantity function and saves it to disk.
 */
class Spectrum {
  public:
    /** @brief Collects metadata.
     * @param   function        Integrand.
     * @param   identification  Unique identification string.
     */
    Spectrum(spectrum_fct function, std::string identification)
        : fct{function}, ident{identification} { };

    void Compute(const int id, const hid_t file_id, Sledgehamr* sim);

    static void Fft(const amrex::MultiFab& field, const int comp,
                    amrex::MultiFab& field_fft_real_or_abs,
                    amrex::MultiFab& field_fft_imag,
                    const amrex::Geometry& geom, bool abs);

    /** @brief Function pointer to integrand.
     */
    spectrum_fct fct;

    /** @brief Unique identification string.
     */
    std::string ident = "None";
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_SPECTRUM_H_
