#ifndef SLEDGEHAMR_IO_MODULE_H_
#define SLEDGEHAMR_IO_MODULE_H_

#include <AMReX_ParmParse.H>

#include "sledgehamr.h"
#include "output_module.h"
#include "projection.h"

namespace sledgehamr {

class Sledgehamr;
class Projection;
class Spectrum;

/** @brief Class that handles all I/O operations besides parsing the inputs
 *         file.
 */
class IOModule {
  public:
    IOModule (Sledgehamr* owner);

    /** @brief Writes output if requested.
     */
    void Write(bool force=false);

    void RestartSim();

    void UpdateOutputModules();

    void WriteBoxArray(amrex::BoxArray& ba);

    /** @brief Vectors containing instructions for projections and spectra.
     */
    std::vector<Projection> projections;
    std::vector<Spectrum> spectra;

    /** @brief Easy access pointers for user-modifications.
     */
    int idx_slices = -1;
    int idx_coarse_box = -1;
    int idx_full_box = -1;
    int idx_slices_truncation_error = -1;
    int idx_coarse_box_truncation_error = -1;
    int idx_full_box_truncation_error = -1;
    int idx_projections = -1;
    int idx_spectra = -1;
    int idx_gw_spectra = -1;
    int idx_performance_monitor = -1;
    int idx_checkpoints = -1;

    /** @brief Vector of output modules
     */
    std::vector<OutputModule> output;

    std::string output_folder;
    std::string alternative_output_folder = "";
    std::string old_checkpoint = "";

  private:
    bool WriteSlices(double time, std::string prefix);
    bool WriteSlicesTruncationError(double time, std::string prefix);
    bool WriteCoarseBox(double time, std::string prefix);
    bool WriteCoarseBoxTruncationError(double time, std::string prefix);
    bool WriteFullBox(double time, std::string prefix);
    bool WriteFullBoxTruncationError(double time, std::string prefix);
    bool WriteProjections(double time, std::string prefix);
    bool WriteSpectra(double time, std::string prefix);
    bool WriteGravitationalWaveSpectrum(double time, std::string prefix);
    bool WritePerformanceMonitor(double time, std::string prefix);
    bool WriteCheckpoint(double time, std::string prefix);
    int FindLatestCheckpoint(std::string folder);

    std::vector<std::string> GetDirectories(const std::string prefix);
    void CheckIfOutputAlreadyExists(std::string folder);
    void CreateOutputFolder(std::string folder);

    void AddOutputModules();
    void ParseParams();

    /** Downsampling factors for coarse/full level output.
     */
    int coarse_box_downsample_factor = 1;
    int coarse_box_truncation_error_downsample_factor = 1;
    int full_box_downsample_factor = 1;
    int full_box_truncation_error_downsample_factor = 1;

    std::string initial_chk = "";
    bool rolling_checkpoints = false;
    bool delete_restart_checkpoint = false;

    /** @brief Pointer to owner on whose data this class operates.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_IO_MODULE_H_
