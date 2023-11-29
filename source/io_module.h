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

/** @brief Class that handles all I/O operations.
 */
class IOModule {
  public:
    IOModule (Sledgehamr* owner);

    void Write(bool force=false);
    void RestartSim();
    void UpdateOutputModules();
    void WriteBoxArray(amrex::BoxArray& ba);

    /** @brief Vectors containing instructions for projections.
     */
    std::vector<Projection> projections;

    /** @brief Vector containing instructions for spectrum computations.
     *         Gravitational wave spectra are handled separately and no
     *         instructions need to be provided.
     */
    std::vector<Spectrum> spectra;

    /** @brief Output module ID's of various output types.
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
    int idx_amrex_plotfile = -1;
    int idx_checkpoints = -1;

    /** @brief Vector of output modules.
     */
    std::vector<OutputModule> output;

    /** @brief Path to output folder.
     */
    std::string output_folder;

    /** @brief Path to alternative output folder if one is provided.
     */
    std::string alternative_output_folder = "";

    /** @brief Path to last checkpoint. Needed so we can delete it if
     *         we need to.
     */
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
    bool WriteAmrexPlotFile(double time, std::string prefix);
    bool WriteCheckpoint(double time, std::string prefix);

    int FindLatestCheckpoint(std::string folder);
    std::vector<std::string> GetDirectories(const std::string prefix);
    void CheckIfOutputAlreadyExists(std::string folder);
    void CreateOutputFolder(std::string folder);

    void AddOutputModules();
    void ParseParams();

    /** @brief Path to initial checkpoint file if any.
     */
    std::string initial_chk = "";

    /** @brief If we are using rolling checkpoints, i.e. only ever keep the
     *         latest one.
     */
    bool rolling_checkpoints = false;

    /** @brief If we want to delete the checkpoint we restarted the sim from.
     *         Will only be deleted once we created a newer one.
     */
    bool delete_restart_checkpoint = false;

    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_IO_MODULE_H_
