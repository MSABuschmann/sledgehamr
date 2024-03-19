#include <filesystem>

#include <AMReX_VisMF.H>
#include <AMReX_PlotFileUtil.H>
#include <AMReX_PlotFileDataImpl.H>

#include "io_module.h"
#include "hdf5_utils.h"
#include "sledgehamr_utils.h"
#include "slices.h"
#include "level_writer.h"
#include "amrex_plotfile.h"

namespace sledgehamr {

/** @brief Constructor will read relevant input parameters, create output
 *         folder(s), and set up all pre-defined output types.
 * @param  owner   Pointer to the simulation.
 */
           
IOModule::IOModule(Sledgehamr* owner) : sim(owner) {
    ParseParams();
    CheckIfOutputAlreadyExists(output_folder);
    CheckIfOutputAlreadyExists(alternative_output_folder);
    amrex::ParallelDescriptor::Barrier();
    CreateOutputFolder(output_folder);
    CreateOutputFolder(alternative_output_folder);
    AddOutputModules();
}

/** @brief Parsing of all relevant input parameters.
 */
void IOModule::ParseParams() {
    amrex::ParmParse pp("");
    pp.get("output.output_folder", output_folder);
    pp.query("output.alternative_output_folder", alternative_output_folder);

    std::string param_name = "output.checkpoints.rolling";
    pp.query(param_name.c_str(), rolling_checkpoints);
    utils::ErrorState validity = rolling_checkpoints ?
                utils::ErrorState::WARNING : utils::ErrorState::OK;
    std::string warning_msg = "Only the latest checkpoint will be kept.";
    utils::AssessParam(validity, param_name, rolling_checkpoints,
                       "", warning_msg, sim->nerrors, sim->do_thorough_checks);

    param_name = "input.delete_restart_checkpoint";
    pp.query(param_name.c_str(), delete_restart_checkpoint);
    validity = delete_restart_checkpoint ?
                utils::ErrorState::WARNING : utils::ErrorState::OK;
    warning_msg = "Restart checkpoint will be deleted!";
    utils::AssessParam(validity, param_name, delete_restart_checkpoint,
                       "", warning_msg, sim->nerrors, sim->do_thorough_checks);
}

/** @brief Will set up all pre-defined output types.
 */
void IOModule::AddOutputModules() {
    idx_slices = output.size();
    output.emplace_back("slices", OUTPUT_FCT(IOModule::WriteSlices));

    idx_coarse_box = output.size();
    output.emplace_back("coarse_box", OUTPUT_FCT(IOModule::WriteCoarseBox));

    idx_full_box = output.size();
    output.emplace_back("full_box", OUTPUT_FCT(IOModule::WriteFullBox));

    idx_slices_truncation_error = output.size();
    output.emplace_back("slices_truncation_error",
                        OUTPUT_FCT(IOModule::WriteSlicesTruncationError));

    idx_coarse_box_truncation_error = output.size();
    output.emplace_back("coarse_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteCoarseBoxTruncationError));

    idx_full_box_truncation_error = output.size();
    output.emplace_back("full_box_truncation_error",
                        OUTPUT_FCT(IOModule::WriteFullBoxTruncationError));

    idx_projections = output.size();
    output.emplace_back("projections", OUTPUT_FCT(IOModule::WriteProjections));

    idx_spectra = output.size();
    output.emplace_back("spectra", OUTPUT_FCT(IOModule::WriteSpectra));

    idx_gw_spectra = output.size();
    output.emplace_back("gw_spectra",
                        OUTPUT_FCT(IOModule::WriteGravitationalWaveSpectrum));

    idx_amrex_plotfile = output.size();
    output.emplace_back("amrex_plotfile",
                        OUTPUT_FCT(IOModule::WriteAmrexPlotFile));

    idx_performance_monitor = output.size();
    output.emplace_back("performance_monitor",
                        OUTPUT_FCT(IOModule::WritePerformanceMonitor));

    // Checkpoint. Always add checkpoints last.
    idx_checkpoints = output.size();
    output.emplace_back("checkpoints", OUTPUT_FCT(IOModule::WriteCheckpoint));

    bool write_at_start = false;
    amrex::ParmParse pp("");
    std::string param_name = "output.write_at_start";
    pp.query(param_name.c_str(), write_at_start);
    utils::AssessParamOK(param_name, write_at_start, sim->do_thorough_checks);
    if (!write_at_start) {
        for (OutputModule& out : output) {
            out.SetLastTimeWritten( sim->t_start );
        }
    }
}

/** @brief Will check if an output folder already exists and throws an error
 *         if we do not explicitly want to restart or rename any existing 
 *         output.
 * @param   folder  Folder to be checked. If string is empty we do nothing.
 */
void IOModule::CheckIfOutputAlreadyExists(std::string folder) {
    if (folder=="")
        return;

    // Determine and create output folder.
    bool rename_old = false;
    amrex::ParmParse pp("");
    std::string param_name = "output.rename_old_output";
    pp.query(param_name.c_str(), rename_old);

    utils::ErrorState validity = (utils::ErrorState)(
            !(amrex::FileExists(folder) &&
              !sim->restart_sim &&
              !rename_old));
    std::string error_msg =
                "Output folder " + folder + " already exists! "
                "If you intended to restart the simulation from the latest "
                "checkpoint within this folder please add 'input.restart = 1' "
                "to your input file. Otherwise please choose a different "
                "directory or set output.rename_old_output = 1";
    utils::AssessParam(validity, param_name, rename_old, error_msg, "",
                       sim->nerrors, sim->do_thorough_checks);

    if (validity == utils::ErrorState::ERROR &&
        !sim->do_thorough_checks) {
        amrex::Abort();
    }
}

/** @brief Will create an output folder if we are not restarting a sim. Will
 *         rename any already existing folder to *.old.xxxxxxxxx.
 * @param   folder  Folder to be created.
 */
void IOModule::CreateOutputFolder(std::string folder) {
    if (folder == "" || sim->do_thorough_checks)
        return;

    if (!sim->restart_sim) {
        std::string tmp = folder;
        while(tmp.back() == '/') {
            tmp.pop_back();
        }

        amrex::Print() << "Create output folder: " << tmp << std::endl;
        amrex::UtilCreateCleanDirectory(tmp, true);
    }
}

/** @brief Will go through all output types and trigger a write if criteria
 *         are met.
 * @param   force   Will force a write even if writing interval has not been
 *                  reached yet.
 */
void IOModule::Write(bool force) {
    // Make sure checkpoints are written last to ensure most up-to-date meta
    // data.
    for(int i = 0; i < output.size(); ++i) {
        if (i != idx_checkpoints) {
            sim->performance_monitor->Start(
                    sim->performance_monitor->idx_output, i);
            output[i].Write(sim->grid_new[0].t, force);
            sim->performance_monitor->Stop(
                    sim->performance_monitor->idx_output, i);
        }
    }

    sim->performance_monitor->Start(
            sim->performance_monitor->idx_output, idx_checkpoints);
    output[idx_checkpoints].Write(sim->grid_new[0].t, force);
    sim->performance_monitor->Stop(
            sim->performance_monitor->idx_output, idx_checkpoints);
}

/** @brief Will write slices on all levels.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteSlices(double time, std::string prefix) {
    Slices slices(sim, prefix, false);
    slices.Write();
    return true;
}

/** @brief Will write slices of truncation error estimates on all levels.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteSlicesTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    Slices slices(sim, prefix, true);
    slices.Write();
    return true;
}

/** @brief Will write scalar fields on the coarse level.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteCoarseBox(double time, std::string prefix) {
    LevelWriter writer(sim, prefix, idx_coarse_box);
    writer.Write();
    return true;
}

/** @brief Will write truncation error estimates on the coarse level.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteCoarseBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_coarse_box_truncation_error);
    writer.Write();
    return true;
}

/** @brief Will write scalar fields on all levels.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteFullBox(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_full_box);
    writer.Write();
    return true;
}

/** @brief Will write truncation error estimates on all levels.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteFullBoxTruncationError(double time, std::string prefix) {
    if (!sim->grid_old[0].contains_truncation_errors)
        return false;

    LevelWriter writer(sim, prefix, idx_full_box_truncation_error);
    writer.Write();
    return true;
}

/** @brief Will write an AMReX plotfile for yt support.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteAmrexPlotFile(double time, std::string prefix) {
    AmrexPlotFile writer(sim, prefix);
    writer.Write();
    return true;
}

/** @brief Will write any projections.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteProjections(double time, std::string prefix) {
    if (projections.empty())
        return false;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/projections.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    for (int p = 0; p < projections.size(); ++p)
        projections[p].Compute(p, file_id, sim);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

/** @brief Will write any spectra.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteSpectra(double time, std::string prefix) {
    if (spectra.empty())
        return false;

    sim->ReadSpectrumKs();

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    for (int p = 0; p < spectra.size(); ++p)
        spectra[p].Compute(p, file_id, sim);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

/** @brief Will write a gravitational wave spectrum.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteGravitationalWaveSpectrum(double time, std::string prefix) {
    if (!sim->with_gravitational_waves)
        return false;

    sim->ReadSpectrumKs();

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/spectra.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    sim->gravitational_waves->ComputeSpectrum(file_id);

    if (amrex::ParallelDescriptor::IOProcessor())
        H5Fclose(file_id);

    return true;
}

/** @brief Will write a checkpoint.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WriteCheckpoint(double time, std::string prefix) {
    Checkpoint chk(sim, prefix);
    chk.Write();

    if (rolling_checkpoints) {
        if (old_checkpoint != "") {
            Checkpoint chk_del(sim, old_checkpoint);
            chk_del.Delete();
        }

        old_checkpoint = prefix;
    }

    return true;
}

/** @brief Prints and logs the latest performance update.
 * @param   time    Current time.
 * @param   prefix  Assigned output folder.
 * @return  Whether the write was successfull or not.
 */
bool IOModule::WritePerformanceMonitor(double time, std::string prefix) {
    if (!sim->performance_monitor->IsActive())
        return false;

    hid_t file_id;
    if (amrex::ParallelDescriptor::IOProcessor()) {
        std::string filename = prefix + "/log.hdf5";
        file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    }

    sim->performance_monitor->Log(file_id);

    if (amrex::ParallelDescriptor::IOProcessor()) {
        H5Fclose(file_id);
    }

    return true;
}

/** @brief Writes an amrex::BoxArray to an hdf5 file.
 * @param   ba  Box array.
 */
void IOModule::WriteBoxArray(amrex::BoxArray& ba) {
    if (!amrex::ParallelDescriptor::IOProcessor())
        return;

    const int nba = ba.size();
    std::vector<int> x0(nba), y0(nba), z0(nba), x1(nba), y1(nba), z1(nba);
    for (int b = 0; b < nba; ++b) {
        x0[b] = ba[b].smallEnd(0);
        y0[b] = ba[b].smallEnd(1);
        z0[b] = ba[b].smallEnd(2);
        x1[b] = ba[b].bigEnd(0);
        y1[b] = ba[b].bigEnd(1);
        z1[b] = ba[b].bigEnd(2);
    }
    int header[1] = {nba};

    std::string filename = output_folder + "/box_layout.h5";
    hid_t file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT,
                            H5P_DEFAULT);
    utils::hdf5::Write(file_id, "header", header, 1);
    utils::hdf5::Write(file_id, "x0", x0.data(), nba);
    utils::hdf5::Write(file_id, "y0", y0.data(), nba);
    utils::hdf5::Write(file_id, "z0", z0.data(), nba);
    utils::hdf5::Write(file_id, "x1", x1.data(), nba);
    utils::hdf5::Write(file_id, "y1", y1.data(), nba);
    utils::hdf5::Write(file_id, "z1", z1.data(), nba);
    H5Fclose(file_id);
}


/** @brief Will restart a sim from a checkpoint after locating the appropiate
 *         one.
 */
void IOModule::RestartSim() {
    amrex::ParmParse pp("input");
    std::string selected_chk = "None Selected";
    pp.query("select_checkpoint", selected_chk);

    if (selected_chk == "None Selected") {
        int latest     = FindLatestCheckpoint(output_folder);
        int latest_alt = FindLatestCheckpoint(alternative_output_folder);

        if (latest > latest_alt) {
            initial_chk = output_folder + "/checkpoints/"
                        + std::to_string(latest);
        } else {
            initial_chk = alternative_output_folder + "/checkpoints/"
                        + std::to_string(latest_alt);
        }

        if (std::max(latest, latest_alt) == -1 &&
            amrex::ParallelDescriptor::IOProcessor()) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "No checkpoint found!";
            amrex::Abort(msg);
        }
   } else {
        if (amrex::is_integer(selected_chk.c_str())) {

            initial_chk = output_folder + "/checkpoints/" + selected_chk;

            if (!amrex::FileExists(initial_chk)) {
                initial_chk = alternative_output_folder + "/checkpoints/"
                            + selected_chk;
            }
        } else {
            initial_chk = selected_chk;
        }

        if (!amrex::FileExists(initial_chk)) {
            const char* msg = "Sledgehamr::IOModule::RestartSim: "
                              "Selected checkpoint not found!";
            amrex::Abort(msg);
        }
    }

    amrex::ParallelDescriptor::Barrier();
    Checkpoint chk(sim, initial_chk);
    chk.Read();

    if (delete_restart_checkpoint) {
        old_checkpoint = initial_chk;
    }
}

/** @brief After we restarted a sim from a checkpoint this procedure will 
 *         update the meta data of all output modules to keep intervals and
 *         counters consistent between runs.
 */
void IOModule::UpdateOutputModules() {
    Checkpoint chk(sim, initial_chk);
    chk.UpdateOutputModules();
}

/** @brief Locates the latest checkpoint within a parent folder.
 * @param   folder  Parent folder containing checkpoints.
 * @return  ID of the latest checkpoint.
 */
int IOModule::FindLatestCheckpoint(std::string folder) {
    if (folder == "")
        return -1;

    std::string prefix = folder + "/checkpoints/";
    std::vector<std::string> folders = GetDirectories(prefix);

    double latest_time = -DBL_MAX;
    int latest_chk = -1;
    for (std::string& str : folders) {
        str.erase(0, prefix.size());

        if (amrex::is_integer(str.c_str())) {
            std::string chk_folder = prefix + str;

            Checkpoint chk(sim, chk_folder);
            if (!chk.ReadHeader())
                continue;

            double time = chk.GetTime();
            if (time > latest_time) {
                latest_time = time;
                latest_chk = std::stoi(str);
            }
        }
    }

    return latest_chk;
}

/** @brief Obtains all directories within a parent directory.
 * @param   prefix  Parent directoy.
 * @return  Vector of all directories within parent directory.
 */
std::vector<std::string> IOModule::GetDirectories(const std::string prefix) {
    std::vector<std::string> res;
    for (auto& p : std::filesystem::recursive_directory_iterator(prefix))
        if (p.is_directory())
            res.push_back(p.path().string());
    return res;
}

}; // namespace sledgehamr
