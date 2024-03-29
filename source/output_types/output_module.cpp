#include <filesystem>

#include <AMReX_ParmParse.H>

#include "output_module.h"

namespace sledgehamr {

/** @brief Constructor that provides the settings.
 * @param   module_name     Name of output (to appear in file path e.g.).
 * @param   function        Function pointer to function that writes the
 *                          actual output.
 * @param   is_forcable     If true, output will be written at the very end
 *                          of the simulation independent of the time
 *                          interval.
 */
OutputModule::OutputModule(std::string module_name, output_fct function,
                           bool is_forceable)
    : fct(function),
      name(module_name),
      forceable(is_forceable) {
    ParseParams();
    CreateParentFolder(prefix);
    if (alternate)
        CreateParentFolder(alt_prefix);
}

/** @brief Parses all parameters related to this output type.
 */
void OutputModule::ParseParams() {
    amrex::ParmParse pp_out("output");
    pp_out.get("output_folder", prefix);
    pp_out.query("alternative_output_folder", alt_prefix);

    std::string pre = "output." + name;
    amrex::ParmParse pp(pre);
    pp.query("interval", interval);
    pp.query("alternate", alternate);
    pp.query("min_t", t_min);
    pp.query("max_t", t_max);

    if (alternate && alt_prefix == "") {
        std::string msg = "sledgehamr::OutputModule::ParseParams: "
                          "Alternating output selected but no alternative "
                          "output folder given";
        amrex::Abort(msg);
    }
}

/** @brief Creates a parent folder for the output.
 * @param this_prefix   parent folder.
 */
void OutputModule::CreateParentFolder(std::string this_prefix) {
    std::string output_folder = this_prefix + "/" + name;
    if (!amrex::UtilCreateDirectory(output_folder, 0755)) {
        std::string msg = "sledgehamr::OutputModule::CreateParentFolder: "
                          "Could not create output folder " + output_folder;
        amrex::Abort(msg);
    }
}

/** @brief Does the actual writing if criteria are met.
 * @param   time    Current time.
 * @param   force   Output will be written independent of the current time
 *                  interval if forceable=true.
 */
void OutputModule::Write(double time, bool force) {
    if (interval < 0) return;

    // Check if it is time to write output.
    double t_now  = time_modifier(time);
    double t_last = time_modifier(last_written);

    if (t_now > t_max || t_now < t_min) return;
    if (t_now - t_last < interval && (!force && forceable)) return;

    std::string this_prefix = (alternate && next_id%2 == 1) ?
                              alt_prefix : prefix;

    // Create output folder.
    std::string folder = this_prefix + "/" + name + "/"
                       + std::to_string(next_id);
    amrex::UtilCreateCleanDirectory(folder, true);
    folder += "/";


    // Attempt to write.
    if (fct(time, folder)) {
        next_id++;
        last_written = time;
        amrex::Print() << "Wrote " << name <<  ": " << folder << std::endl;
    } else {
        std::filesystem::remove(folder);
    }
}

}; // namespace sledgehamr
