#include <filesystem>

#include "output_module.h"

namespace sledgehamr {

OutputModule::OutputModule(std::string output_prefix, std::string folder,
                           output_fct function, double write_interval,
                           bool is_forceable)
    : prefix(output_prefix), fct(function), interval(write_interval),
      name(folder), forceable(is_forceable) {
    std::string output_folder = prefix + "/" + name;
    if (!amrex::UtilCreateDirectory(output_folder, 0755)) {
        std::string msg = "sledgehamr::OutputModule::OutputModule: "
                          "Could not create output folder " + output_folder;
        amrex::Abort(msg);
    }
}

void OutputModule::Write(double time, bool force) {
    if (interval <= 0) return;

    // Check if it is time to write output.
    double t_now  = time_modifier(time);
    double t_last = time_modifier(last_written);

    if( t_now == t_last ) return;
    if( t_now - t_last < interval && (!force && forceable) ) return;

    // Create output folder.
    std::string folder = prefix + "/" + name + "/" + std::to_string(next_id);
    amrex::UtilCreateCleanDirectory(folder, true);
    folder += "/";

    // Attempt to write.
    if (fct(time, folder)) {
        next_id++;
        last_written = time;
    } else {
        std::filesystem::remove(folder);
    }
}

}; // namespace sledgehamr
