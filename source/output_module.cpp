#include <filesystem>

#include "output_module.h"

namespace sledgehamr {

void OutputModule::CreateParentFolder(std::string this_prefix) {
    std::string output_folder = this_prefix + "/" + name;
    if (!amrex::UtilCreateDirectory(output_folder, 0755)) {
        std::string msg = "sledgehamr::OutputModule::CreateParentFolder: "
                          "Could not create output folder " + output_folder;
        amrex::Abort(msg);
    }
}

void OutputModule::Write(double time, bool force) {
    if (interval <= 0) return;

    // Check if it is time to write output.
    double t_now  = time_modifier(time);
    double t_last = time_modifier(last_written);

    if( t_now - t_last < interval && (!force && forceable) ) return;

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
    } else {
        std::filesystem::remove(folder);
    }
}

}; // namespace sledgehamr
