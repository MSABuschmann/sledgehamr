#include "output_module.h"

namespace sledgehamr {

OutputModule::OutputModule(std::string output_prefix, output_fct function,
                           double write_interval, bool is_forceable)
    : prefix(output_prefix), fct(function), interval(write_interval),
      forceable(is_forceable) {
    /* TODO: Figure out at what index to start */

    amrex::UtilCreateDirectory(prefix.c_str(), 0755);
}

void OutputModule::Write(double time, bool force) {
    // Check if it is time to write output.
    /* TODO: Allow for custom interval and triggers */
    if( time == last_written ) return;
    if( time - last_written < interval && (!force && forceable) ) return;

    // Create output folder.
    std::string folder = prefix + "/" + std::to_string(next_id) + "/";
    amrex::UtilCreateDirectory(folder.c_str(), 0755);

    // Write.
    fct(time, folder);
    next_id++;
    last_written = time;
}

}; // namespace sledgehamr
