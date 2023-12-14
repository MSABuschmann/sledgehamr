#ifndef SLEDGEHAMR_OUTPUTMODULE_H_
#define SLEDGEHAMR_OUTPUTMODULE_H_

#include <functional>
#include <string>

#include <AMReX_AmrCore.H>
#include <AMReX_Utility.H>

namespace sledgehamr {

/** @brief Macro and typedef to make the handling of function pointers easier.
 *         The function pointer will point to a function that does the actual
 *         writing of output.
 */
#define OUTPUT_FCT(fct) std::bind(&fct, this, std::placeholders::_1,\
                                  std::placeholders::_2)
typedef std::function<bool(double, std::string)> output_fct;

/** @brief Macro and typedef to make it possible to modify the time interval.
 */
#define TIME_FCT(fct) std::bind(&fct, this, std::placeholders::_1)
typedef std::function<double(double)> time_fct;

/** @brief This class handles the writing of an individual output format
 *         provided through a function pointer. It keeps track of timings to
 *         check if this particular format should be written at the current time
 *         or not.
 */
class OutputModule {
  public:
    OutputModule(std::string module_name, output_fct function,
                 bool is_forceable=true);
    void Write(double time, bool force=false);

    /** @brief Change the time interval to something arbitrary.
     * @param   mod Time modifier function.
     */
    void SetTimeFunction(time_fct mod) {
        time_modifier = mod;
    };

    /** @brief Default time interval is linear in time.
     * @param   time    Current time.
     */
    double DefaultInterval(double time) {
        return time;
    };

    /** @brief Change time interval width.
     */
    void SetInterval(double new_interval) {
        interval = new_interval;
    };

    /** @brief Get id of next output.
     */
    int GetNextId() const {
        return next_id;
    };

    /** @brief Sets the id of the next output.
     */
    void SetNextId(int id) {
        next_id = id;
    };

    /** @brief Sets the time we wrote output last.
     */
    void SetLastTimeWritten(double time) {
        last_written = time;
    };

    /** @brief Sets whether we want to alternate between output paths.
     */
    void Alternate(bool do_alternate) {
        alternate = do_alternate;
    }

    /** @brief Sets and create alternative output folder.
     */
    void SetAlternativePrefix(bool alternative_prefix) {
        alt_prefix = alternative_prefix;
        CreateParentFolder(alt_prefix);
    }

    /** @brief Returns the last time we wrote output.
     */
    double GetLastTimeWritten() const {
        return last_written;
    };

    /** @brief Returns the name of this output type.
     */
    std::string GetName() const {
        return name;
    };

private:
    void ParseParams();
    void CreateParentFolder(std::string this_prefix);

    /** @brief Function pointer to function that does the actual writing.
     */
    output_fct fct;

    /** @brief Time interval modifier function. By default linear in time.
     */
    time_fct time_modifier = TIME_FCT(OutputModule::DefaultInterval);

    /** @brief Next output index. Will be iterated +1 after each writing.
     */
    int next_id = 0;

    /** @brief Time at which this output has been written last.
     */
    double last_written = -DBL_MAX;

    /** @brief Interval at which this output shall be written.
     */
    double interval = -1;

    /** @brief Prefix of the output folder.
     */
    std::string prefix;

    /** @brief Prefix of alternative output folder.
     */
    std::string alt_prefix = "";

    /** @brief Whether we want to alternate between two output folders.
     */
    bool alternate = false;

    /** @brief Minimum time before we start writing.
     */
    double t_min = -DBL_MAX;

    /** @brief Maximum time before we start writing.
     */
    double t_max =  DBL_MAX;

    /** @brief Defines whether we can force a write, e.g. at the end of the
     *         simulation. Default is forceable=true.
     */
    bool forceable;

    /** @brief Output name.
     */
    std::string name = "Unknown";
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUTMODULE_H_
