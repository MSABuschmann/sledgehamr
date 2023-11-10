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

#define TIME_FCT(fct) std::bind(&fct, this, std::placeholders::_1)
typedef std::function<double(double)> time_fct;

/** @brief This class handles the writing of an individual output format
 *         provided through a function pointer. It keeps track of timings to
 *         check if this particular format should be written at the current time
 *         or not.
 */
class OutputModule {
  public:
    /** @brief Constructor that provides the settings.
     * @param   output_prefix   Output path. Output will be written to
     *                          output_prefix/$id/.
     * @param   function        Function pointer to function that writes the
     *                          actual output.
     * @param   write_interval  Interval at which the output shall be written.
     * @param   is_forcable     If true, output will be written at the very end
     *                          of the simulation independent of the time
     *                          interval.
     */
    OutputModule(std::string output_prefix, std::string folder,
                 output_fct function, double write_interval,
                 bool is_forceable=true)
        : prefix(output_prefix),
          fct(function),
          interval(write_interval),
          name(folder),
          forceable(is_forceable) {
        CreateParentFolder(prefix);
    }

    /** @brief Does the actual writing if criteria are met.
     * @param   time    Current time.
     * @param   force   Output will be written independent of the current time
     *                  interval if forceable=true.
     */
    void Write(double time, bool force=false);

    void SetTimeFunction(time_fct mod) {
        time_modifier = mod;
    };

    double DefaultInterval(double time) {
        return time;
    };

    void SetInterval(double new_interval) {
        interval = new_interval;
    };

    int GetNextId() const {
        return next_id;
    };

    void SetNextId(int id) {
        next_id = id;
    };

    void SetLastTimeWritten(double time) {
        last_written = time;
    };

    void Alternate(bool do_alternate) {
        alternate = do_alternate;
    }

    void SetAlternativePrefix(bool alternative_prefix) {
        alt_prefix = alternative_prefix;
        CreateParentFolder(alt_prefix);
    }

    double GetLastTimeWritten() const {
        return last_written;
    };

    std::string GetName() const {
        return name;
    };

private:
    /** @brief Function pointer to function that does the actual writing.
     */
    output_fct fct;

    time_fct time_modifier = TIME_FCT(OutputModule::DefaultInterval);

    /** @brief Next output index. Will be iterated +1 after each writing.
     */
    int next_id = 0;

    /** @brief Time at which this output has been written last.
     */
    double last_written = -DBL_MAX;

    /** @brief Interval at which this output shall be written.
     */
    double interval;

    /** @brief Prefix of the output.
     */
    std::string prefix;

    std::string alt_prefix = "";

    bool alternate = false;

    /** @brief Defines whether we can force a write, e.g. at the end of the
     *         simulation. Default is forceable=true.
     */
    bool forceable;

    std::string name = "Unknown";
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUTMODULE_H_
