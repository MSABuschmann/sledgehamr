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
    OutputModule(std::string output_prefix, output_fct function,
                 double write_interval, bool is_forceable=true);

    /** @brief Does the actual writing if criteria are met.
     * @param   time    Current time.
     * @param   force   Output will be written independent of the current time
     *                  interval if forceable=true.
     */
    void Write(double time, bool force=false);

private:
    /** @brief Function pointer to function that does the actual writing.
     */
    output_fct fct;

    /** @brief Next output index. Will be iterated +1 after each writing.
     */
    int next_id = 0;

    /** @brief Time at which this output has been written last.
     */
    double last_written = -1e99;

    /** @brief Interval at which this output shall be written.
     */
    double interval;

    /** @brief Prefix of the output.
     */
    std::string prefix;

    /** @brief Defines whether we can force a write, e.g. at the end of the
     *         simulation. Default is forceable=true.
     */
    bool forceable;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUTMODULE_H_
