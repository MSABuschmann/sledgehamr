#ifndef SLEDGEHAMR_FILL_LEVEL_H_
#define SLEDGEHAMR_FILL_LEVEL_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief This class will fill an entire level with data from either a
 *         checkpoint, an hdf5 file, an array or using a constant value,
 *         depending on the chosen routine. Operates on a single level given in
 *         the constructor unless we are trying to initialize a run via a
 *         checkpoint.
 */
class FillLevel {
  public:
    /** @brief Constructor. Doesn't do any work except collect meta data.
     * @param   owner   Pointer to the simulation.
     * @param   level   Number of level to be filled with data (unless we
     *                  initialize from a checkpoint).
     */
    FillLevel(Sledgehamr* owner, const int level) : sim(owner), lev(level) {};

    void FromInitialStateFile();
    void FromCheckpoint(std::string folder);
    void FromHdf5File(std::string initial_state_file);
    void FromArray(const int comp, double* data, const long long dimN);
    void FromArrayChunks(const int comp, double* data);
    void FromConst(const int comp, const double c);

  private:
    /** @brief Pointer to the simulation.
     */
    Sledgehamr* sim;

    /** @brief Number of level to be filled with data (unless we initialize from
     *         a checkpoint).
     */
    const int lev;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_FILL_LEVEL_H_
