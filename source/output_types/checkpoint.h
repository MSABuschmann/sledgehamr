#ifndef SLEDGEHAMR_CHECKPOINT_H_
#define SLEDGEHAMR_CHECKPOINT_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Writes or reads a checkpoint.
 */
class Checkpoint {
  public:
    /** @brief Constructor.
     * @param   owner       Pointer to simulation.
     * @param   chk_folder  Checkpoint directoy.
     */
    Checkpoint(Sledgehamr* owner, std::string chk_folder)
        : sim(owner), folder(chk_folder) {};

    void Write();
    void Read();
    bool ReadHeader();
    void UpdateOutputModules();
    void Delete();

    /** @brief Returns the time of the checkpoint.
     */
    double GetTime() const {
        return time;
    }

  private:
    static void GotoNextLine(std::istream& is);
    void UpdateLevels();

    /** @brief Returns the full path to the meta data header file.
     */
    std::string GetHeaderName() const {
        return folder + "/Meta.hdf5";
    }

    /** @brief Returns the full path to the box array data file.
     */
    std::string GetBoxArrayName() const {
        return folder + "/BoxArrays";
    }

    /** @brief Returns the full path to the individual level data folder.
     * @param   lev Level.
     */
    std::string GetLevelDirName(const int lev) const {
        return folder + "/Level_" + std::to_string(lev);
    }

    /** @brief Pointer to Simulation.
     */
    Sledgehamr* sim;

    /** @brief Checkpoint directory.
     */
    std::string folder;

    /** @brief Time of checkpoint.
     */
    double time;

    /** @brief Number of MPI ranks with which checkpoint was written.
     */
    int MPIranks;

    /** @brief Finest level in checkpoint.
     */
    int finest_level;

    /** @brief Coarse level grid cells along one axis in checkoint.
     */
    int dim0;

    /** @brief Number of ghost cells in checkpoint.
     */
    int nghost;

    /** @brief Number of scalar fields in checkpoint.
     */
    int nscalars;

    /** @brief Total number of output types with meta data.
     */
    int noutput;

    /** @brief Number of pre-defined output types with meta data.
     */
    int npredefoutput;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_CHECKPOINT_H_
