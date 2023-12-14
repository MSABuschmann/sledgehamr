#ifndef SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
#define SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_

#include "sledgehamr.h"

namespace sledgehamr {

/** @brief Writes level data to disk.
 */
class LevelWriter {
  public:
    LevelWriter(Sledgehamr* owner, std::string prefix, int output_type);
    void Write();

  private:
    void DetermineSetup();
    void ParseParams();
    void CheckDownsampleFactor();
    void WriteSingleLevel(const LevelData* state, int lev, hid_t file_id,
                          std::string ident, bool is_truncation_error);

    /** @brief Pointer to simulation.
     */
    Sledgehamr* sim;

    /** @brief Local output folder.
     */
    std::string folder;

    /** @brief Mininum level to write (inclusive).
     */
    int level_min;

    /** @brief Maximum level to write (inclusive).
     */
    int level_max;

    /** @brief Local unique output id. Must be either
     *          - io_module->idx_coarse_box
     *          - io_module->idx_coarse_box_truncation_error
     *          - io_module->idx_full_box
     *          - io_module->idx_full_box_truncation_error
     */
    const int output_id;

    /** @brief Display name of output type.
     */
    std::string name;

    /** @brief Information to print to screen for this output type.
     */
    std::string info;

    /** @brief Whether we are including truncation errors estimates.
     */
    bool with_truncation_errors;

    /** @brief Downsample data by this factor before writing.
     */
    int downsample_factor = 1;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
