#ifndef SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
#define SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_

#include "sledgehamr.h"

namespace sledgehamr {

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

    Sledgehamr* sim;
    std::string folder;
    int level_min;
    int level_max;
    const int output_id;
    std::string name;
    std::string info;
    bool with_truncation_errors;
    int downsample_factor = 1;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_LEVEL_WRITE_H_
