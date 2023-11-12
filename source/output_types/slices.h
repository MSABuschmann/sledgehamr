#ifndef SLEDGEHAMR_OUTPUT_TYPES_SLICES_H_
#define SLEDGEHAMR_OUTPUT_TYPES_SLICES_H_

#include "sledgehamr.h"

namespace sledgehamr {

class Slices {
  public:
    Slices(Sledgehamr* owner, std::string prefix,
           bool write_with_truncation_errors)
        : sim(owner),
          folder(prefix),
          with_truncation_errors(write_with_truncation_errors) {};

    void Write();

  private:
    void WriteSingleSlice(const LevelData* state, int lev, hid_t file_id,
                          std::string ident, int d1, int d2, int d3,
                          bool is_truncation_error);

    std::string folder;
    bool with_truncation_errors;
    Sledgehamr* sim;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_OUTPUT_TYPES_SLICES_H_
