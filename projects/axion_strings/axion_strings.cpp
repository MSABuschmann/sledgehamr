#include "axion_strings.h"

namespace axion_strings{

void axion_strings::Init() {
    cosmo.Init(this);
}

bool axion_strings::CreateLevelIf(const int lev, const double time) {
    return cosmo.CreateLevelIf(lev, time);
}

}; // namespace axion_strings
