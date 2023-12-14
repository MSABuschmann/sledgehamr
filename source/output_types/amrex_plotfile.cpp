#include "amrex_plotfile.h"
#include <AMReX_PlotFileUtil.H>

namespace sledgehamr {

/** @brief Writes the AMReX output file.
 */
void AmrexPlotFile::Write() {
    std::string plotfilename = folder + "/output";
    int nlevels = sim->GetFinestLevel() + 1;

    amrex::Vector<const amrex::MultiFab*> r;
    amrex::Vector<int> level_steps;
    amrex::Vector<amrex::IntVect> ref_ratio;
    double time;
    for (int lev = 0; lev < nlevels; ++lev) {
        LevelData* ld = &sim->GetLevelData(lev);
        r.push_back(ld);
        level_steps.push_back(ld->istep);
        ref_ratio.emplace_back(2,2,2);
        time = ld->t;
    }

    amrex::Vector<std::string> varnames;
    for (int c = 0; c < r[0]->nComp(); ++c) {
        varnames.push_back(sim->GetScalarFieldName(c));
    }

    amrex::WriteMultiLevelPlotfile(plotfilename, nlevels, r, varnames,
                                   sim->GetGeometry(), time, level_steps,
                                   ref_ratio);
}

}; // namespace sledgehamr
