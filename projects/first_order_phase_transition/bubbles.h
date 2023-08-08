#ifndef PROJECTS_FIRST_ORDER_PHASE_TRANSITION_BUBBLE_H_
#define PROJECTS_FIRST_ORDER_PHASE_TRANSITION_BUBBLE_H_

namespace first_order_phase_transition {

class Bubble {
  public:
    double GetPos(const double D) const {
        return D >= p_bubble->L ? -1. : D * p_bubble->inv_dx;
    }

    int GetFinestLevel() const {
        return p_bubble->finest_level;
    }

    double GetVal(const int n, const int ind, const double frac) const {
        return p_bubble->profile[n][ind]   * (1.-frac) 
             + p_bubble->profile[n][ind+1] *     frac;
    }

    double GetLevel(const int ind) const {
        return p_bubble->level[ind];
    }
    
    int GetNBins() const {
        return p_bubble->profile.size();
    }

    bool operator<(const Bubble& other) const {
        return t < other.t;
    } 

    double x, y, z; // bubble location
    double t; // injection time
    std::vector<std::vector<double>> profile; // bubble profile
    std::vector<int> level; // bubble level
    double inv_dx; // inverse dx of profile binning
    double L; // profile cut-off (assume profile=0 beyond that)
    int finest_level; // finest level at which the profile is to be injected
    Bubble* p_bubble; // pointer to bubble which has the correct profile (can be this)
};

}

#endif // PROJECTS_FIRST_ORDER_PHASE_TRANSITION_BUBBLE_H_ 
