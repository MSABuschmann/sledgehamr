#ifndef SLEDGEHAMR_UNIQUE_LAYOUT_H_
#define SLEDGEHAMR_UNIQUE_LAYOUT_H_

#include "local_regrid.h"

namespace sledgehamr {

class LocalRegrid;

typedef unsigned short uit;
typedef std::set<uit> row;
typedef std::unordered_map<uit,row> plane;

class UniqueLayout {
  public:
    UniqueLayout(LocalRegrid* local_regrid, const int N);
    ~UniqueLayout();

    void Add(const uit i, const uit j, const uit l) const;
    void Merge(std::vector<std::unique_ptr<UniqueLayout> >& uls);
    void Distribute();
    void Clear();

    bool Contains(const uit i, const uit j, const uit k) const;
    int Size();
    int SizeAll();
    amrex::BoxList BoxList(const int blocking_factor);
    amrex::BoxArray BoxArray(const int blocking_factor);

  private:
    void MergePlane(const uit cp, plane* pm);
    inline int Wrap(const int i);
    bool Owns(const uit cp);

    void SendDistribution(const int op);
    void RecvDistribution(const int op);
    inline void SendVector(const int op, const std::vector<uit>& v);
    inline std::vector<uit> RecvVector(const int op);
    void IncorporatePlanes();

    const uit Np;
    uit Np_this;
    std::vector< std::vector<uit> > owner_of;
    std::vector<plane> nps;
    plane* p;

    const int mpi_n;
    const int mpi_mp;

    LocalRegrid* lr;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_UNIQUE_LAYOUT_H_
