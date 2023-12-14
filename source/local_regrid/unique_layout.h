#ifndef SLEDGEHAMR_UNIQUE_LAYOUT_H_
#define SLEDGEHAMR_UNIQUE_LAYOUT_H_

#include "local_regrid.h"

namespace sledgehamr {

class LocalRegrid;

typedef unsigned short uit;
typedef std::set<uit> row;
typedef std::unordered_map<uit,row> plane;

/** @brief This class ensures the uniqueness of a grid by reducing the problem
 *         to sets of touples. Each core can work independently on its own
 *         UniqueLayout which can then be efficiently merged across cores and
 *         nodes to a global unique grid. Memory consumption scales with the
 *         number of boxes added, not with the potential number of boxes.
 */
class UniqueLayout {
  public:
    UniqueLayout(LocalRegrid* local_regrid, const int N);

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
    void IncorporatePlanes();
    inline int Wrap(const int i);
    bool Owns(const uit cp);

    void SendDistribution(const int op);
    void RecvDistribution(const int op);
    inline void SendVector(const int op, const std::vector<uit>& v);
    inline std::vector<uit> RecvVector(const int op);

    /** Total number of planes.
     */
    const uit Np;

    /** Number of planes owned by this node.
     */
    uit Np_this;

    /** List of owners of each plane.
     */
    std::vector< std::vector<uit> > owner_of;

    /** Vector of planes that require incorporation.
     */
    std::vector<plane> nps;

    /** Array of Np planes.
     */
    std::unique_ptr<plane[]> p;

    /** Number of nodes.
     */
    const int mpi_n;

    /** Own MPI rank.
     */
    const int mpi_mp;

    /** Pointer to LocalRegrid module.
     */
    LocalRegrid* lr;
};

}; // namespace sledgehamr

#endif // SLEDGEHAMR_UNIQUE_LAYOUT_H_
