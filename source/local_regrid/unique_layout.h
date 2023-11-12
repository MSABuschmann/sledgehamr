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
    /** @brief
     * @param   N   Number of potential boxes along one dimension
     *              = dimN / blocking_factor.
     */
    UniqueLayout(LocalRegrid* local_regrid, const int N);
    ~UniqueLayout();

    /** @brief Add a single box.
     * @param   i   i-th index of box.
     * @param   j   j-th index of box.
     * @param   k   k-th index of box.
     */
    void Add(const uit i, const uit j, const uit l) const;

    /** @brief Merges a vector of UniqueLayouts on this node. Assumes
     *         uls[0] = this.
     */
    void Merge(std::vector<std::unique_ptr<UniqueLayout> >& uls);

    /** @brief Merges all UniqueLayouts across nodes.
     */
    void Distribute();

    /** @brief Remove all boxes.
     */
    void Clear();

    /** @brief Returns whether a box of given index has been added.
     */
    bool Contains(const uit i, const uit j, const uit k) const;

    /** @brief Returns number of added boxes within the chunk this node has
     *         has been assigned to.
     */
    int Size();

    /** @brief Returns number of all added boxes.
     */
    int SizeAll();

    /** @brief Returns an amrex::BoxList / amrex::BoxArray containing all boxes
     *         within the assigned chunk.
     * @param   blocking_factor The blocking_factor to be assued to convert from
     *          box index type to cell indices.
     */
    amrex::BoxList BoxList(const int blocking_factor);
    amrex::BoxArray BoxArray(const int blocking_factor);

  private:
    /** @brief Adds an entire plane.
     * @param   cp  Insertion index of plane.
     * @param   pm  Pointer to plane.
     */
    void MergePlane(const uit cp, plane* pm);

    /** @brief Wraps an MPI rank index.
     */
    inline int Wrap(const int i);

    /** @brief Returns whether the plane cp is owned by this node.
     */
    bool Owns(const uit cp);

    /** Sends boxes added to planes owned by node op to op.
     */
    void SendDistribution(const int op);

    /** Receives boxes added to planes owned by this node from node op.
     */
    void RecvDistribution(const int op);

    /** Sends a vector v to node op.
     */
    inline void SendVector(const int op, const std::vector<uit>& v);

    /** Receives a vector from node op.
     */
    inline std::vector<uit> RecvVector(const int op);

    /** Adds a set of planes received by other nodes.
     */
    void IncorporatePlanes();

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
    plane* p;

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
