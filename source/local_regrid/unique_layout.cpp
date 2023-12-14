#include "unique_layout.h"

namespace sledgehamr {

/** @brief  Set up empty layout structure by creating lookup tables.
 * @param   local_regrid    Pointer to local regrid module.
 * @param   N               Number of potential boxes along one dimension
 *                          = dimN / blocking_factor.
 */
UniqueLayout::UniqueLayout(LocalRegrid* local_regrid, const int N)
    : Np(N),
      lr(local_regrid),
      mpi_n(amrex::ParallelDescriptor::NProcs()),
      mpi_mp(amrex::ParallelDescriptor::MyProc()) {
    p = std::make_unique<plane[]>(Np);
    uit Npn = mpi_n > Np ? 0 : Np/mpi_n;

    // Figure out who owns what.
    owner_of.resize(mpi_n);
    for (uit op = 0; op < mpi_n; ++op) {
        if (Npn == 0 && op < Np) {
            owner_of[op].push_back(op);
        } else {
            for (uit cp = 0; cp < Npn; ++cp)
                owner_of[op].push_back(cp + Npn*op);
        }
    }
    Np_this = owner_of[mpi_mp].size();
}

/** @brief  Add location.
 * @param   i   i-th index.
 * @param   j   j-th index.
 * @param   k   k-th index.
 */
void UniqueLayout::Add(const uit i, const uit j, const uit k) const {
    if (!p[i].count(j))
        p[i][j] = std::set<uit>();

    p[i][j].insert(k);
}

/** @brief Creates one unique layout structure from an array of layouts.
 * @param   uls Array of layouts.
 */
void UniqueLayout::Merge(std::vector<std::unique_ptr<UniqueLayout> >& uls) {
    // Merge all planes. Inner loop parallelized to maintain thread-safety.
    for (int l = 1; l < uls.size(); ++l) {
#pragma omp parallel for
        for (uit cp = 0; cp < Np; ++cp) {
            MergePlane(cp, &(uls[l]->p[cp]));
        }
    }
}

/** @brief Distributes any chunks that are not owned by this MPI rank but where
 *         added to it to the owning MPI rank. Will also receive any chunks by
 *         other ranks as well.
 */
void UniqueLayout::Distribute() {
    // Send/receive relevant planes according to communication matrix.
    for (int c = 1; c < mpi_n; ++c) {
        int op = lr->comm_matrix[mpi_mp][c];

        if (op < mpi_mp) {
            SendDistribution(op);
            RecvDistribution(op);
        } else {
            RecvDistribution(op);
            SendDistribution(op);
        }

        IncorporatePlanes();
    }

    // Clean up.
    for(uit cp = 0; cp < Np; ++cp){
        if( !Owns(cp) ){
            p[cp].clear();
        }
    }
}

/** @brief Empties the layout structure without destroying look-up tables.
 */
void UniqueLayout::Clear() {
    for (uit cp = 0; cp < Np; ++cp) {
        p[cp].clear();
    }
}

/** @brief Checks if location has been added to layout structure.
 * @param   i   i-th index of location.
 * @param   j   j-th index of location.
 * @param   k   k-th index of location.
 * @return  Whether location has been added.
 */
bool UniqueLayout::Contains(const uit i, const uit j, const uit k) const {
    if (p[i].count(j))
        return p[i][j].count(k);

    return false;
}

/** @brief Returns the number of locations added to the local chunks.
 */
int UniqueLayout::Size() {
    int size = 0;

    for (uit cp = 0; cp < Np_this; ++cp) {
        int i = owner_of[mpi_mp][cp];
        for (const std::pair<const uit,row>& n : p[i]) {
            size += n.second.size();
        }
    }

    return size;
}

/** @brief Returns the number of location added to all chunks even if not owned
 *         by this MPI ranks.
 */
int UniqueLayout::SizeAll() {
    int size = 0;

    for (uit cp = 0; cp < Np; ++cp) {
        for (const std::pair<const uit,row>& n : p[cp]) {
            size += n.second.size();
        }
    }

    return size;
}

/** @brief Constructs BoxList from layout structure.
 * @param   blocking_factor Blocking factor for BoxList.
 * @return BoxList.
 */
amrex::BoxList UniqueLayout::BoxList(const int blocking_factor) {
    amrex::BoxList bl;
    int i,j,k0,km;

    for (uit cp=0; cp<Np_this; ++cp) {
        i = owner_of[mpi_mp][cp];
        for (const std::pair<const uit, row>& n : p[i]) {
            j = n.first;

            k0=-1;
            for (const int &k : n.second) {
                if (k0 == -1) {
                    k0 = k;
                    km = k;
                } else if (k == km + 1) {
                     km = k;
                } else {
                    amrex::IntVect sm(i,   j,   k0);
                    amrex::IntVect bg(i+1, j+1, km+1);
                    bl.push_back(amrex::Box(sm*blocking_factor,
                                            bg*blocking_factor-1));
                    k0 = k;
                    km = k;
                }
            }

            if (k0 != -1) {
                amrex::IntVect sm(i,   j,   k0);
                amrex::IntVect bg(i+1, j+1, km+1);
                bl.push_back(amrex::Box(sm*blocking_factor,
                                        bg*blocking_factor-1));
            }
        }
    }

    return bl;
}

/** @brief Constructs BoxArray from layout structure.
 * @param   blocking_factor Blocking factor for BoxArray.
 * @return BoxArray.
 */
amrex::BoxArray UniqueLayout::BoxArray(const int blocking_factor) {
    return amrex::BoxArray(BoxList(blocking_factor));
}

/** @brief Merges individual planes to ensure uniqueness.
 * @param   cp  Current plane.
 * @param   pm  Plane to merge.
 */
void UniqueLayout::MergePlane(const uit cp, plane* pm) {
    p[cp].merge(*pm);

    for (const std::pair<const uit,row>& n : *pm) {
        row r;
        std::merge(p[cp][n.first].begin(), p[cp][n.first].end(),
                   n.second.begin(), n.second.end(),
                   std::inserter(r, r.begin()));
        p[cp][n.first] = r;
    }
}

/** @brief Wrap MPI rank index.
 * @param   i   MPI rank to warp.
 * @return Wrapped index.
 */
inline int UniqueLayout::Wrap(const int i) {
    return i < 0 ? i+mpi_n : i%mpi_n;
}

/** @brief Returns whether plane is owned by current rank.
 * @param   cp  Plane to check.
 * @return Ownership.
 */
bool UniqueLayout::Owns(const uit cp) {
    return (std::find(owner_of[mpi_mp].begin(), owner_of[mpi_mp].end(), cp) !=
            owner_of[mpi_mp].end());
}

/** @brief Sends everthing we have of a chunk owned by another MPI rank to this
 *         rank.
 * @param   op  Other MPI rank.
 */
void UniqueLayout::SendDistribution(const int op) {
    for (uit cpo = 0; cpo < owner_of[op].size(); ++cpo) {
        uit cp = owner_of[op][cpo];

        // For each plane determine number of rows and their length.
        std::vector<uit> r, length;
        uit total_length = 0;
        for (const std::pair<const uit,row>& n : p[cp]) {
            r.push_back( n.first );
            length.push_back( n.second.size() );
            total_length += n.second.size();
        }

        // Collapse all data.
        std::vector<uit> buf;
        buf.reserve(total_length);
        for (const std::pair<const uit,row>& n : p[cp]) {
            buf.insert(buf.end(), std::make_move_iterator(n.second.begin()),
                       std::make_move_iterator(n.second.end()));
        }

        SendVector(op, r);
        SendVector(op, length);
        SendVector(op, buf);
    }
}

/** @brief Receives everthing another MPI rank has about our local chunk.
 * @param   op  Other MPI rank.
 */
void UniqueLayout::RecvDistribution(const int op) {
    nps.clear();

    for (uit cpo=0; cpo<Np_this; ++cpo) {
        std::vector<uit> r = RecvVector(op);
        std::vector<uit> length = RecvVector(op);
        std::vector<uit> buf = RecvVector(op);

        // Unravel data according to metadata
        plane np;
        uit offset = 0;
        for (int cr = 0; cr < r.size(); ++cr) {
            std::set<uit> s{buf.begin() + offset,
                            buf.begin() + offset + length[cr]};
            np[r[cr]] = s;
            offset += length[cr];
        }

        nps.push_back(np);
    }
}

/** @brief  Sends data to another MPI rank.
 * @param   op  Other MPI rank.
 * @param   v   Data vector.
 */
inline void UniqueLayout::SendVector(const int op, const std::vector<uit>& v) {
    int size = v.size();
    MPI_Send(&size, 1, MPI_UNSIGNED, op, 501,
             amrex::ParallelDescriptor::m_comm);

    if (size > 0) {
        MPI_Send(&(v[0]), size, MPI_UNSIGNED_SHORT, op, 502,
                 amrex::ParallelDescriptor::m_comm);
    }
}

/** @brief Recevies data from another MPI rank.
 * @param   op  Other MPI rank.
 * @return  Received data vector.
 */
inline std::vector<uit> UniqueLayout::RecvVector(const int op) {
    int size;
    std::vector<uit> v;

    MPI_Recv(&size, 1, MPI_UNSIGNED, op, 501,
             amrex::ParallelDescriptor::m_comm, MPI_STATUS_IGNORE);

    if (size > 0) {
        v.resize(size);
        MPI_Recv(&(v[0]), size, MPI_UNSIGNED_SHORT, op, 502,
                 amrex::ParallelDescriptor::m_comm, MPI_STATUS_IGNORE);
    }

    return v;
}

/** @brief Merge received planes with our own local structure.
 */
void UniqueLayout::IncorporatePlanes() {
#pragma omp parallel for
    for (uit cp = 0; cp < Np_this; ++cp){
        MergePlane(owner_of[mpi_mp][cp], &(nps[cp]));
    }

    nps.clear();
}

}; // namespace sledgehamr
