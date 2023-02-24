#ifndef LevelData_H_
#define LevelData_H_

/** @brief Struct that holds the MultiFab data while 
 *  	   keeping track of time
 */
struct LevelData
    : public amrex::MultiFab
{
public:
	double t;
}

#endif // LevelData_H_
