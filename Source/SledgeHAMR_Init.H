#ifndef SledgeHAMR_Init_H_
#define SledgeHAMR_Init_H_

#include<SledgeHAMR.H>

class SledgeHAMR_Init
{
public:
	SledgeHAMR_Init ();
	SledgeHAMR* createInstance();
	
private:
	void determine_project_id ();
	void finish_AMReX_setup ();

	std::string project_name;
};

#endif // SledgeHAMR_Init_H_
