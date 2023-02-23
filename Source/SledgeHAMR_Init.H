#ifndef SledgeHAMR_Init_H_
#define SledgeHAMR_Init_H_

class SledgeHAMR_Init
{
public:
	SledgeHAMR_Init ();
	
private:
	void determine_project_id ();
	void finish_AMReX_setup ();

	std::string project_name;
};

#endif // SledgeHAMR_Init_H_
