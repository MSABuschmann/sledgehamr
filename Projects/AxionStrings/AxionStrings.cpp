#include <AxionStrings.H>

AxionStrings::AxionStrings ()
{
	amrex::Print() << "Starting AxionStrings project..." << std::endl;
	amrex::Print() << "Added " << scalar_fields.size() << " fields:" << std::endl;
	for(int i=0;i<scalar_fields.size();++i){
		amrex::Print() << scalar_fields[i]->name << std::endl;
	}
	amrex::Print() << std::endl;
}
