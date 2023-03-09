template <typename T> 
void IOModule::ReadFromHDF5(std::string filename, std::vector<std::string> dnames,
					T* data)
{
	// Identify datatype
	hid_t mem_type_id;
	if( std::is_same<T, float>::value ){
		mem_type_id = H5T_NATIVE_FLOAT;
	}else if( std::is_same<T, double>::value ){
		mem_type_id = H5T_NATIVE_DOUBLE;
	}else if( std::is_same<T, int>::value ){
		mem_type_id = H5T_NATIVE_INT;
	}

	// Try and open HDF5 file
	hid_t file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
	if( file_id == H5I_INVALID_HID ){
		amrex::Abort("#error: Could not open initial state file: " + filename);
	}

	// Try and find dataset. Iterate over vector, use first to be found.
	std::string dname_conc = " |";
	std::string dname_found = "";
	for(std::string dname : dnames){
		htri_t exists = H5Lexists(file_id, dname.c_str(), H5P_DEFAULT);

		if( exists > 0 ){
			dname_found = dname;
			break;
		}

		dname_conc += " " + dname;
	}

	if( dname_found  == "" ){
		H5Fclose(file_id);
		amrex::Abort("#error: Could not find correct dataset in initial state file: " + filename + dname_conc);
	}

	// Read dataset
	hid_t dataset_id = H5Dopen2(file_id, dname_found.c_str(), H5P_DEFAULT);
	H5Dread(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(dataset_id);
	H5Fclose(file_id);
}

template <typename T>
void IOModule::WriteToHDF5(hid_t file_id, std::string dset, T* data, 
				 unsigned long long size)
{
	// Identify datatype.
	hid_t mem_type_id, dset_type_id;
	if( std::is_same<T, float>::value ){
		mem_type_id = H5T_NATIVE_FLOAT;
		dset_type_id = H5T_IEEE_F32LE;
	}else if( std::is_same<T, double>::value ){
		mem_type_id = H5T_NATIVE_DOUBLE;
		dset_type_id = H5T_IEEE_F64LE;
	}else if( std::is_same<T, int>::value ){
		mem_type_id = H5T_NATIVE_INT;
		dset_type_id = H5T_IEEE_F64LE;
	}else{
		amrex::Abort("#error: Writing of dataset " + dset + " failed due to unknown datatype."); 
	}

	// Create and write dataset.
	hsize_t dims[1] = {size};
	hid_t space = H5Screate_simple(1,dims, NULL);
	hid_t dataset_id = H5Dcreate(file_id, dset.c_str(), dset_type_id, space,H5P_DEFAULT,H5P_DEFAULT,H5P_DEFAULT);
	H5Dwrite(dataset_id, mem_type_id, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
	H5Dclose(dataset_id);
}
