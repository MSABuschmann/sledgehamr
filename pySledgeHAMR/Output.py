from os import path
import numpy as np
import h5py

## This class reads all SledgeHAMR output
class Output:

	## Crawls the output and loads headers.
	def __init__(self, output_folder):
		self._prefix = output_folder
		self.__ParseFolderStructure()

	## Returns array of times at which slices have been written.
	def GetTimesOfSlices(self):
		print( self._slice_headers[:,0] )

	## Returns a slice.
	# @param	i		Number of slice to be read.
	# @param	direction	Direction of slice, e.g. 'x'.
	# @param	level		Which level should be returned.
	# @param	fields		List of scalar field names.
	# @return	d		Dictionary containing the time and slices.
	def GetSlices(self, i, direction, level, fields):
		# Get relevant parameters from header
		t = self._slice_headers[i,0]
		ranks = int(self._slice_headers[i,1])
		nlevels = int(self._slice_headers[i,2])
		dim0 = int(self._slice_headers[i,3])

		dim = dim0 * 2**level
		folder = self._prefix + '/slices'

		# Start dictionary	
		d = dict();
		d['t'] = t

		# Loop over fields, files, and boxes to reassemble data
		for s in fields:
			field = np.zero((dim,dim), dtype=np.float32)
			for f in range(ranks):
				file = folder + '/' + str(i)\
					+ '/Level_'+str(level)\
					+ '/' + str(f) + '.hdf5'

				fin = h5py.File(file,'r')
				le1 = np.array(fin['le1_'+direction], dtype='int')
				le2 = np.array(fin['le2_'+direction], dtype='int')
				he1 = np.array(fin['he1_'+direction], dtype='int')
				he2 = np.array(fin['he2_'+direction], dtype='int')

				for b in range(len(le1)):
					dset = s+'_'+direction+'_'+str(b+1)
					field[le1[b]:he1[b],le2[b]:he2[b]] =\
						(fin[dset][:]).reshape((he1[b]-le1[b], he2[b]-le2[b]))

			# Add result to dict.
			d[s] = field

		return d

	## Determines what output exists.
	def __ParseFolderStructure(self):
		self.__ParseSlices()

	## Determines how many slices have been written and when.
	## Reads header files.
	def __ParseSlices(self):
		folder = self._prefix + '/slices'
		i = 0

		# iterate over batches and read header of first files
		while True:
			file = folder + '/' + str(i) + '/Level_0/0.hdf5'
			if not path.exists( file ):
				break
			fin = h5py.File(file,'r')
			self._slice_headers.append( fin['Header_x'][:] )
			i = i + 1
		
		self._slice_headers = np.array( self._slice_headers )

		print('Number of slices found: ', len(self._slice_headers))
	
	## Path of output folder		
	_prefix = "."

	## List of headers files from slices
	_slice_headers = []
