from os import path
import numpy as np
import h5py

## This class reads all SledgeHAMR output
class Output:
    ## Crawls the output and loads headers.
    # @param    output_folder   Folder containing the simulation output.
    def __init__(self, output_folder):
        self._prefix = output_folder

        self._slice_headers = []
        self._coarse_box_headers = []
        self._full_box_headers = []
        self._slice_truncation_error_headers = []
        self._coarse_box_truncation_error_headers = []
        self._full_box_truncation_error_headers = []
        self._projection_headers = []
        self._spectrum_headers = []
        self._gravitational_wave_spectrum_headers = []

        self.__ParseFolderStructure()

    ## TODO check for empty
    ## Returns array of times at which slices have been written.
    def GetTimesOfSlices(self):
        return self.__GetTimes(self._slice_headers)

    ## Returns array of times at which coarse boxes have been written.
    def GetTimesOfCoarseBoxes(self):
        return self.__GetTimes(self._coarse_box_headers)

    ## Returns array of times at which full boxes have been written.
    def GetTimesOfFullBoxes(self):
        return self.__GetTimes(self._full_box_headers)

    ## Returns array of times at which slices of truncation errors have been
    #  written.
    def GetTimesOfSlicesTruncationError(self):
        return self.__GetTimes(self._slice_truncation_error_headers)

    ## Returns array of times at which a coarse box of truncation errors has
    #  been written.
    def GetTimesOfCoarseBoxesTruncationError(self):
        return self.__GetTimes(self._coarse_box_truncation_error_headers)

    ## Returns array of times at which a full box of truncation errors has been
    #  written.
    def GetTimesOfFullBoxesTruncationError(self):
        return self.__GetTimes(self._full_box_truncation_error_headers)

    ## Returns array of times at which projections have been written.
    def GetTimesOfProjections(self):
        return self.__GetTimes(self._projection_headers)

    ## Returns array of times at which projections have been written.
    def GetTimesOfSpectra(self):
        return self.__GetTimes(self._spectrum_headers)

    ## Returns array of times at which projections have been written.
    def GetTimesOfGravitationalWaveSpectra(self):
        return self.__GetTimes(self._gravitational_wave_spectrum_headers)

    def __GetTimes(self, array):
        if len(array) > 0:
            return array[:,0]

    ## Returns a slice.
    # @param    i           Number of slice to be read.
    # @param    direction   Direction of slice, e.g. 'x'.
    # @param    level       Which level should be returned.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetSlice(self, i, direction, level, fields):
        # Get relevant parameters from header
        t = self._slice_headers[i,0]
        ranks = int(self._slice_headers[i,1])
        nlevels = int(self._slice_headers[i,2])
        dim0 = int(self._slice_headers[i,3])

        dim = dim0 * 2**level
        folder = self._prefix + '/slices/' + str(i) + '/Level_'+str(level)

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s] = self.__Read2dField(folder, dim, direction, ranks, s, 1)

        return d

    ## Returns a slice of truncation errors.
    # @param    i           Number of slice to be read.
    # @param    direction   Direction of slice, e.g. 'x'.
    # @param    level       Which level should be returned.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetSliceTruncationError(self, i, direction, level, fields):
        # Get relevant parameters from header
        t = self._slice_truncation_error_headers[i,0]
        ranks = int(self._slice_truncation_error_headers[i,1])
        nlevels = int(self._slice_truncation_error_headers[i,2])
        dim0 = int(self._slice_truncation_error_headers[i,3])

        dim = dim0 * 2**level // 2
        folder = self._prefix + '/slices_truncation_error/' + str(i)\
               + '/Level_'+str(level)

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s+'_truncation_error'] =\
                   self.__Read2dField(folder, dim, 'te_'+direction, ranks, s, 2)
            d[s] = self.__Read2dField(folder, dim*2, direction, ranks, s, 1)

        return d

    ## Returns a coarse box.
    # @param    i           Number of slice to be read.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetCoarseBox(self, i, fields):
        # Get relevant parameters from header
        t = self._coarse_box_headers[i,0]
        ranks = int(self._coarse_box_headers[i,1])
        dim0 = int(self._coarse_box_headers[i,3])
        downsample = int(self._coarse_box_headers[i,4])

        dim = int(dim0 / downsample)
        folder = self._prefix + '/coarse_box/' + str(i) + '/Level_0/'

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s] = self.__Read3dField(folder, dim, ranks, s, downsample, 'data')

        return d

    ## Returns a coarse box of truncation errors
    # @param    i           Number of slice to be read.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetCoarseBoxTruncationError(self, i, fields):
        # Get relevant parameters from header
        t = self._coarse_box_truncation_error_headers[i,0]
        ranks = int(self._coarse_box_truncation_error_headers[i,1])
        dim0 = int(self._coarse_box_truncation_error_headers[i,3])
        downsample = int(self._coarse_box_truncation_error_headers[i,4])

        dim = int(dim0 / downsample) // 2
        folder = self._prefix + '/coarse_box_truncation_error/' + str(i) \
               + '/Level_0/'

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s+'_truncation_error'] =\
                   self.__Read3dField(folder, dim, ranks, s, downsample*2, 'te')
            d[s] = self.__Read3dField(folder, dim*2, ranks, s, downsample,\
                                      'data')

        return d

    ## Returns a full box.
    # @param    i           Number of slice to be read.
    # @param    level       Level.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetFullBox(self, i, level, fields):
        # Get relevant parameters from header
        t = self._full_box_headers[i,0]
        ranks = int(self._full_box_headers[i,1])
        dim0 = int(self._full_box_headers[i,3])
        downsample = int(self._full_box_headers[i,4])

        dim = int(dim0 * 2**level / downsample)
        folder = self._prefix + '/full_box/'+str(i)+'/Level_' + str(level) + '/'

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s] = self.__Read3dField(folder, dim, ranks, s, downsample, 'data')

        return d

    ## Returns a full box.
    # @param    i           Number of slice to be read.
    # @param    level       Level.
    # @param    fields      List of scalar field names.
    # @return   d           Dictionary containing the time and slices.
    def GetFullBoxTruncationError(self, i, level, fields):
        # Get relevant parameters from header
        t = self._full_box_truncation_error_headers[i,0]
        ranks = int(self._full_box_truncation_error_headers[i,1])
        dim0 = int(self._full_box_truncation_error_headers[i,3])
        downsample = int(self._full_box_truncation_error_headers[i,4])

        dim = int(dim0 * 2**level / downsample) // 2
        folder = self._prefix + '/full_box_truncation_error/'+str(i)+'/Level_'\
                + str(level) + '/'

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s+'_truncation_error'] =\
                   self.__Read3dField(folder, dim, ranks, s, downsample*2, 'te')
            d[s] = self.__Read3dField(folder, dim*2, ranks, s, downsample,\
                                      'data')

        return d

    ## Returns a projection.
    # @param    i           Number of projection to be read.
    # @para     names       List of projection names
    # @return   d           Dictionary containing the time and projection.
    def GetProjection(self, i, names):
        # Get relevant parameters from header
        t = self._projection_headers[i,0]
        dim = int(self._projection_headers[i,1])
        file = self._prefix + '/projections/'+str(i)+'/projections.hdf5'

        # Start dictionary
        d = dict();
        d['t'] = t

        fin = h5py.File(file,'r')
        for s in names:
            d[s] = fin[s+'_data'][:].reshape((dim, dim))
        fin.close()
        return d

    ## Returns projections.
    # @param    i           Number of projection to be read.
    # @para     names       List of projection names
    # @return   d           Dictionary containing the time and projection.
    def GetProjectionN(self, i, names):
        # Get relevant parameters from header
        t = self._projection_headers[i,0]
        dim = int(self._projection_headers[i,1])
        file = self._prefix + '/projections/'+str(i)+'/projections.hdf5'

        # Start dictionary
        d = dict();
        d['t'] = t

        fin = h5py.File(file,'r')
        for s in names:
            d[s] = fin[s+'_n'][:].reshape((dim, dim))
        fin.close()
        return d

    ## Returns sprectra.
    # @param    i           State number.
    # @para     names       List of spectrum names.
    # @return   d           Dictionary containing the time and spectra.
    def GetSpectrum(self, i, names):
        # Get relevant parameters from header
        t = self._spectrum_headers[i,0]
        file = self._prefix + '/spectra/'+str(i)+'/spectra.hdf5'

        fin = h5py.File(file,'r')

        # Start dictionary
        d = dict();
        d['t'] = t
        d['k_sq'] = fin['k_sq'][:]

        for s in names:
            d[s] = fin[s][:]

        fin.close()
        return d

    ## Returns gravitational wave sprectra.
    # @param    i           State number.
    # @para     names       List of spectrum names.
    # @return   d           Dictionary containing the time and spectra.
    def GetGravitationalWaveSpectrum(self, i, spectype='gw_spectra'):
        # Get relevant parameters from header
        t = self._gravitational_wave_spectrum_headers[i,0]
        file = self._prefix + '/'+spectype+'/'+str(i)+'/spectra.hdf5'

        fin = h5py.File(file,'r')

        # Start dictionary
        d = dict();
        d['t'] = t
        d['k_sq'] = fin['k'][:]
        d['spectrum'] = fin['Spectrum'][:]

        fin.close()
        return d

    ## Helper function parsing a 2D field.
    # @param    folder      Folder containing the chunks.
    # @param    dim         Number of cells in each dimension.
    # @param    direction   'x', 'y', or 'z'.
    # @param    ranks       Number of MPI ranks = number of chunks.
    # @param    ident       Name of component.
    # @param    downsample  In case the field has been downsampled prior to
    #                       writing.
    # @return   2D field.
    def __Read2dField(self, folder, dim, direction, ranks, ident, downsample):
        field = np.ones((dim,dim), dtype=np.float32) * np.nan
        for f in range(ranks):
            file = folder + '/' + str(f) + '.hdf5'

            fin = h5py.File(file,'r')
            if 'le1_'+direction in fin.keys():
                le1 = np.array(fin['le1_'+direction], dtype='int') // downsample
                le2 = np.array(fin['le2_'+direction], dtype='int') // downsample
                he1 = np.array(fin['he1_'+direction], dtype='int') // downsample
                he2 = np.array(fin['he2_'+direction], dtype='int') // downsample

                for b in range(len(le1)):
                    dset = ident+'_'+direction+'_'+str(b+1)
                    field[le1[b]:he1[b],le2[b]:he2[b]] =\
                            (fin[dset][:]).reshape(\
                                    (he1[b]-le1[b], he2[b]-le2[b]))
            fin.close()
        return field

    ## Helper function parsing a 3D field.
    # @param    folder      Folder containing the chunks.
    # @param    dim         Number of cells in each dimension.
    # @param    ranks       Number of MPI ranks = number of chunks.
    # @param    ident       Name of component.
    # @param    downsample  In case the field has been downsampled prior to
    #                       writing.
    # @param    ident2      Dataset name, either "data" or "te" for truncation
    #                       error estimates.
    # @return   3D field.
    def __Read3dField(self, folder, dim, ranks, ident, downsample, ident2):
        field = np.ones((dim, dim, dim), dtype=np.float64) * np.nan
        for f in range(ranks):
            file = folder + '/' + str(f) + '.hdf5'

            fin = h5py.File(file,'r')
            if 'lex_data' in fin.keys():
                lx = np.array(fin['lex_'+ident2], dtype='int') // downsample
                ly = np.array(fin['ley_'+ident2], dtype='int') // downsample
                lz = np.array(fin['lez_'+ident2], dtype='int') // downsample
                hx = np.array(fin['hex_'+ident2], dtype='int') // downsample
                hy = np.array(fin['hey_'+ident2], dtype='int') // downsample
                hz = np.array(fin['hez_'+ident2], dtype='int') // downsample

                for b in range(len(lx)):
                    dset = ident+'_'+ident2+'_'+str(b+1)
                    field[lx[b]:hx[b], ly[b]:hy[b], lz[b]:hz[b]] = \
                            (fin[dset][:]).reshape(\
                                    (hx[b]-lx[b], hy[b]-ly[b], hz[b]-lz[b]))
            fin.close()
        return field

    ## Determines what output exists.
    def __ParseFolderStructure(self):
        self.__ParseSlices()
        self.__ParseCoarseBoxes()
        self.__ParseFullBoxes()
        self.__ParseProjections()
        self.__ParseSpectra()
        self.__ParseGravitationalWaveSpectra()
        self.__ParseSlicesTruncationError()
        self.__ParseCoarseBoxesTruncationError()
        self.__ParseFullBoxesTruncationError()

    ## Determines how many slices have been written and when.
    ## Reads header files.
    def __ParseSlices(self):
        folder = self._prefix + '/slices/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._slice_headers.append( fin['Header_x'][:] )
            i = i + 1

        self._slice_headers = np.array( self._slice_headers )

        print('Number of slices found:', len(self._slice_headers))

    ## Determines how many coarse level boxes have been written and when.
    ## Reads header files.
    def __ParseCoarseBoxes(self):
        folder = self._prefix + '/coarse_box/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._coarse_box_headers.append( fin['Header_data'][:] )
            i = i + 1

        self._coarse_box_headers = np.array( self._coarse_box_headers )

        print('Number of coarse boxes found:', len(self._coarse_box_headers))

    ## Determines how many full boxes have been written and when.
    ## Reads header files.
    def __ParseFullBoxes(self):
        folder = self._prefix + '/full_box/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._full_box_headers.append( fin['Header_data'][:] )
            i = i + 1

        self._full_box_headers = np.array( self._full_box_headers )

        print('Number of full boxes found:', len(self._full_box_headers))

    ## Determines how many slices of truncation errors have been written and
    ## when. Reads header files.
    def __ParseSlicesTruncationError(self):
        folder = self._prefix + '/slices_truncation_error/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._slice_truncation_error_headers.append( fin['Header_te_x'][:] )
            i = i + 1

        self._slice_truncation_error_headers =\
                np.array( self._slice_truncation_error_headers )

        print('Number of slices of truncation errors found:',\
                len(self._slice_truncation_error_headers))

    ## Determines how many coarse level boxes have been written and when.
    ## Reads header files.
    def __ParseCoarseBoxesTruncationError(self):
        folder = self._prefix + '/coarse_box_truncation_error/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._coarse_box_truncation_error_headers.append(\
                    fin['Header_data'][:] )
            i = i + 1

        self._coarse_box_truncation_error_headers =\
                np.array( self._coarse_box_truncation_error_headers )

        print('Number of coarse boxes of truncation errors found:',\
              len(self._coarse_box_truncation_error_headers))

    ## Determines how many full boxes have been written and when.
    ## Reads header files.
    def __ParseFullBoxesTruncationError(self):
        folder = self._prefix + '/full_box_truncation_error/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/Level_0/0.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._full_box_truncation_error_headers.append(\
                    fin['Header_data'][:] )
            i = i + 1

        self._full_box_truncation_error_headers =\
                np.array( self._full_box_truncation_error_headers )

        print('Number of full boxes of truncation errors found:',\
              len(self._full_box_truncation_error_headers))

    ## Determines how many projections have been written and parses their header
    ## files.
    def __ParseProjections(self):
        folder = self._prefix + '/projections/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/projections.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._projection_headers.append(fin['Header'][:])
            i = i + 1

        self._projection_headers =\
                np.array( self._projection_headers )

        print('Number of projections found:',\
              len(self._projection_headers))

    ## Determines how many projections have been written and parses their header
    ## files.
    def __ParseSpectra(self):
        folder = self._prefix + '/spectra/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/spectra.hdf5'
            if not path.exists( file ):
                break
            fin = h5py.File(file,'r')
            self._spectrum_headers.append(fin['Header'][:])
            i = i + 1

        self._spectrum_headers =\
                np.array( self._spectrum_headers )

        print('Number of spectra found:',\
              len(self._spectrum_headers))

    ## Determines how many projections have been written and parses their header
    ## files.
    def __ParseGravitationalWaveSpectra(self):
        folder = self._prefix + '/gw_spectra/'
        i = 0

        # iterate over batches and read header of first files
        while True:
            file = folder + str(i) + '/spectra.hdf5'
            if not path.exists( file ):
                break
            try:
                fin = h5py.File(file,'r')
            except:
                print("Warning could not open file", file)
            self._gravitational_wave_spectrum_headers.append(fin['Header'][:])
            i = i + 1

        self._gravitational_wave_spectrum_headers =\
                np.array( self._gravitational_wave_spectrum_headers )

        print('Number of gravitational wave spectra found:',\
              len(self._gravitational_wave_spectrum_headers))


    ## Path of output folder
    _prefix = "."

    ## List of headers files from slices
    _slice_headers = []
    _coarse_box_headers = []
    _full_box_headers = []
    _slice_truncation_error_headers = []
    _coarse_box_truncation_error_headers = []
    _full_box_truncation_error_headers = []
    _projection_headers = []
    _spectrum_headers = []
    _gravitational_wave_spectrum_headers = []
