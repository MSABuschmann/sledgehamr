from os import path
import numpy as np
import h5py

## This class reads all SledgeHAMR output
class Output:
    ## Crawls the output and loads headers.
    def __init__(self, output_folder):
        self._prefix = output_folder
        self._slice_headers = []
        self._coarse_box_headers = []
        self._full_box_headers = []

        self.__ParseFolderStructure()

    ## Returns array of times at which slices have been written.
    def GetTimesOfSlices(self):
        return self._slice_headers[:,0]

    ## Returns array of times at which coarse boxes have been written.
    def GetTimesOfCoarseBoxes(self):
        return self._coarse_box_headers[:,0]

    ## Returns array of times at which full boxes have been written.
    def GetTimesOfFullBoxes(self):
        return self._full_box_headers[:,0]

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
        folder = self._prefix + '/slices/'

        # Start dictionary
        d = dict();
        d['t'] = t

        # Loop over fields, files, and boxes to reassemble data
        for s in fields:
            field = np.ones((dim,dim), dtype=np.float32) * np.nan
            for f in range(ranks):
                file = folder + str(i) + '/Level_'+str(level)\
                     + '/' + str(f) + '.hdf5'

                fin = h5py.File(file,'r')
                if 'le1_'+direction in fin.keys():
                    le1 = np.array(fin['le1_'+direction], dtype='int')
                    le2 = np.array(fin['le2_'+direction], dtype='int')
                    he1 = np.array(fin['he1_'+direction], dtype='int')
                    he2 = np.array(fin['he2_'+direction], dtype='int')

                    for b in range(len(le1)):
                        dset = s+'_'+direction+'_'+str(b+1)
                        field[le1[b]:he1[b],le2[b]:he2[b]] =\
                                (fin[dset][:]).reshape(\
                                        (he1[b]-le1[b], he2[b]-le2[b]))
                fin.close()

            # Add result to dict.
            d[s] = field

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
        folder = self._prefix + '/coarse_box/' + str(i)

        # Start dictionary
        d = dict();
        d['t'] = t

        for s in fields:
            d[s] = self.__Read3dField(folder, dim, ranks, s, downsample)

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
            d[s] = self.__Read3dField(folder, dim, ranks, s, downsample)

        return d

    def __Read3dField(self, folder, dim, ranks, ident, downsample):
        field = np.ones((dim, dim, dim), dtype=np.float32) * np.nan
        for f in range(ranks):
            file = folder + '/' + str(f) + '.hdf5'

            fin = h5py.File(file,'r')
            if 'lex_data' in fin.keys():
                lx = np.array(fin['lex_data'], dtype='int') // downsample
                ly = np.array(fin['ley_data'], dtype='int') // downsample
                lz = np.array(fin['lez_data'], dtype='int') // downsample
                hx = np.array(fin['hex_data'], dtype='int') // downsample
                hy = np.array(fin['hey_data'], dtype='int') // downsample
                hz = np.array(fin['hez_data'], dtype='int') // downsample

                for b in range(len(lx)):
                    dset = ident+'_data_'+str(b+1)
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
            file = folder + str(i) + '/0.hdf5'
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

    ## Path of output folder
    _prefix = "."

    ## List of headers files from slices
    _slice_headers = []
    _coarse_box_headers = []
    _full_box_headers = []
