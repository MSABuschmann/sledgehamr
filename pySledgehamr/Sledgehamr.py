import numpy as np
import h5py

## Generates the unique k values (combinations of i^2 + j^2 + k^2).
# @param    dim Maximum value for i-1, j-1, k-1.
# @return   Unique values.
def CreateUniqueVals(dim):
    kVals = np.array(np.fft.fftfreq(dim)*dim, dtype = 'int')
    unique_mags = np.unique(np.unique(np.unique(kVals**2)[:, None]\
                + np.unique(kVals**2)[None, :])[:, None]\
                + np.unique(kVals[None, :]**2)[None, :])
    return unique_mags

## Add unique k values to file.
# @param    data_dir    Sledgehamr data directory.
# @param    dim         Dimension for which to generate k values.
def CreateSpectrumBins(data_dir, dim):
    file = 'spectra_ks.hdf5'
    filename = data_dir +'/' + file

    archive = h5py.File(filename, 'a')
    if str(dim) not in archive.keys():
        vals = CreateUniqueVals(dim)
        archive.create_dataset(str(dim), data = vals)
        archive.create_dataset(str(dim)+'_nks', data = np.array([len(vals)]))
        archive.close()
        print('Binning for a '+str(dim)+'^3 spectrum added!')
    else:
        print('Binning already exists in file.')
