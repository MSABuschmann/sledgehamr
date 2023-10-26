import numpy as np
import h5py

def CreateUniqueVals(dim):
    kVals = np.array(np.fft.fftfreq(dim)*dim, dtype = 'int')
    unique_mags = np.unique(np.unique(np.unique(kVals**2)[:, None]\
                + np.unique(kVals**2)[None, :])[:, None]\
                + np.unique(kVals[None, :]**2)[None, :])
    return unique_mags

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
