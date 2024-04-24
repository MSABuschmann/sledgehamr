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

## Compute the effective momentum k from an index for different projection types.
# @param    a               Index.
# @param    N               Maximum index.
# @param    projection_type Projection type.
def IndexToK(a, N, projection_type):
    n_tilde = a
    n_tilde[np.where(a - N > -N / 2 - 1)] -= N

    two_pi_n_tilde = 2. * np.pi / N * n_tilde
    if projection_type == 1:
        return np.sin(two_pi_n_tilde)
    elif projection_type == 2:
        return (8. * np.sin(two_pi_n_tilde) - np.sin(2. * two_pi_n_tilde)) / 6.
    return 0.

## Add unique k values to file.
# @param    data_dir    Sledgehamr data directory.
# @param    dim         Dimension for which to generate k values.
def CreateSpectrumBins(data_dir, dim):
    file = 'spectra_ks.hdf5'
    filename = data_dir +'/' + file

    archive = h5py.File(filename, 'a')
    if str(dim) not in archive.keys():
        vals = CreateUniqueVals(dim)
        i = np.linspace(0, dim-1, dim)
        archive.create_dataset(str(dim)+'_bins', data = vals)
        archive.create_dataset(str(dim)+'_nks', data = np.array([len(vals)]))
        archive.create_dataset(str(dim)+'_k1', data = IndexToK(np.copy(i), dim, 1))
        archive.create_dataset(str(dim)+'_k2', data = IndexToK(np.copy(i), dim, 2))
        archive.close()
        print('Binning for a '+str(dim)+'^3 spectrum added!')
    else:
        print('Binning already exists in file.')
