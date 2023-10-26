from os import path
from scipy import stats
import numpy as np
import h5py
import matplotlib.pyplot as plt

## Class that handles everything special for the Axion string project
#  such as creating initial states and computing the string length and
#  axion spectra.
class AxionStrings:    
    ## Creates an initial state for Psi1, Psi2, Pi1, and Pi2.
    # @param    L               Box length in units of 1/(a_1 H_1).
    # @param    N               Number of coarse level grid sites.
    # @param    k_max           Maximum wave number to be included.
    # @param    output_file     Name of output file (must include prefix).
    # @param    lambda_param    Coupling parameter lambda.
    def CreateInitialState(self, L, N, k_max, t_start, output_file,\
                           lambda_param = 1., sim_output_with_box_layout=''):
        if sim_output_with_box_layout != '':
            box_layout = self.__GetBoxLayout(sim_output_with_box_layout)
        else:
            box_layout = np.array([])

        # The temperature at eta = 1 in units of f_a.
        T1 = np.sqrt(1.687) 

        # The value of eta at the initial state.
        eta = t_start 

        # The temperature at the initial state.
        T = T1 / eta 

        # Grid spacing.
        dx = L/N 

        # The effective thermal mass.
        mEffSq = lambda_param * (T**2 / 3 - 1) 

        # These possible values of components of the 3-momentum. There are two
        # arrays to take advantage of that our fields are real to reduce the
        # memory overhead.
        k = 2*np.pi * np.fft.fftfreq(N, dx)
        k_real = 2*np.pi * np.fft.rfftfreq(N, dx)

        # Compute the magnitude of the 3-momentum.
        kMags = np.sqrt(k[:, None, None]**2\
                      + k[None, :, None]**2\
                      + k_real[None, None, :]**2)

        # Compute \omega_k and n_k as defined.
        omegaK = np.sqrt(kMags**2 + mEffSq)
        nK = 1 / (np.exp(omegaK/ T) - 1)
        kMax = k_real[k_max]

        # Generate and write files.
        archive = h5py.File(output_file, 'w')
        self.__SaveField(archive, 'Psi1', box_layout, True, N,\
                         omegaK, kMags, nK, k_max, eta)
        self.__SaveField(archive, 'Psi2', box_layout, True, N,\
                         omegaK, kMags, nK, k_max, eta)
        self.__SaveField(archive, 'Pi1', box_layout, True, N,\
                         omegaK, kMags, nK, k_max, eta)
        self.__SaveField(archive, 'Pi2', box_layout, True, N,\
                         omegaK, kMags, nK, k_max, eta)
        archive.close()
        print('Done.')

    def PlotAxionAndRadialModeSlice(self, axion, radial_mode, save_name=""):
        fig, ax = plt.subplots(figsize=(20,10),ncols=2)

        im = ax[0].imshow(axion, cmap='twilight', vmin=-np.pi, vmax=np.pi,\
                          interpolation="nearest")
        cb = fig.colorbar(im, ax=ax[0])
        cb.ax.tick_params(labelsize=15)
        cb.set_label(label='axion',fontsize=20)
        ax[0].set_xticks([])
        ax[0].set_yticks([])

        im = ax[1].imshow(radial_mode, cmap='gist_heat_r')
        cb = fig.colorbar(im, ax=ax[1])
        cb.ax.tick_params(labelsize=15)
        cb.set_label(label='radial mode',fontsize=20)
        ax[1].set_xticks([])
        ax[1].set_yticks([])

        plt.tight_layout()
        if save_name != "":
            plt.save_fig(save_name)

    ## Returns an indiviual field.
    # @param    field   Boolean value. True: Return Psi, False: Return Pi.
    # @param    N       Number of coarse level grid sites.
    # @param    omegaK  \omega_k.
    # @param    kMags   Magnitude of 3-momentum.
    # @param    nK      n_k.
    # @param    k_max   Maximum wave number to be included.
    # @param    eta     Starting eta of simulation.
    def __GetField(self, field, N, omegaK, kMags, nK, k_max, eta):
        # Compute the norm and phase.
        norm = (N / 2 / np.pi)**3
        phase = np.random.uniform(size = kMags.shape)*2*np.pi*1j

        # Compute magnitude depending of if we want Psi or Pi.
        magnitude = nK * stats.chi2.rvs(df = 2, size = omegaK.shape) 
        if field:
            magnitude = np.nan_to_num(magnitude  / omegaK)
        else:
            magnitude = np.nan_to_num(magnitude  * omegaK)

        magnitude[kMags > k_max] = 0
        magnitude[kMags == 0] = 0

        # Generate spectra.
        field_spectrum = np.sqrt(magnitude)*np.exp(phase)*norm
   
        # Return inverse FFT of spectra.
        if field:
            return np.fft.irfftn(field_spectrum)
        else:
            return eta * np.fft.irfftn(field_spectrum)

    def __SaveField(self, file, name, box_layout, field, N, omegaK, kMags, nK,
                    k_max, eta):
        print('Generate '+name+' ...')
        field_state = self.__GetField(field, N, omegaK, kMags, nK, k_max, eta)

        if len(box_layout) == 0:
            file.create_dataset(name, data = field_state)
        else:
            for i in range(len(box_layout[0])):
                file.create_dataset(name+'_'+str(i),
                        data = field_state[box_layout[0][i]:box_layout[3][i]+1,
                                           box_layout[1][i]:box_layout[4][i]+1,
                                           box_layout[2][i]:box_layout[5][i]+1])

        del field_state

    def __GetBoxLayout(self, sim_output_with_box_layout):
        filename = sim_output_with_box_layout + '/box_layout.h5'
        file = h5py.File(filename,'r')
        x0 = np.array(file['x0'][:], dtype='int')
        y0 = np.array(file['y0'][:], dtype='int')
        z0 = np.array(file['z0'][:], dtype='int')
        x1 = np.array(file['x1'][:], dtype='int')
        y1 = np.array(file['y1'][:], dtype='int')
        z1 = np.array(file['z1'][:], dtype='int')
        file.close()
        return np.array([x0,y0,z0,x1,y1,z1])

