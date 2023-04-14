from os import path
from scipy import stats
import numpy as np
import h5py

## Class that handles everything special for the Axion string project
#  such as creating initial states and computing the string length and
#  axion spectra.
class AxionStrings:
	
	## Creates an initial state for Psi1, Psi2, Pi1, and Pi2.
	# @param	L		Box length in units of 1/(a_1 H_1).
	# @param	N		Number of coarse level grid sites.
	# @param	k_max		Maximum wave number to be included.
	# @param	output_file	Name of output file (must include prefix).
	# @param	lambda_param	Coupling parameter lambda.
	def CreateInitialState(self, L, N, k_max, t_start, output_file, lambda_param = 1.):
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

		# These possible values of components of the 3-momentum. There are two arrays
		# to take advantage of that our fields are real to reduce the memory overhead.
		k = 2*np.pi * np.fft.fftfreq(N, dx)
		k_real = 2*np.pi * np.fft.rfftfreq(N, dx)

		# Compute the magnitude of the 3-momentum.
		kMags = np.sqrt(k[:, None, None]**2 + k[None, :, None]**2 + k_real[None, None, :]**2)

		# Compute \omega_k and n_k as defined.
		omegaK = np.sqrt(kMags**2 + mEffSq)
		nK = 1 / (np.exp(omegaK/ T) - 1)
		kMax = k_real[k_max]

		# Generate and write files.
		archive = h5py.File(output_file, 'w')

		print('Generate Psi1...')
		Psi1 = self.__GetField(True, N, omegaK, kMags, nK, k_max, eta)
		archive.create_dataset('Psi1', data = Psi1)
		del Psi1

		print('Generate Psi2...')
		Psi2 = self.__GetField(True, N, omegaK, kMags, nK, k_max, eta)
		archive.create_dataset('Psi2', data = Psi2)
		del Psi2

		print('Generate Pi1...')
		Pi1 = self.__GetField(False, N, omegaK, kMags, nK, k_max, eta)
		archive.create_dataset('Pi1', data = Pi1)
		del Pi1

		print('Generate Pi2...')
		Pi2 = self.__GetField(False, N, omegaK, kMags, nK, k_max, eta)
		archive.create_dataset('Pi2', data = Pi2)
		del Pi2

		archive.close()	
		print('Done.')

	## Returns an indiviual field.
	# @param	field	Boolean value. True if Psi is to be return, false for Pi.
	# @param	N	Number of coarse level grid sites.
	# @param	omegaK	\omega_k.
	# @param	kMags	Magnitude of 3-momentum.
	# @param	nK	n_k.
	# @param	k_max	Maximum wave number to be included.
	# @param	eta	Starting eta of simulation.
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
