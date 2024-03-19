# Next-to-Minimal Example 

This example is based on the same scenario as the Minimal Example, however, uses some of the features sledgehamr offers:
* We introduce an extra run-time parameter $\lambda$ which is loaded from the inputs file and used in the equations of motion.
* We output a projection of the axion energy density $\dot{a}^2$.
* We output the power-spectrum of $\dot{a}^2$.
* We implement a custom output type where we save the time-evolution of the average value of the radial mode $r$.

## How to run
1.  Ensure prerequisites are installed: Boost, FFTW3, and HDF5.
2.  Clone the AMReX repository (https://github.com/AMReX-Codes/amrex)
3.  Make sure the paths to the sledgehamr repository (SLEDGEHAMR_HOME) and AMReX
    repository (AMREX_HOME) are set in 'Makefile'. Adjust other compiler flags
    if needed.
4.  Compile sledgehamr: make -j 6
5.  Use the Jupyter Notebook notebooks/MinimalExample.ipynb to generate an
    initial state if it hasn't already been generated for the Minimal Example.
6.  Run the executable. An example slurm submission script is provided: run.sh
5.  Use the example notebook notebooks/NextToMinimalExample.ipynb to plot the results

## Simulation output
Below are illustrations of the simulation output:

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/next_to_minimal_example_projection.png">
</p>

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/next_to_minimal_example_spectrum.png">
</p>

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/next_to_minimal_example_vev.png">
</p>
