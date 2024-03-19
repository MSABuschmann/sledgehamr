# Minimal Example 

This minimal example illustrates how to run sledgehamr. We use axion strings as an example, which are based on a complex scalar field $\Phi$ with dimensionless components $(\psi_1, \psi_2)^T$.

The equations of motion are 
* $\psi_i' = \Pi_i$
  
and

* $\Pi_i' = -\frac{2}{\eta}\Pi_i + \nabla^2\psi_i - \eta^2\psi_i(\psi_1^2+\psi_2^2-1) - a\psi_i$,
  
where $\Pi_i$ is the conjugate momenta of $\psi_i$, $\eta$ is time, and $a$ some constant.

## How to run
1.  Ensure prerequisites are installed: Boost, FFTW3, and HDF5.
2.  Clone the AMReX repository (https://github.com/AMReX-Codes/amrex)
3.  Make sure the paths to the sledgehamr repository (SLEDGEHAMR_HOME) and AMReX
    repository (AMREX_HOME) are set in 'Makefile'. Adjust other compiler flags
    if needed.
4.  Compile sledgehamr: make -j 6
5.  Use the Jupyter Notebook notebooks/MinimalExample.ipynb to generate an
    initial state.
6.  Run the executable. An example slurm submission script is provided: run.sh
5.  Use the example notebook notebooks/MinimalExample.ipynb to plot the results

## Simulation output
Below are illustrations of the simulation output for the two scalar degrees of freedom, the axion $a=\mathrm{arctan2}(\psi2, \psi1)$ and the radial mode $r=\sqrt{\psi_1^2+\psi_2^2}$.

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/minimal_example_slice_level_0.png">
</p>

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/minimal_example_slice_level_1.png">
</p>

<p align="left">
  <img src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/minimal_example_box.png">
</p>
