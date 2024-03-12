# Sledgehamr: Simulating scalar fields with adaptive mesh refinement
Sledgehamr (**S**ca**L**ar fi**E**ld **D**ynamics **G**etting solv**E**d wit**H** **A**daptive **M**esh **R**efinement) is an AMReX-based code to simulate the dynamics of coupled scalar fields on a 3-dimensional mesh. Adaptive mesh refinement (AMR) can boost performance if spatially localized regions of the scalar field require high resolution. Compatible with both GPU and CPU clusters, sledgehamr offers a flexible and customizable framework. This framework enables various applications, such as the generation of gravitational wave spectra.

<p align="left">
  <img width="512" height="512" src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/axion.gif">
</p>

For a detailed description of the code please consult the accompanying paper:
https://arxiv.org/abs/23xx.xxxxx

For questions, comments, or bug reports please use the GitHub issues feature, or contact the author:
Malte Buschmann (m.s.a.buschmann@uva.nl)

## Installation

### Prerequisites
* AMReX (```git clone https://github.com/AMReX-Codes/amrex```)
* FFTW3
* Boost

### Create a project
Sledgehamr comes with a few different physics scenarios already implemented such as axion strings and a first order phase transition of a single scalar field, as well as two example projects ```MinimalExample``` and ```NextToMinimalExample```. The accompanying paper describes in detail how other scenarios can be implemented.

### Running the minimal examples
The minimal examples can be run by following the instructions in the Jupyter notebooks ```notebooks/MinimalExample.ipynb``` and ```notebooks/NextToMinimalExample.ipynb```.

## Code documentation
* https://msabuschmann.github.io/sledgehamr/
* https://arxiv.org/abs/23xx.xxxxx 

## How to cite
If you use sledgehamr, please cite its accompanying paper:

* Malte Buschmann, "Sledgehamr: Simulation scalar fields with adaptive mesh refinement",
arXiv:23xx.xxxxx [hep-ph].

BibTex:

## Publications using this code

### 2022
* "Dark matter from axion strings with adaptive mesh refinement", M. Buschmann et al., Nature Commun. 13 (2022) 1, 1049, https://arxiv.org/abs/2108.05368

### 2023
* "The Cosmological Dynamics of String Theory Axion Strings", J. Benabou et al., https://arxiv.org/abs/2312.08425

### 2024
* "Signatures of primordial energy injection from axion strings", J. Benabou et al., Phys.Rev.D 109 (2024) 5, 055005, https://arxiv.org/abs/2308.01334
* "Thick and Thin Wall Collisions with Adaptive Mesh Refinement", M. Buschmann et al., to appear
* "Axion Mass Prediction from adaptive mesh refinement cosmological lattice simulations", J. Benabou et al., to appear


