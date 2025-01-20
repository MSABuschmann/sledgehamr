[![arXiv](https://img.shields.io/badge/arXiv-2404.02950%20-green.svg)](https://arxiv.org/abs/2404.02950)

# Sledgehamr: Simulating scalar fields with adaptive mesh refinement
Sledgehamr (**S**ca**L**ar fi**E**ld **D**ynamics **G**etting solv**E**d wit**H** **A**daptive **M**esh **R**efinement) is an AMReX-based code to simulate the dynamics of coupled scalar fields on a 3-dimensional mesh. Adaptive mesh refinement (AMR) can boost performance if spatially localized regions of the scalar field require high resolution. Compatible with both GPU and CPU clusters, sledgehamr offers a flexible and customizable framework. This framework enables various applications, such as the generation of gravitational wave spectra.

<p align="left">
  <img width="512" height="512" src="https://github.com/MSABuschmann/sledgehamr/blob/main/assets/axion.gif">
</p>

For a detailed description of the code please consult the accompanying paper:
https://arxiv.org/abs/arXiv:2404.02950

For questions, comments, or bug reports please use the GitHub issues feature, or contact me via email:
m.s.a.buschmann@uva.nl

## Installation

### Prerequisites
* AMReX (```git clone https://github.com/AMReX-Codes/amrex```)
* FFTW3

### Running the minimal examples
The minimal examples can be run by following the instructions in the Jupyter notebooks ```notebooks/MinimalExample.ipynb``` and ```notebooks/NextToMinimalExample.ipynb```.

### Create a project
Sledgehamr comes with a few different physics scenarios already implemented such as axion strings and a first-order phase transition of a single scalar field + gravitational waves, as well as two example projects ```MinimalExample``` and ```NextToMinimalExample```. The accompanying paper describes in detail how other scenarios can be implemented.

## Code documentation
* https://msabuschmann.github.io/sledgehamr/
* https://arxiv.org/abs/2404.02950

## How to cite
If you use sledgehamr, please cite its accompanying paper:

* Malte Buschmann, "Sledgehamr: Simulation Scalar Fields with Adaptive Mesh Refinement",
arXiv:2404.02950

BibTex:
```
@misc{buschmann2024sledgehamr,
      title={Sledgehamr: Simulating Scalar Fields with Adaptive Mesh Refinement}, 
      author={Malte Buschmann},
      year={2024},
      eprint={2404.02950},
      archivePrefix={arXiv},
      primaryClass={hep-ph}
}
```

## Publications using this code

### 2022
* "Dark matter from axion strings with adaptive mesh refinement", M. Buschmann et al., Nature Commun. 13 (2022) 1, 1049, https://arxiv.org/abs/2108.05368

### 2024
* "Signatures of primordial energy injection from axion strings", J. Benabou et al., Phys.Rev.D 109 (2024) 5, 055005, https://arxiv.org/abs/2308.01334
* "Axion mass prediction from adaptive mesh refinement cosmological lattice simulations", J. Benabou et al., https://arxiv.org/abs/2412.08699

### 2025
* "Thick and Thin Wall Collisions with Adaptive Mesh Refinement", M. Buschmann et al., to appear


