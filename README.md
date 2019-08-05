		ECMech

	8888888888 .d8888b.  888b     d888                   888
	888       d88P  Y88b 8888b   d8888                   888
	888       888    888 88888b.d88888                   888
	8888888   888        888Y88888P888  .d88b.   .d8888b 88888b.
	888       888        888 Y888P 888 d8P  Y8b d88P"    888 "88b 
	888       888    888 888  Y8P  888 88888888 888      888  888 
	888       Y88b  d88P 888   "   888 Y8b.     Y88b.    888  888 
	8888888888 "Y8888P"  888       888  "Y8888   "Y8888P 888  888 
                                                              
BACKGROUND
======

ECMech is a GPU-friendly library of constitutive models. The models are based on standard continuum mechanics concepts. Crystal-mechanics-based and porosity-mechanics-based models are a principal focus. Models are meant for standard crystalline metallic materials deforming under quasi-static conditions. 

Constitutive model response is a main ingredient in the simulation of deformation of material. In the solution of the governing equations, the constitutive model sub-problem is typically solved at every discretization point in the problem. In a finite element setting for solving the balance of momentum (stress equilibrium in the quasi-static case), these would be quadrature points in the finite element integrals. The constitutive model has two main jobs: it provides the stress tensor that goes into the balance of linear momentum, and it updates the evolving state of the material. This state can be tracked by variables for the stress, dislocation density, orientation of crystal lattices, grain size, and so forth.

Implementations here are meant for quasi-static applications, and include material tangent stiffness computations for incorporation into finite element stiffness matrices for quasi-static applications.

A common feature of advanced constitutive models is the need to solve coupled systems of non-linear equations. These solutions typically occur at each time step at every integration point in the host code, but with computational work at the various integration points being independent and well-suited to parallel computation. These multi-dimensional systems of equations arise in a variety of models including those for crystal-mechanics, porous-plasticity, and phase kinetics. The systems of equations are of comparatively small size; but they are often numerically stiff, requiring robust iterative solution schemes to produce numerically stable results. ECMech makes use of the open-source library snls to solve these systems of equations.

Examples of work that has made use of similar models and algorithms:
  * [Journal publication](https://dx.doi.org/10.1088/1361-651X/aa841c) on application to an alpha-beta titanium alloy
  * [Journal publication](https://dx.doi.org/10.1016/j.actamat.2017.02.042) for comparison to diffraction data

BUILDING
======

The build system is cmake-based.

Dependencies:
* blt -- required
  - https://github.com/LLNL/blt
  - use `git submodule` commands, or specify with `-DBLT_SOURCE_DIR=...`
* snls -- required (currently, needs to be the `templated` branch)
  - https://github.com/LLNL/SNLS
  - use `git submodule` commands, or specify with `-DSNLS_DIR=...`
* raja -- required
  - https://github.com/LLNL/RAJA
  - specify with `-DRAJA_DIR=...`

TESTING
======

Run cmake with `-DENABLE_TESTS=ON` and do `make test`

DEVELOPMENT
======

The develop branch is the main development branch. Changes to develop are by pull request.

AUTHORS
======

The principal developer of ECMech is Nathan Barton, nrbarton@llnl.gov. Robert Carson has also made contributions. 

LICENSE
======

License is under the BSD-3-Clause license. See [LICENSE](LICENSE) file for details. And see also the [NOTICE](NOTICE) file. 

`SPDX-License-Identifier: BSD-3-Clause`

``LLNL-CODE-784997``
