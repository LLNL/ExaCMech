# JAX ECMech

This is a port of the C++ ECMech code over to the python Google JAX framework. The purpose of this is to provide users with a framework in which they can quickly iterate on new model forms for crystal plasticity models without having to worry about all the fun numerics of getting things right the first time. In other words, users don't have to worry about getting the derivative terms right during the initial exploration stage of things. So, model forms should be quicker to iterate on as we use JAX's automatic differentiation capabilities for our derivative terms. The added bonus of this is that we can now easily check that our later hand-derived or potentially AD derived C++ derivatives are at least consistent with what JAX is providing.

In order to provide an example of how to use the code, we have provided a replicant of the `pyecmech/example.py` in this new framework in the `jax_ecmech/example.py` file. This simple example takes a copper-like material and uniaxially deforms it. It should be noted that the answers from these calculations are more akin to a material point simulation (MPS) and do not accound for how a material's body's deformation might affect the applied velocity gradient / deformation gradient. So, the answer won't necessarily be consistent with those from a finite element (FE) code like ExaConstit. However, one can still get an idea for how the material might behave based on these simple MPS runs.

This project is very much a work in-progress as it's a new addition to the code-base. Some basic testing has been conducted so far to ensure that the results between the C++ and JAX version of things match, but mistakes can happen. So, if we messed up translating things from the C++ code please let us know. Alternatively, if you have improvements on anything done in here please feel to open up a PR for it, and we'll take a look to see if it should be brought back in.

Nonetheless, we hope you enjoy this new tool and if you find it useful in developing a new model please give us a citation as it helps us in a number of different ways including networking with our community. 

## Required Software

* [Google's JAX software](https://github.com/google/jax/) - see their install notes but this is pip installable which is a plus
* [Optimistix - nonlinear solver](https://github.com/patrick-kidger/optimistix) - see their install notes but this can also be installed using pip
* [ExaCMech](https://github.com/LLNL/ExaCMech) - You'll need to install ExaCMech and it's python bindings as we make use of some of the constants provided by ECMech. Note, we are working on providing a pure python version of things if you don't want to deal with installing ExaCMech and just want to run the python code.


## Programming notes

So, JAX is quite fussy about what it will allow as valid code if we want things to just-in-time (JIT) compile. JIT is an important feature to strive for as it typically leads to much faster simulation runs if the calculations are expensive, or if they're running for a long period of time. The most notable thing one has to avoid for this aspect of things is you can't use `if-else` statements in your code. Instead, you need to use the JAX logic control flow that they provide. This might be things like their version of `np.where` or one of the several tools in their `jax.lax.*` module.

### Note: citation and license are the same as ECMech as it's a part of the ECMech library... below are just for reference though

## CITATION

ECMech may be cited using the following `bibtex` entry:

```
@Misc{ecmech,
title = {{ECMech}},
author = {Barton, Nathan R. and Carson, Robert A. and Wopschall, Steven R and USDOE National Nuclear Security Administration},
abstractNote = {{ECMech} is a {GPU}-friendly library of constitutive models. Crystal-mechanics-based and porosity mechanics-based models are a principal focus.},
url = {https://github.com/LLNL/ExaCMech},
doi = {10.11578/dc.20190809.2},
year = {2018},
month = {12},
annote = {
   https://www.osti.gov//servlets/purl/1550790
   https://www.osti.gov/biblio/1550790-llnl-exacmech
}
}
```

## LICENSE

License is under the BSD-3-Clause license. See [LICENSE](LICENSE) file for details. And see also the [NOTICE](NOTICE) file. 

`SPDX-License-Identifier: BSD-3-Clause`

``LLNL-CODE-784997``