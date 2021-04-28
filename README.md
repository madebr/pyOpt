# pyOpt
PYthon OPTimization Framework
Copyright (c) 2008-2014, pyOpt Developers

pyOpt is an object-oriented framework for formulating and solving
nonlinear constrained optimization problems.

Some of the features of pyOpt:

*   Object-oriented development maintains independence between
    the optimization problem formulation and its solution by
    different optimizers
*   Allows for easy integration of gradient-based, gradient-free,
    and population-based optimization algorithms
*   Interfaces both open source as well as industrial optimizers
*   Ease the work required to do nested optimization and provides
    automated solution refinement
*   On parallel systems it enables the use of optimizers when
    running in a mpi parallel environment, allows for evaluation
    of gradients in parallel, and can distribute function
    evaluations for gradient-free optimizers
*   Optimization solution histories can be stored during the
    optimization process. A partial history can also be used
    to warm-restart the optimization

see QUICKGUIDE.md for further details.

## Building and installing

### Requirements:

- python
- numpy and numpy-ext
- fortran compiler
- swig


### Build commands

Build default pyOpt

```sh
python setup.py build_ext --inplace
```

Build debug pyOpt with no optimization

```sh
python setup.py config_fc --debug --noopt build_ext --inplace
```

Get information about the available compilers

```sh
python setup.py config_fc --help-fcompiler
```


## Licensing
Distributed using the GNU Lesser General Public License (LGPL); see
the LICENSE file for details.

Please cite pyOpt and the authors of the respective optimization
algorithms in any publication for which you find it useful.
(This is not a legal requirement, just a polite request.)


## Contact and Feedback

If you have questions, comments, problems, want to contribute to the
framework development, or want to report a bug, please contact the
main developers:

-   [Dr. Ruben E. Perez](mailto:Ruben.Perez@rmc.ca)
-   [Peter W. Jansen](mailto:Peter.Jansen@rmc.ca)
