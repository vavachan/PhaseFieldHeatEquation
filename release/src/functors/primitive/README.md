This folder contains _primitive_ functors that implement common computations and patterns usable in a variety of situtations.

Examples for such functors include:

- LBM primitives such as functions to compute equilibria, recover moments and so on
- (inter,extra)polation schemes
- Basic ODE solvers
- Geometry primitives such as signed distance functions

Such functors can be simple free (template) functions or non-virtual classes.

All functionality contained in this folder must be declared `any_platform` and be executable both on CPUs and GPUs.
