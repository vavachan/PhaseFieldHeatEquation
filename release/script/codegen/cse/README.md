# OpenLB Automatic Code Generation

Common subexpression elimination (CSE) is an established approach to reducing
the arithmetic complexity of computation kernels in LBM. OpenLB currently
provides CSE-optimized versions of common LBM collision operators.

## Prerequisites:

- python 
  - Mako
  - sympy
- gcc
- make

(When using Nix Flake: instantiating evironment for cse via 'nix develop .#env-code-generation')

## Execute in order

[1]: make prepare_extraction
[2]: make extraction
[3]: make optimize_cse

## Directory tree

- root: script/codegen/cse/
  - source/            # contains python source files
  - extractions/       # performing expression extraction
  - expressions/       # cse optimizing the expression tree
  - expressions/cse/   # contains the optimized collision operators
