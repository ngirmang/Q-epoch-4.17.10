# Quantum EPOCH

The EPOCH code with quantum add-on modules for modelling condensed matter electromagnetic response. See the paper for details.

Compilation is the same as the usual code, please see the old Readme.

## Relevant Defines

As seen in the `Makefile`, the following `#define`'s have been added that allow turning on new features and deck options:

* `BOUND_HARMONIC` adds a linear dispersive response via a particle species. This adds the `harmonic_omega` and `harmonic_gamma` options to the species block.

* `CONSTEPS` adds a constant dielectric constant in the `eps` option in the fields block.

* `NEWPML` adds an absorbing layer that might be more stable than the CPML.

A number of others will be documented elsewhere.


