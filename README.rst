PINC
====

*Particle-In-Cell* (PINC) is an open source scientific program for simulating plasmas using the *Particle-In-Cell* (PIC) method on a structured mesh. The field quantities are solved using the *Geometric Multigrid Method*. The focus is on plasma-object interaction, collisions in Farley-Buneman instabilities, and blob instability in tokamaks.

A user guide is unfortunately unavailable.

Contributors
------------

Authors:

- `Sigvald Marholm`_: architecture and data structures, pushers, weighting schemes, spectral solver, parallelization
- Gullik Killie: architecture and data structures, multigrid solver, parallelization
- Vigdis Holta (see separate branch): Neumann boundaries, blob instability for tokamak simulations
- Steffen Brask (see separate branch): Collision module, Farley-Buneman instability
- Jan Deca (see separate branches): Capacitance matrix method for plasma-object interactions

.. _`Sigvald Marholm`: mailto:sigvald@marebakken.com
