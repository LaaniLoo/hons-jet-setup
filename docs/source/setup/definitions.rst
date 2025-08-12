.. _Configuration:

Configuration
=============

The PLUTO code is configured through definitions.
These control which features are enabled, and alter how certain features behave.
As this code is purpose-built for AGN jet simulations, several new configuration options have been added, in addition to the definitions found in ``definitions.h`` in the original PLUTO distribution.
These reside in the ``definitions_usr.h`` file, while the option definitions can be found in ``pluto_usr.h``.

The main reference for the ``definitions.h`` file is found in the `PLUTO User Guide <http://plutocode.ph.unito.it/userguide.pdf>`_, while the documentation for the new definitions is below.

.. doxygenfile:: definitions_usr.h
