Unit values
===========

The unit system used for the simulations is defined in ``definitions.h`` with the ``UNIT_DENSITY``, ``UNIT_LENGTH``, and ``UNIT_VELOCITY`` definitions.

By default, the following unit system is defined:

.. code-block:: c

    #define  MU_NORM                        0.60364
    #define  UNIT_DENSITY                   ((CONST_amu)*MU_NORM)
    #define  UNIT_LENGTH                    ((CONST_pc)*1e3)
    #define  UNIT_VELOCITY                  CONST_c

where ``CONST_*`` are the constant values defined in cgs.

This produces the following unit values:

Length
    :math:`1\,\textrm{kpc}`

Density
    :math:`0.60364 * 1\,\textrm{amu} = 1.0023678\times10^{-24}\,\textrm{g cm}^{-3}`

Speed
    :math:`1\,\textrm{c}`, speed of light

Time
    :math:`0.0032615638\,\textrm{Myr}`

Pressure
    :math:`9.0088324\times10^{-5}\,\textrm{Pa}`

Mass
    :math:`2.9449555\times10^{37}\,\textrm{kg}`

Energy
    :math:`2.646794\times10^{54}\,\textrm{J}`

General notes
-------------

The grid setup (in ``pluto.ini``) is defined in length units.
This means that for the default length unit of 1 kpc, the grid is defined in kpc.

The data output, particle output, and simulation times are in simulation time units.
