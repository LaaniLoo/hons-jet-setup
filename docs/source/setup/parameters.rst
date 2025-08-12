.. _Parameters:

Parameters
==========

This page will mainly focus on the jet setup specific parts of ``pluto.ini``. The main reference for the ``pluto.ini`` file is found in the `PLUTO User Guide <http://plutocode.ph.unito.it/userguide.pdf>`_.

Particles
---------

Nparticles
    The first parameter is ignored.
    The second parameter should be set to 1 or higher to enable particle injection.
    If :any:`PARTICLES_INJECT_IN_BATCHES` is no, then the second parameter controls how many particles are injected per PART_INJECT_FREQUENCY.

.. code-block:: none

    Nparticles          -1    1.0

particles_{dbl,flt,vtk,tab}
    Controls the output intervals for the double, float, vtk, or tabulated (ascii) particle data files.
    Same syntax as the ``[Static Grid Output]`` block. Set the first field to specify the time interval (in code units) between consecutive outputs. For example, to output double particle data every 0.01 Myr:

.. code-block:: none

    particles_dbl        3.06601394  -1

User Parameters
---------------

JET_PWR
    Jet power, in :math:`\textrm{erg s}^{-1}` for a **one-sided** jet (default: :math:`1\times10^{45}`)

JET_SPD
    Jet velocity, in units of c (default: 0.9)

JET_CHI
    The value of :math:`\chi` used to calculate the jet density for a relativistic jet when :any:`JET_RELATIVISTIC_DENSITY` is set to ``YES``.
    Otherwise this is taken to be the jet injection region overpressure factor with respect to the environment.
    Default is 100.0

JET_OA_PRIMARY, JET_OA_SECONDARY
    Opening angle of the first and second jet respectively, in degrees. ``JET_OA_PRIMARY`` is used for both jets if :any:`UNIQUE_JETS` is set to ``NO`` (default: :math:`15.0`)

JET_INJECTION_HEIGHT
    Height of jet injection region (for internal boundary injection), in kpc (default: :math:`0.5`)

JET_INITIAL_RADIUS
    Initial radius of jet injection region at the base of the cone, in kpc (default: :math:`0.1`)

JET_X1O, JET_X2O, JET_X3O
    Jet injection region offsets for X1, X2, X3 respectively, in kpc (default: :math:`0.0`)

JET_ROTATION_ANGLE
    Jet injection region rotation angle (in XZ plane) in degrees (default: :math:`0.0`)

JET_START_TIME
    Start time of the jet, in Myr (default: :math:`0.0`)

JET_END_TIME
    End time of the jet, in Myr (default: :math:`310`)

JET_ACTIVE_TIME
    Total active time of the jet, in Myr (default: :math:`310`)

JET_EPISODES
    Number of jet episodes (default: :math:`1.0`)

ENV_RHO_0
    Core density of the environment, in :math:`\textrm{g cm}^{-3}`. Only used for Makino & King profiles (default: :math:`2.4\times10^{-27}`)

R_SCALING
    Scaling radius of the environment, in kpc. Only used for Makino & King profiles (default: :math:`144`)

DELTA_NFW
    NFW parameter, only used for Makino profile (default: :math:`0.0`)

B_EXPONENT
    Beta parameter, only used for King profile (default: :math:`0.38`)

ENV_TEMP
    Environment temperature, in units of Kelvin (default: :math:`3.462\times 10^{7}`)

ENV_X1O, ENV_X2O, ENV_X3O
    Environment centre offsets for X1, X2, X3 respectively, in kpc (default: :math:`0.0`)

PART_INJECT_FREQUENCY
    Particle injection frequency in Myr (default: :math:`0.01`)

STELLAR_MASS
    Stellar mass parameter for supernovae feedback, in units of solar mass. UNTESTED (default: :math:`0.0`)
