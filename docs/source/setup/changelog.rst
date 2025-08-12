Changelog
=========

Unreleased / in development
---------------------------

* Add the :any:`EXT_ENV_VELOCITY_FACTOR` definition.
  This is used to scale external environment velocities if :any:`EXT_ENV_VELOCITY` is ``YES``.
* The :any:`JET_RELATIVISTIC_DENSITY` definition is introduced, to allow a choice between relativistic or classical jet power relation.
  This is irrespective of the actual physics module used.
* Added documentation using doxygen, sphinx, breathe, sphinx-multiversion, GitLab CI, and the build-the-docs theme.
* New spherical injection region option
* Improve conical injection region smoothing
* Fix bugs in conical injection density calculations

  * Account for jet injection cone offset
  * Remove clamping of injected density w.r.t. ambient density

* Restart with particles even if time is too different
* Improve handling of lower radial boundary in spherical coordinates
* Add Makino gravitational acceleration

  * Add test for Makino gravitational acceleration

* Treat ``JET_SPD`` parameter as lorentz factor if greater than speed of light (1.0)
* Reset integration timestep at jet injection
* Only inject particles if the jet is active
* Add support for a moving injection region
* Refactor particles and jet injection
* Add support for injecting particles as a stream, rather than in batches
* Refactor jet intermittency code, clean up ``init.c``

6.3
---

**Minor changes:**

* Add new particle shock parameters (introduced in 6.1) to the default definitions file.

6.2
---

**Incompatiblities:**

.. warning::
    The particle restart functionality has changed in this release.
    Restarting with outputs from a previous version will not work.
 
**New features:**

* Grid tracer values are now saved in particle outputs, alongside density and pressure information.

**Fixes:**

* Fixed an crash when restarting from a simulation that didn't include particle outputs.

6.1
---

**New features:**

* Support for multiple particle shock thresholds.
  They are controlled through the following new configuration definitions:

``PARTICLES_LP_SHK_BINS``
    Number of thresholds (default 3)

``PARTICLES_LP_SHK_THRESH_MIN``
    Minimum shock threshold (default 0.05)

``PARTICLES_LP_SHK_THRESH_MAX``
    Maximum shock threshold (default 5.0)

**Fixes:**

* Particle restart functionality has been improved.
* Particle fields with multiple dimensions are now handled correctly.

6.0
---

**New features:**

* Added support for specifying hot and cold jets using the ``JET_CHI`` parameter
* Add support for restarting with particles
* Improved tool scripts


**Fixes:**

* Revert addition of 1/2 in non-relativistic jet density-power equation.

**Compatibility:**

* PLUTO = 4.3.6

5.0
---

**New features:**

* Use the new method of specifying environment.ini on the command line (introduced in PLUTO 4.3.5)
* External environment velocity can now be enabled or disabled through user configuration (see :ref:`Configuration`).
* Timing information is shown by default
* Added support for the TAUB EOS (manually using non-relativistic gamma when relevant)
* Added ``start-jupyter`` script to the tools directory.
  This will start a jupyterlab session on an HPC site (kunanyi or gadi).
  You can specify a bunch of job options as arguments to the script.

**Fixes:**

* The ``create-build-script`` now defaults to using  8 cpus for compiling, rather than everything available.
* Fixed a missing factor of 1/2 in the non-relativistic jet-power equation.
