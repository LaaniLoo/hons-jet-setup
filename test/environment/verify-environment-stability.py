#!/usr/bin/env python3

import os
import sys
if os.path.exists(os.path.expanduser('~/plutokore')):
    sys.path.append(os.path.expanduser('~/plutokore'))
else:
    sys.path.append(os.path.expanduser('~/uni/plutokore'))
import plutokore as pk
import argparse

# utilities
from numba import jit
from pathlib import Path
import h5py
from IPython.display import display,clear_output
from tabulate import tabulate
import colorful as cf

# science imports
import numpy as np
import scipy.interpolate
import scipy.integrate
from scipy.integrate import trapz

# matplotlib imports
import matplotlib as mpl
import matplotlib.pyplot as plot

# astropy imports
from astropy.table import Table
from astropy import units as u # Astropy units
from astropy import cosmology as cosmo # Astropy cosmology
from astropy import constants as const # Astropy constants
from astropy.convolution import convolve, Gaussian2DKernel # Astropy convolutions

from plutokore.jet import UnitValues

unit_length = 1 * u.kpc
unit_density = (0.60364 * u.u / (u.cm ** 3)).to(u.g / u.cm ** 3)
unit_speed = const.c
unit_time = (unit_length / unit_speed).to(u.Myr)
unit_pressure = (unit_density * (unit_speed ** 2)).to(u.Pa)
unit_mass = (unit_density * (unit_length ** 3)).to(u.kg)
unit_energy = (unit_mass * (unit_length**2) / (unit_time**2)).to(u.J)

uv = UnitValues(
    density=unit_density,
    length=unit_length,
    time=unit_time,
    mass=unit_mass,
    pressure=unit_pressure,
    energy=unit_energy,
    speed=unit_speed,
)

def verify_environment(sim_dir):
    times = pk.simulations.get_times(sim_dir)
    initial_output = pk.simulations.load_timestep_data(0, sim_dir)

    headers = [['Output', 'Time (Myr)', 'Equiv. Power (W)', 'Density % Diff', 'Mach', 'Vel (km/s)', 'Passed']]
    results = []

    comparison_kinetic_energy = 0
    comparison_power = 0

    for output_number, time in enumerate(times):
        output = pk.simulations.load_timestep_data(output_number, sim_dir)
        cur_res = []

        passed = True

        cur_res.append(output_number)
        cur_res.append(output.SimTime * uv.time.value)

        # power check
        (equiv_power, passed, comparison_power) = check_power(sim_output = output, comparison_power = comparison_power)

        if passed:
            cur_res.append(f'{equiv_power:0.03g}')
        else:
            cur_res.append(cf.red(f'{equiv_power:0.03g}'))

        # central density check
        (density_diff, passed) = check_density(sim_output = output, initial_output = initial_output)
        if passed:
            cur_res.append(f'{density_diff:0.03g}')
        else:
            cur_res.append(cf.red(f'{density_diff:0.03g}'))

        # mach number
        (mach, vel, passed) = check_velocity(sim_output = output)

        if passed:
            cur_res.append(f'{mach:0.03g}')
            cur_res.append(f'{vel:0.03g}')
        else:
            cur_res.append(cf.red(f'{mach:0.03g}'))
            cur_res.append(cf.red(f'{vel:0.03g}'))

        # check passed
        if passed:
            cur_res.append(cf.bold_green('Yes'))
        else:
            cur_res.append(cf.bold_red('No'))

        # append to main results list
        results.append(cur_res)

    return headers + results

def check_density(*, sim_output, initial_output):
    if (sim_output.geometry == 'CARTESIAN'):
        initial_dens = initial_output.rho[initial_output.n1 // 2, initial_output.n2 // 2, initial_output.n3 // 2]
        current_dens = sim_output.rho[sim_output.n1 // 2, sim_output.n2 // 2, sim_output.n3 // 2]
    elif (sim_output.geometry == 'SPHERICAL'):
        initial_dens = initial_output.rho[0, 0]
        current_dens = sim_output.rho[0, 0]

    density_diff = (initial_dens - current_dens) / initial_dens * 100

    passed = True
    if density_diff > 1:
        passed = False

    return (density_diff, passed)

def check_velocity(*, sim_output):
    if 'vx3' in sim_output.vars:
        vel = np.sqrt(sim_output.vx1 ** 2 + sim_output.vx2 ** 2 + sim_output.vx3 ** 2) * uv.speed
    else:
        vel = np.sqrt(sim_output.vx1 ** 2 + sim_output.vx2 ** 2) * uv.speed

    sound_speed = np.sqrt((3.0 / 5.0) * (sim_output.prs * uv.pressure) / (sim_output.rho * uv.density)).to(u.km / u.s)
    max_mach = np.max((vel / sound_speed).si)
    max_vel = np.max(vel).to(u.km / u.s)

    passed = True
    if max_mach > 0.1:
        passed = False

    return (max_mach.value, max_vel.value, passed)

def check_power(*, sim_output, comparison_power):
    if sim_output.SimTime == 0:
        return (0, True, 0)
    vol = pk.simulations.calculate_cell_volume(sim_output) * (uv.length ** 3)
    if 'vx3' in sim_output.vars:
        vel = np.sqrt(sim_output.vx1 ** 2 + sim_output.vx2 ** 2 + sim_output.vx3 ** 2) * uv.speed
    else:
        vel = np.sqrt(sim_output.vx1 ** 2 + sim_output.vx2 ** 2) * uv.speed
    mass = sim_output.rho * uv.density * vol
    kin_energy = 0.5 * mass * (vel ** 2)
    kin_energy_density = 0.5 * sim_output.rho * uv.density * (vel ** 2)

    max_kin_energy_density = kin_energy_density.max().to(u.J / (u.kpc ** 3)).value

    total_energy = kin_energy.sum().to(u.J)
    equiv_power = (total_energy / (sim_output.SimTime * uv.time)).to(u.W)

    passed = True

    if equiv_power > 0:
        if comparison_power == 0:
            comparison_power = equiv_power
        elif equiv_power > 10*comparison_power:
            passed = False
    return (equiv_power.value, passed, comparison_power)

def verify(sim_dir):
    results = verify_environment(sim_dir)

    print(tabulate(results, headers = 'firstrow', tablefmt = 'github'))

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('simulation_directory')

    args = ap.parse_args()

    verify_environment(args.simulation_directory)

if __name__ == "__main__":
    main()
