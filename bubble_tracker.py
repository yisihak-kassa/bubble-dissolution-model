"""
This module simulates the dissolutionof a gas bubble rising through a water column.
It includes functions to interpolate temperature profiles, generate pressure profiles, calculate
solubilities using Henry's law, compute terminal velocities of bubbles, and calculate mass transfer
coefficients. The main function, `track_bubble`, integrates the differential equation for bubble
dissolution using Euler's method over space to simulate the bubble's behavior as it ascends.
Functions:
    interpolate_temperature_profile(depth, dz, profile) -> tuple[np.ndarray, np.ndarray]:
        Generate temperature values at midpoint of each depth interval from input temperature profile.
    generate_pressure_profile(points, temperatures, atmospheric_pressure=Q_(1.01325, "bar"), density_function=None) -> Q_:
        Compute pressure profile based on depth intervals and temperatures.
    calculate_solubilities(t) -> Q_:
        Calculate solubilities using Henry's law.
    calculate_terminal_velocity(diameter, rho_water, mu_water) -> Q_:
        Compute terminal velocity of a rising bubble.
    calculate_mass_transfer_coefficient(d, v) -> Q_:
        Compute the liquid-side mass transfer coefficient.
    track_bubble(initial_diameter, release_depth, initial_composition_by_gas, temperature_profile, water_density, water_viscosity, dz=Q_(0.05, "m")) -> dict:
        Simulate the dissolution and tracking of a gas bubble rising through a water column.
"""

import numpy as np
import pint

ureg = pint.UnitRegistry()
Q_ = ureg.Quantity

# -------------------------------
# Physical Constants
# -------------------------------
g = Q_(9.80665, "m/s²")
R = Q_(0.0831446261815324, "L*bar/(mol*K)")

# -------------------------------
# Gas Properties
# -------------------------------
GASES = ["methane", "oxygen", "nitrogen", "argon", "carbon_dioxide", "helium"]

# Henry's Law Coefficients for a taylor series expansion of the
# Clausius-Clapeyron equation as suggested by Boehrer et al. (2021)
U1 = np.log(986.9 / 22391)
U2 = np.log(1 / 1.01325)
U3 = np.log(3.96e-4)

HENRY_COEFFS = np.array(
    [
        [-68.8862, 101.4956, 28.7314, U1],  # Methane
        [-58.3877, 85.8079, 23.8439, U1],  # Oxygen
        [-59.6274, 85.7661, 24.3696, U1],  # Nitrogen
        [-55.66578, 82.0262, 22.5929, U1],  # Argon
        [-58.0931, 90.5069, 22.2940, U2],  # CO2
        [0, 0, 0, U3],  # Helium
    ]
)

DIFFUSION_COEFFS = Q_(
    0.95 * np.array([1.67e-5, 2.29e-5, 1.99e-5, 1.98e-5, 1.92e-5, 6.7e-5]), "cm²/s"
)

ATMOSPHERIC_MOLE_FRACTIONS = Q_(
    [0.0000019, 0.2095, 0.7808, 0.0093, 0.00041, 5.24e-6], "dimensionless"
)
ATMOSPHERIC_PARTIAL_PRESSURES = Q_(
    [0.0000019, 0.2095, 0.7808, 0.0093, 0.00041, 5.24e-6], "atm"
)

# Mass transfer coefficients derived from that of oxygen
K = Q_((4e-4 * (DIFFUSION_COEFFS / DIFFUSION_COEFFS[1]) ** 0.5).m, "m/s")


def interpolate_temperature_profile(
    depth: Q_, dz: Q_, profile
) -> tuple[np.ndarray, np.ndarray]:
    """
    Generate temperature values at midpoint of each depth interval of length dz
    from input temperature profile.

    Parameters:
    depth (Q_): Total depth (units: length)
    dz (Q_): Depth interval step (units: length)
    temperature_profile (Q_, callable, list): Temperature specification

    Returns:
    tuple[np.ndarray, np.ndarray]: (intervals, temperatures)
        - intervals (np.ndarray): Interval lengths (units: length)
        - temperatures (np.ndarray): Midpoint temperatures (units: temperature)
    """
    # Convert units
    depth = depth.to("m").magnitude
    dz = dz.to("m").magnitude

    num_points = int(depth // dz)
    remainder = depth % dz

    if remainder > 1e-9:  # If there's a remainder, we need an extra point
        points = np.linspace(0, depth, num_points + 2)
    else:
        points = np.linspace(0, depth, num_points + 1)

    if isinstance(
        profile, list
    ):  # if temperature profile is a set of discrete measurements
        depths, temps = zip(*sorted((d.to("m").m, t.to("degC").m) for d, t in profile))
        interpolated_temp = np.interp(
            points, depths, temps, left=temps[0], right=temps[-1]
        )

    elif callable(profile):  # if temperature profile is a function
        interpolated_temp = np.array([profile(z).to("degC").m for z in points])

    else:  # if temperature profile is a constant value
        if not isinstance(profile, Q_):
            profile = Q_(profile, "degC")
        interpolated_temp = np.full_like(points, profile.to("degC").m)

    return Q_(points, "m"), Q_(interpolated_temp, "degC")


def generate_pressure_profile(
    points: Q_,
    temperatures: Q_,
    atmospheric_pressure: Q_ = Q_(1.01325, "bar"),
    density_function: callable = None,
) -> Q_:
    """
    Compute pressure profile.

    Parameters:
    atmospheric_pressure (Q_): Surface pressure (units: pressure)
    intervals (Q_): Depth intervals (units: length)
    temperatures (Q_): Midpoint temperatures (units: temperature)
    density_function (callable): Temperature -> density (units: mass/volume)

    Returns:
    Q_: Total pressures (units: pressure)
    """
    t = temperatures.to("degC").m
    densities = Q_(
        (
            999.84847
            + 6.337563e-2 * t
            - 8.523829e-3 * t**2
            + 6.943248e-5 * t**3
            - 3.821216e-7 * t**4
        ),
        "kg/m³",
    )

    hydrostatic = densities * g * points

    return points, (atmospheric_pressure + hydrostatic)


def calculate_solubilities(t):
    """Calculate solubilities using Henry's law."""
    t = t.to("K").magnitude
    ones = np.ones_like(t)
    inv_t = 100 / t
    log_t = np.log(t / 100)
    f_vec = np.vstack([ones, inv_t, log_t, ones])

    return Q_(np.exp(HENRY_COEFFS @ f_vec).T, "mol/(L.bar)")


def calculate_terminal_velocity(diameter, rho_water, mu_water):
    r = diameter.to("m").magnitude / 2
    """Compute terminal velocity of a rising bubble."""
    if r < 7 * 10**-4:
        v = 4474 * r**1.357
    elif r <= 5.1 * 10**-3:
        v = 0.23
    else:
        v = 4.202 * r**0.547
    return Q_(v, "m/s")


import math


def calculate_mass_transfer_coefficient(d, v):
    """Compute the liquid-side mass transfer coefficient"""

    d = d.to("cm").magnitude
    v = v.to("cm/s").magnitude
    D = DIFFUSION_COEFFS.to("cm²/s").magnitude

    c = 0.2 + 0.8 * d
    k = 0.0113 * np.sqrt(v / c) * D**0.5

    return Q_(k, "m/s")  # k


def track_bubble(
    initial_diameter,
    release_depth,
    initial_composition_by_gas,
    temperature_profile,
    water_density,
    water_viscosity,
    dz=Q_(0.05, "m"),
):

    rho_water = water_density.to("kg/m³").magnitude
    mu_water = water_viscosity.to("Pa·s").magnitude
    intervals, temperatures = interpolate_temperature_profile(
        release_depth, dz, temperature_profile
    )
    temperatures = temperatures.to("K")

    depths, pressures = generate_pressure_profile(intervals, temperatures)
    solubility_profile_by_gas = calculate_solubilities(temperatures)

    bubble_diameter = initial_diameter.to("mm")

    bubble_volume = 4 / 3 * np.pi * (bubble_diameter / 2) ** 3
    bubble_composition_by_gas = np.array(initial_composition_by_gas)
    max_steps = len(temperatures)
    diameters = np.zeros(max_steps)
    n = Q_(np.zeros(max_steps), "mol")
    compositions = np.zeros((max_steps, len(bubble_composition_by_gas)))

    diameters[0] = initial_diameter.magnitude
    compositions[0] = initial_composition_by_gas

    temperature = temperatures[-1]
    pressure = pressures[-1]
    bubble_mass = (pressure * bubble_volume / (R * temperature)).to("mol")
    n[0] = bubble_mass

    bubble_mass_by_gas = bubble_mass * bubble_composition_by_gas
    for step in range(1, max_steps):
        solubility_by_gas = solubility_profile_by_gas[-step]
        bubble_surface_area = 4 * np.pi * (bubble_diameter / 2) ** 2

        bubble_pressure_by_gas = pressure * bubble_composition_by_gas
        bubble_velocity = calculate_terminal_velocity(
            bubble_diameter, rho_water, mu_water
        )

        k_by_gas = calculate_mass_transfer_coefficient(bubble_diameter, bubble_velocity)

        pressure_gradient = bubble_pressure_by_gas - ATMOSPHERIC_PARTIAL_PRESSURES
        concentration_gradient_by_gas = solubility_by_gas * pressure_gradient
        mass_transferred_by_gas = (
            -k_by_gas
            * concentration_gradient_by_gas
            * bubble_surface_area
            / bubble_velocity
            * dz
        )

        bubble_mass_by_gas = np.maximum(bubble_mass_by_gas + mass_transferred_by_gas, 0)
        bubble_mass = sum(bubble_mass_by_gas)
        n[step] = bubble_mass.to("mol")
        bubble_composition_by_gas = bubble_mass_by_gas / bubble_mass

        pressure = pressures[-(step + 1)]
        bubble_volume = bubble_mass * R * temperature / pressure
        bubble_diameter = 2 * (3 * bubble_volume.to("m³") / (4 * np.pi)) ** (1 / 3)

        diameters[step] = bubble_diameter.to("mm").magnitude
        compositions[step, :] = bubble_composition_by_gas.m

        step = step + 1

    return {
        "diameters": diameters,
        "depths": depths.m,
        "pressures": pressures.m,
        "compositions": compositions,
        "bubble_masses": n,
    }


if __name__ == "__main__":
    import time

    initial = time.time()
    result = track_bubble(
        initial_diameter=Q_(2.13, "mm"),
        release_depth=Q_(4.69, "m"),
        initial_composition_by_gas=[
            0,
            0,
            1,
            0,
            0,
            0,
        ],  # Gases in the order "methane", "oxygen", "nitrogen", "argon", "carbon_dioxide"
        temperature_profile=Q_(19, "°C"),
        water_density=Q_(1000, "kg/m³"),
        water_viscosity=Q_(0.001, "Pa·s"),
    )
    final = time.time()
    print(f"Simulation took {final - initial:.2f} seconds.")

    print("Simulation Completed Successfully!")

    print(f"Final Diameter: {result['diameters'][-1]:.2f} mm")
    print(f"Depth: {result['depths'][-1]:.3f} m")
    print("Final Composition: ", end="")
    print(", ".join([f"{comp*100:.2f}%" for comp in result["compositions"][-1]]))
