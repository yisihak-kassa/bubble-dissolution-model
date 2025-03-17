"""Module: data_utils
This module provides utility functions for calculating various properties related to gas dissolution in water,
including adjustments for oxygen content, depth-adjusted gas volumes, partial pressures, and Henry's Law constants.
It also includes functions for converting temperature units, calculating water density, and determining travel distances
and capture depths in experimental setups.
Functions:
    - adjust_oxygen_content(oxygen_content, gas_volume, testing_temperature): Adjusts the oxygen content of samples based on gas volume and testing temperature.
    - calculate_depth_adjusted_volume(measured_volume, target_depth, reference_pressure, reference_depth, water_density): Calculates the depth-adjusted gas volume based on a known volume at a reference depth.
    - get_gas_volume(exetainer_mass_with_water, sample_mass, water_density): Calculates the gas volume in samples.
    - calculate_partial_pressure(gas_fraction, total_pressure): Calculates the partial pressure of a gas.
    - calculate_pressure(depth, water_density, atmospheric_pressure): Calculates the total pressure at a given depth in water.
    - calculate_k_H_cp(gas, T): Calculates Henry's Law constant K_H_cp in mol/L/bar.
    - calculate_k_H_cc(gas, T): Calculates dimensionless Henry's Law constant K_H_cc.
    - calculate_water_density(temperature, saturated): Calculates the density of water at a given temperature.
    - convert_celsius_to_kelvin(celsius): Converts temperature from Celsius to Kelvin.
    - calculate_diameter(volume): Calculates the diameter of a sphere given its volume.
    - calculate_travel_distance(depth, free_board): Calculates the travel distance of a gas bubble in an experimental setup.
    - calculate_capture_depth(nominal_depth, free_board): Calculates the capture depth of a gas bubble in an experimental setup.
"""

import math
import units

# Constants
ATMOSPHERIC_OXYGEN_FRACTION = 0.2095
EXETAINER_VOLUME = 12.0  # [mL]
GRAVITATIONAL_ACCELERATION = 9.80665  # Gravitational constant in m/s²
R = 0.0831446261815324  # ideal gas constant [L⋅bar/mol.K]
DYNAMIC_VISCOSITY = 0.001002  # Dynamic viscosity of water at 19°C in Pa·s

# Coefficients for the calculation of Henry's Law constant K_H, according to Boehrer et. al. ()
GASES = ["nitrogen", "oxygen", "argon", "methane", "carbon_dioxide"]
GAS_COEFFICIENTS = {
    "nitrogen": {"A1": -59.6274, "A2": 85.7661, "A3": 24.3696, "u": 986.9 / 22391},
    "oxygen": {"A1": -58.3877, "A2": 85.8079, "A3": 23.8439, "u": 986.9 / 22391},
    "argon": {"A1": -55.66578, "A2": 82.0262, "A3": 22.5929, "u": 986.9 / 22391},
    "methane": {"A1": -68.8862, "A2": 101.4956, "A3": 28.7314, "u": 986.9 / 22391},
    "carbon_dioxide": {"A1": -58.0931, "A2": 90.5069, "A3": 22.2940, "u": 1 / 1.01325},
}

GAS_FRACTIONS = {
    "nitrogen": 0.7808,
    "oxygen": 0.2095,
    "argon": 0.0093,
    "methane": 0.00019,
    "carbon_dioxide": 0.00041,
}


def adjust_oxygen_content(oxygen_content, gas_volume, testing_temperature):
    """
    Adjust the oxygen content of samples based on gas volume and testing temperature.
    Takes account of equilibration of between gas phase and water.

    Parameters:
        oxygen_content (float): Oxygen content in the gas phase (percentage).
        gas_volume (float): Volume of gas present (mL).
        testing_temperature (float): Temperature at which the measurement is taken (°C).

    Returns:
        float: Adjusted oxygen content as a percentage.
    """

    equilibrium_oxygen_content = oxygen_content / 100  # Convert percentage to fraction

    # Calculate water volume based on gas volume
    water_volume = EXETAINER_VOLUME - gas_volume  # [mL]
    k_h = calculate_k_H_cc(
        "oxygen", testing_temperature
    )  # Henry's law constant for oxygen

    # Calculate gas phase and dissolved oxygen volumes
    gas_phase_oxygen_volume = equilibrium_oxygen_content * gas_volume
    equilibrium_dissolved_oxygen_volume = (
        k_h * equilibrium_oxygen_content * water_volume
    )

    # Calculate initial dissolved oxygen volume based on equilibrium with atmospheric gases
    initial_dissolved_oxygen_volume = k_h * GAS_FRACTIONS["oxygen"] * water_volume

    # the flux of oxygen between the gas phase and the water (convention: water to gas flux considered as a positive flux)
    flux = initial_dissolved_oxygen_volume - equilibrium_dissolved_oxygen_volume

    # Adjusted oxygen volume accounting for initial dissolved oxygen
    adjusted_oxygen_volume = gas_phase_oxygen_volume - flux

    # Calculate adjusted oxygen content as a percentage of the gas volume
    adjusted_oxygen_content = adjusted_oxygen_volume / gas_volume

    return adjusted_oxygen_content * 100


def calculate_depth_adjusted_volume(
    measured_volume,
    target_depth,
    reference_pressure=1.01325 * units.bar,
    reference_depth=0,
    water_density=0.997,
):
    """
    This function calculates the depth-adjusted gas volume based on a known volume at a reference depth.

    Parameters:
    - measured_volume: The gas volume at the reference depth (in ml)
    - target_depth: The depth at which you want to calculate the gas volume (in meters)
    - reference_depth: The depth where the initial volume was measured (default is 0 meters, i.e., at the surface)
    - reference_pressure: reference pressure at the reference depth (default is 1013.25 mbar for a reference depth of 0 meters)
    - water_density: The density of the water in g/ml (default is 0.997 g/ml, assuming fresh water)

    Returns:
    - The adjusted gas volume at the target depth.
    """
    reference_pressure = calculate_pressure(reference_depth, water_density)
    # Calculate pressure at the target depth
    target_pressure = calculate_pressure(
        target_depth, water_density
    )  # Target pressure in mbar

    # Calculate new volume using Boyle's Law
    adjusted_volume = (
        reference_pressure / target_pressure * measured_volume
    )  # New volume in mL
    # print(target_pressure)
    return adjusted_volume


# calculate gas volume in samples
def get_gas_volume(exetainer_mass_with_water, sample_mass, water_density):
    v_gas = (exetainer_mass_with_water - sample_mass) / water_density
    return v_gas


def calculate_partial_pressure(gas_fraction, total_pressure=1.01325):
    return gas_fraction * total_pressure


# Function to calculate total pressure underwater
def calculate_pressure(depth, water_density=1.0, atmospheric_pressure=1013.25):
    """
    Calculate the total pressure at a given depth in water in mbar.

    Parameters:
        depth (float): Depth of water in meters.
        water_density (float): Density of water in g/ml (default is 1.0 g/ml for pure water).

    Returns:
        float: Total pressure in mbar.
    """
    g = GRAVITATIONAL_ACCELERATION
    # Convert density from g/mL to kg/m³ (1 g/mL = 1000 kg/m³)
    density_kg_m3 = water_density * 1000
    hydrostatic_pressure = (
        density_kg_m3 * g * depth / 100
    )  # convert Pa to mbar by dividing by 100

    total_pressure = hydrostatic_pressure + atmospheric_pressure
    return total_pressure


# calculate Henry's Law constant K_H_cc in mol/L/bar
def calculate_k_H_cp(gas="oxygen", T=25.0):
    coefficients = GAS_COEFFICIENTS.get(gas)
    A1 = coefficients["A1"]
    A2 = coefficients["A2"]
    A3 = coefficients["A3"]
    u = coefficients["u"]

    T_abs = convert_celsius_to_kelvin(T)

    k_H = math.exp(A1 + A2 * 100 / T_abs + A3 * math.log(T_abs / 100)) * u

    return k_H


# calculate dimensionless Henry's Law constant K_H_cc
def calculate_k_H_cc(gas="oxygen", T=25.0):
    T_abs = convert_celsius_to_kelvin(T)
    k_H = calculate_k_H_cp(gas, T) * R * T_abs

    return k_H


def calculate_water_density(temperature, saturated=True):
    """
    Calculates the density of water (in g/ml) at a given temperature.

    Parameters:
    - temperature: Temperature of water in degrees Celsius (°C)

    Returns:
    - Tuple of density of air free and saturated water respectively in g/ml
    """
    if temperature < 0 or temperature > 100:
        raise ValueError(
            "Temperature should be in the range 0-100°C for accurate density calculation."
        )

    # If the water is assummed to be saturated with air
    if saturated:
        saturated_density = (
            999.84847
            + 6.337563 * 10**-2 * temperature
            - 8.523829 * 10**-3 * temperature**2
            + 6.943248 * 10**-5 * temperature**3
            - 3.821216 * 10**-7 * temperature**4
        ) / 1000
        return saturated_density
    else:
        air_free_density = (
            999.85308
            + 6.32693 * 10**-2 * temperature
            - 8.523829 * 10**-3 * temperature**2
            + 6.943248 * 10**-5 * temperature**3
            - 3.821216 * 10**-7 * temperature**4
        ) / 1000
        return air_free_density


def convert_celsius_to_kelvin(celsius):
    """
    Convert temperature from Celsius to Kelvin.

    Parameters:
        celsius (float): Temperature in degrees Celsius.

    Returns:
        float: Temperature in Kelvin.
    """
    if celsius < -273.15:
        raise ValueError(
            "Temperature in Celsius cannot be below absolute zero (-273.15°C)."
        )

    kelvin = celsius + 273.15
    return kelvin


def calculate_diameter(volume):
    """
    Calculate the diameter of a sphere given its volume.

    Parameters:
    - volume (float): The volume of the sphere.

    Returns:
    - float: The diameter of the sphere.
    """
    if volume < 0:
        raise ValueError("Volume must be a non-negative value.")

    diameter = (6 * volume / math.pi) ** (1 / 3)  # Calculate radius
    return diameter


def calculate_travel_distance(depth, free_board=0.1):
    column_depth = 5.0
    nominal_travel_distance = column_depth - depth
    exetainer_height = 0.04
    funnel_height = 0.17
    diffuser_height = 0.08
    diffuser_level = 0.04
    release_level = funnel_height + diffuser_height - diffuser_level + exetainer_height
    if depth == 0:
        return nominal_travel_distance + exetainer_height - release_level - free_board

    elif depth == 5:
        return 0
    else:
        return nominal_travel_distance + exetainer_height - release_level


def calculate_capture_depth(nominal_depth, free_board=0.05):
    column_depth = 5.0
    exetainer_height = 0.04
    funnel_height = 0.17
    diffuser_height = 0.08

    travel_dist = calculate_travel_distance(nominal_depth)
    if nominal_depth == 5:
        return (
            column_depth
            - free_board
            - funnel_height
            - diffuser_height
            - exetainer_height
        )
    elif nominal_depth == 0:
        return 0.00
    else:
        return column_depth - free_board - travel_dist
