# imports
import math
import units
import constants

# Constants
ATMOSPHERIC_OXYGEN_FRACTION = 0.2095
EXETAINER_VOLUME = 12.0  # [mL]
GRAVITATIONAL_ACCELERATION = constants.g  # Gravitational constant in m/s²
R = constants.R  # ideal gas constant [L⋅bar/mol.K]
DYNAMIC_VISCOSITY = 0.001002  # Dynamic viscosity of water at 19°C in Pa·s


GASES = ["nitrogen", "oxygen", "argon", "methane", "carbon_dioxide"]
# Coefficients for the calculation of Henry's Law constant K_H, according to Boehrer et. al. (2021)
GAS_COEFFICIENTS = {
    "nitrogen": {"A1": -59.6274, "A2": 85.7661, "A3": 24.3696, "u": 986.9 / 22391},
    "oxygen": {"A1": -58.3877, "A2": 85.8079, "A3": 23.8439, "u": 986.9 / 22391},
    "argon": {"A1": -55.66578, "A2": 82.0262, "A3": 22.5929, "u": 986.9 / 22391},
    "methane": {"A1": -68.8862, "A2": 101.4956, "A3": 28.7314, "u": 986.9 / 22391},
    "carbon_dioxide": {"A1": -58.0931, "A2": 90.5069, "A3": 22.2940, "u": 1 / 1.01325},
    # "helium": {"A1": -59.6274, "A2": 85.7661, "A3": 24.3696, "u": 986.9 / 22391},
}

GAS_FRACTIONS = {
    "nitrogen": 0.7808,
    "oxygen": 0.2095,
    "argon": 0.0093,
    "methane": 0.0000019,
    "carbon_dioxide": 0.00041,
}

DIFFUSION_COEFFICIENTS = {
    "nitrogen": 1.99,
    "oxygen": 2.29,
    "argon": 1.98,
    "methane": 1.67,
    "carbon_dioxide": 1.92,
}

MOLAR_VOLUMES = {
    "nitrogen": 24.5,
    "oxygen": 27.9,
    "argon": 29.2,
    "methane": 37.7,
    "carbon_dioxide": 37.3,
}

MOLAR_MASSES = {
    "nitrogen": 28.0134,
    "oxygen": 31.9988,
    "argon": 39.948,
    "carbon_dioxide": 44.0095,
    "methane": 16.04,
    "helium": 4.0026,
    "neon": 20.1797,
    "krypton": 83.798,
    "xenon": 131.293,
    "hydrogen": 2.01588,
    "water_vapor": 18.01528,
}


def dict_as_string(d, s=False, unit=True, multiplier=1):
    string = ""

    if s:
        for key in d.keys():
            string += f"{key}: {d[key]*multiplier:.2e}\t"
    else:
        for key in d.keys():
            string += f"{key}: {d[key]*multiplier:.2f}\t"

    return string


def adjust_oxygen_content(oxygen_content, gas_volume, testing_temperature):
    """
    Adjust the oxygen content based on gas volume and testing temperature.

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
    k_h = calculate_dimensionless_solubility(
        "oxygen", testing_temperature
    )  # Henry's law constant for oxygen

    # Calculate gas phase and dissolved oxygen volumes
    gas_phase_oxygen_volume = equilibrium_oxygen_content * gas_volume
    equilibrium_dissolved_oxygen_volume = (
        k_h * equilibrium_oxygen_content * water_volume
    )

    # Calculate initial dissolved oxygen volume based on equilibrium with atmospheric gases
    initial_dissolved_oxygen_volume = k_h * GAS_FRACTIONS["oxygen"] * water_volume

    # the flux of oxygen between the gas phase and the water (convention: water to gas considered as a positive flux)
    flux = initial_dissolved_oxygen_volume - equilibrium_dissolved_oxygen_volume

    adjusted_oxygen_volume = gas_phase_oxygen_volume - flux

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


# calculate Henry's Law constant K_H_cc in mol/L/bar
def calculate_solubility(gas="oxygen", T=units.in_celsius(25.0)):
    coefficients = GAS_COEFFICIENTS.get(gas)
    A1 = coefficients["A1"]
    A2 = coefficients["A2"]
    A3 = coefficients["A3"]
    u = coefficients["u"]

    T_abs = T.to(units.K).magnitude

    k_H = math.exp(A1 + A2 * 100 / T_abs + A3 * math.log(T_abs / 100)) * u

    return k_H * units.molplbar


# calculate dimensionless Henry's Law constant K_H_cc
def calculate_dimensionless_solubility(gas="oxygen", T=units.in_celsius(25)):
    T_abs = T.to(units.K)
    k_H = calculate_solubility(gas, T) * constants.R * T_abs
    return k_H


def calculate_dynamic_viscosity(T: float) -> float:
    """
    Calculate the dynamic viscosity of water as a function of temperature
    using the Vogel-Fulcher-Tammann equation.

    Parameters:
    - T (float): Temperature in degree celsius.

    Returns:
    - float: Dynamic viscosity in N/ms.
    """
    # Constants
    A = 0.02939 / 1000  # Convert from mPa·s to kg/m·s # mPa·s
    B = 507.88  # K
    C = 149.3  # K

    T = T + 273.15
    # Calculate dynamic viscosity
    viscosity = A * math.exp(B / (T - C))
    return viscosity


# calculate diffusion coefficient in cm2/s
def calculate_diffusion_coefficient(mu, V_i):
    mu = mu.to(units.cp)
    V_i = V_i.to(units.cm3pmol)
    return 13.26 * 10 ** (-5) / (mu**1.14 * V_i**0.589)


# calculate drag coefficient
def calculate_drag_coefficient(Re):
    C_D = 24 / Re + 3 / Re**0.5 + 0.34
    return C_D


# calculate Reynolds number
def calculate_reynolds_number(rho: float, v, d: float, mu: float) -> float:
    return rho * v * d / mu


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


def calculate_travel_distance(depth, free_board=0.08):
    column_depth = 5.0

    exetainer_height = 0.04
    funnel_height = 0.17
    diffuser_height = 0.08
    diffuser_level = 0.04

    if depth == 5:
        return funnel_height + diffuser_height - diffuser_level + exetainer_height
    else:
        return column_depth - depth - diffuser_level + exetainer_height


def calculate_capture_depth(nominal_depth, free_board=0.08):
    ass = {0.0: 0, 1.0: 0.92, 2.0: 1.92, 3.0: 2.92, 4.0: 3.92, 5.0: 4.67}
    column_depth = 5.0
    # calculate the level of water in the exetainer
    exetainer_height = 0.04
    funnel_height = 0.16
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
    return ass[nominal_depth]
