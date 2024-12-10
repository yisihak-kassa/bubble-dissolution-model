import math
import constants
from bubble import Bubble
import utils
from water import Water
import units

GAS_FRACTIONS = constants.GAS_FRACTIONS
GASES = constants.GASES

def calculate_concentration_gradient(bubble, temperature, atmospheric_pressure):
    """Calculate the concentration gradient for each gas."""
    return {
        gas: utils.calculate_solubility(gas, temperature)
        * (
            bubble.gas_fractions.get(gas, 0) * bubble.pressure
            - GAS_FRACTIONS.get(gas, 0) * atmospheric_pressure
        )
        for gas in GASES
    }

def calculate_mass_transfer(bubble, concentration_gradient, dz):
    """Calculate the change in mass for each gas."""
    surface_area = 4 * math.pi * (bubble.diameter / 2) ** 2
    return {
        gas: (
            -bubble.k_liquid.get(gas, 0)
            * concentration_gradient.get(gas, 0)
            * surface_area
            / bubble.terminal_velocity
            * dz
        ).to(units.mol)
        for gas in GASES
    }

def update_bubble_properties(bubble, dz, temperature, atmospheric_pressure, water):
    """Update bubble properties as it ascends."""
    concentration_gradient = calculate_concentration_gradient(
        bubble, temperature, atmospheric_pressure
    )
    change_in_mass_by_gas = calculate_mass_transfer(bubble, concentration_gradient, dz)
    change_in_mass = sum(change_in_mass_by_gas.values())
    new_number_of_moles = bubble.number_of_moles + change_in_mass

    mass_by_gas = bubble.mass_by_gas
    new_mass_by_gas = {
        gas: (mass_by_gas.get(gas, 0) + change_in_mass_by_gas.get(gas, 0))
        for gas in change_in_mass_by_gas.keys()
    }
    total_mass = sum(new_mass_by_gas.values())
    gas_fractions = {gas: new_mass_by_gas.get(gas, 0) / total_mass for gas in GASES}

    pressure = water.calculate_pressure(bubble.depth - dz, atmospheric_pressure)

    try:
        volume = (
            new_number_of_moles * constants.R * bubble.water.temperature.to(units.K) / pressure
        )
        if volume < 0:
            raise ValueError("Volume must be a non-negative value.")

        diameter = utils.calculate_bubble_diameter(volume)
    except ValueError as e:
        raise RuntimeError(f"Error calculating volume or diameter: {e}")

    return diameter, gas_fractions

def bubble_history(initial_diameter, initial_composition, temperature, release_depth, atmospheric_pressure):
    """
    Simulate the history of a bubble as it ascends through water.

    Args:
        initial_diameter (float): Initial bubble diameter.
        initial_composition (dict): Initial gas composition.
        temperature (float): Water temperature.
        release_depth (float): Depth of bubble release.
        atmospheric_pressure (float): Atmospheric pressure.

    Returns:
        dict: A dictionary mapping depth to Bubble objects.
    """
    water = Water(temperature)
    z = release_depth
    dz = 0.1 * units.m

    history = {z: Bubble(initial_diameter, z, initial_composition, water)}

    while z > 0:
        bubble = history[z]
        try:
            diameter, gas_fractions = update_bubble_properties(
                bubble, dz, temperature, atmospheric_pressure, water
            )
        except RuntimeError as e:
            print(e)
            print("Terminating bubble history calculation.")
            return history

        z -= dz
        history[z] = Bubble(diameter, z, gas_fractions, water)

    return history

if __name__ == "__main__":
    fractions = {
        "nitrogen": 0,
        "oxygen": 1,
        "argon": 0,
        "methane": 0,
        "carbon_dioxide": 0,
    }

    temp = units.in_celsius(19)
    d0 = 2.3 * units.mm
    z0 = 4.75 * units.m
    p_atm = 1.01325 * units.bar

    history = bubble_history(
        initial_diameter=d0,
        initial_composition=fractions,
        temperature=temp,
        release_depth=z0,
        atmospheric_pressure=p_atm,
    )

    output_str = "\n".join(
        f"d: {bubble.diameter:.2f} depth:{depth:.2f} {utils.dict_as_string(bubble.gas_fractions)}"
        for depth, bubble in history.items()
    )
    print(output_str)
