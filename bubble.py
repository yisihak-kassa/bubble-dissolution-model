import math
from typing import Dict
from scipy.optimize import newton

import utils
from water import Water
import units
import constants

R = constants.R
MOLAR_MASSES = utils.MOLAR_MASSES
g = constants.g


class Bubble:
    def __init__(
        self,
        diameter,
        depth,
        gas_fractions,
        water: Water,
    ):
        self.diameter = diameter
        self.depth = depth
        self.gas_fractions = gas_fractions
        self.water: Water = water

    @property
    def k_liquid(self):
        n = self.water.n  # diffusion exponent
        v_b = self.terminal_velocity.to(units.cmps).magnitude
        D = {
            gas: self.diffusion_coefficients.get(gas).to(units.cm2ps).magnitude
            for gas in utils.GASES
        }
        d_e = self.diameter.to(units.cm).magnitude
        if d_e < 0.5:
            k_l = {
                gas: 0.0113
                * math.sqrt(v_b / (0.45 + 0.2 * d_e))
                * D[gas] ** n
                * units.mps
                for gas in utils.GASES
            }
        if d_e >= 0.5 and d_e < 1.3:
            k_l = {gas: 6.5 * D[gas] ** n * units.mps for gas in utils.GASES}
        if d_e >= 1.3:
            k_l = {gas: 6.94 * d_e**-0.25 * units.mps for gas in utils.GASES}
        return k_l

    @property
    def terminal_velocity(self) -> float:
        """
        Calculates the terminal velocity of the bubble using Newton's method.
        """
        rho_L = self.water.density.to(units.kgpm3).magnitude  # Water density [kg/m³]
        rho_G = self.density.to(units.kgpm3).magnitude  # Gas bubble density [kg/m^3]
        d_e = self.diameter.to(units.m).magnitude  # Bubble diameter [m]
        mu = self.water.viscosity.to(units.kgpms).magnitude

        def governing_equation(v_b) -> float:
            """
            Governing equation for terminal velocity of a gas bubble.
            """
            Re = utils.calculate_Re(rho_L, v_b, d_e, mu)
            # print(f"Reynold's {str(Re)}")
            C_D = utils.calculate_C_D(Re)
            return v_b**2 - (4 * d_e * g.magnitude * (1 - rho_G / rho_L)) / (3 * C_D)

        # Use Newton's method to solve for terminal velocity (v_b)
        return newton(governing_equation, 0.5) * units.mps

    @property
    def diffusion_coefficients(self):
        mu = self.water.viscosity.to(units.cP).magnitude  # convert from Pa·s to cP
        V = utils.MOLAR_VOLUMES
        diffusion_coefficients = {
            gas: 13.26 * 10**-5 / (mu**1.14 * V.get(gas, 0) ** 0.589) * units.cm2ps
            for gas in utils.GAS_FRACTIONS
        }
        # print(f"diffusion coefficients: {self._diffusion_coefficients}")
        # return constants.DIFFUSION_COEFFICIENTS
        return diffusion_coefficients

    @property
    def mass_by_gas(self):
        return {
            gas: self.gas_fractions.get(gas, 0) * self.number_of_moles
            for gas in self.gas_fractions
        }

    @property
    def number_of_moles(self) -> float:
        mass = (
            self.water.calculate_pressure(self.depth)
            * self.bubble_volume
            / constants.R
            / self.water.temperature.to(units.K)
        )
        return mass

    @property
    def bubble_volume(self) -> float:
        """
        Returns:
        - float: The volume of the gas bubble [cm^3].

        """
        bubble_volume = math.pi * self.diameter**3 / 6
        return bubble_volume

    @property
    def pressure(self) -> float:
        pressure = self.water.calculate_pressure(self.depth)
        return pressure

    @property
    def density(self) -> Dict[str, float]:
        """
        Calculates and caches the density of the gas bubble using the ideal gas law.
        """
        return self.calculate_density()

    @property
    def molar_volume(self):
        """
        Calculate the molar volume of an ideal gas using the ideal gas law.

        Parameters:
        - T (float): Temperature [Kelvin]
        - P (float): Pressure [mbar]

        Returns:
        - float: Molar volume [cm³/mol]
        """
        T = self.water.temperature.to(units.K)  # Convert temperature to Kelvin
        P = self.pressure
        # print("P:", P, "bar")  # Pressure in mbar
        # Calculate molar volume
        Vm = (R * T) / P * 10**3  # m^3 to cm^3
        return Vm

    def calculate_density(self) -> float:
        """
        Calculates the gas mixture density using the ideal gas law.

        Returns:
        - float: The density of the gas mixture [kg/m^3].
        """
        T = self.water.temperature.to(units.K)  # Convert temperature to Kelvin

        P_tot = self.water.calculate_pressure(self.depth)  # Pressure in bar
        # print(f"ptotal {P_tot}")

        # Calculate gas density using the ideal gas law
        density = sum(
            fraction * P_tot * MOLAR_MASSES[gas] * units.gpmol
            for gas, fraction in self.gas_fractions.items()
        ) / (constants.R * T)
        # print(f"bubble density: {density}")
        return density


if __name__ == "__main__":
    import numpy as np

    dia = (
        np.array([1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]) * units.mm
    )  # Bubble diameter in mm
    surrounding_water = Water(
        temperature=units.in_celsius(19.0)
    )  # Water object with temperature
    z = 5 * units.m  # Depth in meters
    mix = {"oxygen": 1}  # Gas mixture composition (e.g., 100% oxygen)
