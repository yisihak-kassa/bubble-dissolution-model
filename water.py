import math
import units
import constants


class Water:
    """
    A class to represent water properties, including temperature, density, and viscosity.
    """

    def __init__(self, temperature):
        """
        Initialize the Water object with temperature.

        :param temperature: Water temperature in degrees Celsius.
        :param pressure: Atmospheric pressure in atm.
        """
        self.temperature = temperature
        #: diffussion exponent (0.5 for clean water)
        self.n = 0.5
        self._density = None
        self._saturated_density = None
        self._air_free_density = None
        self._viscosity = None

        # Initialize calculated properties
        self._update_properties()

    def _update_properties(self):
        """
        Update all dependent properties based on the current state.
        """
        self._density = self.calculate_density()
        self._saturated_density = self.calculate_density(saturated=True)
        self._air_free_density = self.calculate_density(saturated=False)
        self._viscosity = self.calculate_dynamic_viscosity()

    @property
    def density(self) -> float:
        """Get the water density [kg/m^3]"""
        return self._density

    @property
    def saturated_density(self) -> float:
        """Get the density of saturated water in g/ml."""
        return self._saturated_density

    @property
    def air_free_density(self) -> float:
        """Get the density of air-free water in g/ml."""
        return self._air_free_density

    @property
    def viscosity(self) -> float:
        """Get the dynamic viscosity of water in kg/(m·s)."""
        return self._viscosity

    def calculate_dynamic_viscosity(self) -> float:
        """
        Calculate the dynamic viscosity of water as a function of temperature
        using the Vogel-Fulcher-Tammann equation.

        :return: Dynamic viscosity in kg/(m·s).
        """
        A = 0.02939 * units.mPa_s  # Convert from mPa·s to kg/m·s
        B = 507.88 * units.K  # K
        C = 149.3 * units.K  # K

        T = self.temperature.to(units.K)
        viscosity = A * math.exp(B / (T - C))
        print(str(viscosity))
        #     f"Dynamic viscosity of water at {self.temperature.celsius} oC: {viscosity} Pa·s"
        # )
        return viscosity

    def calculate_density(self, saturated: bool = True) -> float:
        """
        Calculate the density of water at a given temperature.

        :param saturated: Whether the water is saturated with air.

        Returns:
        - float: The density of water [kg/m^3].
        """
        temp = self.temperature.to(units.C)
        temperature = temp.magnitude
        if temperature < 0 or temperature > 100:
            raise ValueError(
                "Temperature should be in the range 0-100°C for accurate density calculation."
            )

        density = None
        if saturated:
            density = (
                999.84847
                + 6.337563e-2 * temperature
                - 8.523829e-3 * temperature**2
                + 6.943248e-5 * temperature**3
                - 3.821216e-7 * temperature**4
            )
        else:
            density = (
                999.85308
                + 6.32693e-2 * temperature
                - 8.523829e-3 * temperature**2
                + 6.943248e-5 * temperature**3
                - 3.821216e-7 * temperature**4
            )

        return density * units.kgpm3

    def calculate_pressure(
        self, depth: float, atmospheric_pressure: float = 1.01325 * units.bar
    ) -> float:
        """
        Calculate the total pressure at a given depth in water in mbar.

        :param depth: Depth of water in meters.
        :param atmospheric_pressure: Atmospheric pressure in mbar.
        :return: Total pressure in mbar.
        """
        hydrostatic_pressure = self.density * constants.g * depth
        # print(str(self.density))
        return (hydrostatic_pressure + atmospheric_pressure).to(units.bar)
