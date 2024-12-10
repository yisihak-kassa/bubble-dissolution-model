import units

g = 9.80665 * units.mps2  # m/s²
R = 0.0831446261815324 * units.LbarpmolK
AVOGADRO_NUMBER = 6.02214076e23

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

DIFFUSION_COEFFICIENTS = {
    "nitrogen": 1.99 * 10**-5 * units.cm2ps,
    "oxygen": 2.29 * 10**-5 * units.cm2ps,
    "argon": 1.98 * 10**-5 * units.cm2ps,
    "methane": 1.67 * 10**-5 * units.cm2ps,
    "carbon_dioxide": 1.92 * 10**-5 * units.cm2ps,
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