import pint

units = pint.UnitRegistry()

m = units.meter
cm = units.cm
mm = units.mm
s = units.second
m3 = m**3
cm3 = cm**3
mm3 = mm**3
liter = units.liter
ml = cm3
m2 = m**2
cm2 = cm**2
mm2 = mm**2
kg = units.kilogram
g = units.gram
N = units.newton
mol = units.mol

# speed
mps = m / s
cmps = cm / s

# acceleration
mps2 = m / s**2
# temperature
C = units.celsius
K = units.kelvin

# pressure
Pa = units.pascal
bar = units.bar
millibar = units.millibar
# molar mass
gpmol = g / mol
# density
kgpm3 = kg / m**3
gpcm3 = g / cm**3
gpml = gpcm3

# molar concentration
mpl = mol/liter
mpcm3 = mol / cm3
molpm3bar = mol / m3 / bar
molpcm3bar = mol / cm3 / bar

# molar volume
cm3pmol = cm3 / mol
# viscosity
P = units.poise
cP = units.centipoise
Pa_s = Pa * s
mPa_s = Pa * s / 1000
kgpms = kg / m / s

# surface tension
Npm = N / m

# Henry's law constant
molplbar = mol / liter / bar

# ideal gas constant
# kg_m2_per_s2_mol_K = kg * m2 / s**2 / mol / K
LbarpmolK = liter * bar / mol / K

# Avogadro's number
Avogadro = 6.02214076e23 / mol

#


# units.load_definitions("unit_definitions.txt")
def in_celsius(temp):
    Q_ = units.Quantity
    return Q_(temp, units.degC)


# diffusivity
cm2ps = cm2 / s
v = 1.0 * cm2ps
v.units
