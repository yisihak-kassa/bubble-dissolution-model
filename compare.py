import pandas as pd
import numpy as np

from bubble_tracker import track_bubble
from bubble_tracker import Q_

short_names = {
    "oxygen": "O",
    "nitrogen": "N",
    "argon": "Ar",
    "helium": "H",
    "methane": "M",
}

GASES = ["methane", "oxygen", "nitrogen", "argon", "carbon_dioxide", "helium"]
ATMOSPHERIC_MOLE_FRACTIONS = [0.0000019, 0.2095, 0.7808, 0.0093, 0.00041, 5.24e-6]
DATA = pd.read_csv("data/grouped_analyzed.csv")
GROUPED = DATA.groupby(["id"])
GROUP_NAMES = GROUPED.groups.keys()


def calculate_fractions(tracer_content, main_gas="methane", tracer="oxygen"):

    if main_gas == "oxygen":

        return [0.0, tracer_content, 1 - tracer_content, 0.0, 0.0, 0.0]

    fractions = dict(zip(GASES, ATMOSPHERIC_MOLE_FRACTIONS))

    total = 1 / fractions[tracer] * tracer_content

    initial_fractions = {gas: total * fractions[gas] for gas in fractions}

    if main_gas == "nitrogen":
        initial_fractions["nitrogen"] = 1 - tracer_content
    else:
        initial_fractions[main_gas] = 1 - total

    return list(initial_fractions.values())


def get_groups(main_gas="methane"):
    groups = []
    for group in GROUP_NAMES:
        gas, initial_diameter = group.split()
        if gas == short_names[main_gas]:
            groups.append((main_gas, float(initial_diameter)))
    return groups


def get_experimental_set(main_gas, initial_diameter):

    group_id = f"{short_names[main_gas.lower()]} {initial_diameter:.2f}"

    experimental_set = DATA[DATA["id"] == group_id]

    if experimental_set.empty:

        raise ValueError(
            "No bubble with these properties found, consider changing gas type or checking the dataset."
        )

    return experimental_set


def compare(
    type,
    initial_diameter,
    initial_depth=4.69,
    tracer="oxygen",
    main_gas="methane",
    temperature_profile=19.0,
    atmospheric_pressure=1.01325,
):

    experimental_data = get_experimental_set(main_gas, initial_diameter)
    initial_oxygen_content = experimental_data[
        experimental_data["travel_distance"]
        == min(experimental_data["travel_distance"])
    ]["adjusted_oxygen_content"].values[0]

    observed_travel_distance = experimental_data["travel_distance"]
    observed_oxygen_content = experimental_data["adjusted_oxygen_content"]

    observed_size = np.array(experimental_data["bubble_diameter"])
    observed_n = experimental_data["n"]

    # Compute model prediction
    fractions = calculate_fractions(initial_oxygen_content / 100, main_gas, tracer)

    history = track_bubble(
        initial_diameter=Q_(initial_diameter, "mm"),
        release_depth=Q_(initial_depth, "m"),
        initial_composition_by_gas=fractions,
        temperature_profile=Q_(temperature_profile, "°C"),
        water_density=Q_(1000, "kg/m³"),
        water_viscosity=Q_(0.001, "Pa·s"),
    )

    modeled_travel_distance = [depth for depth in history["depths"]]

    modeled_size = [diameter for diameter in history["diameters"]]

    modeled_oxygen_content = [
        composition[1] * 100 for composition in history["compositions"]
    ]
    modeled_n = history["bubble_masses"]

    if type == "size":
        observation = (observed_travel_distance, observed_size)
        model_prediction = (modeled_travel_distance, modeled_size)
    elif type == "oxygen":
        observation = (observed_travel_distance, observed_oxygen_content)
        model_prediction = (modeled_travel_distance, modeled_oxygen_content)
    else:
        observation = (observed_travel_distance, observed_n)
        model_prediction = (modeled_travel_distance, modeled_n)

    return observation, model_prediction


def compare_sizes(type, main_gas="methane"):

    groups = get_groups(main_gas)
    comparisons = []
    for group in groups:
        main_gas, initial_diameter = group
        observation, model_prediction = compare(
            type, initial_diameter, main_gas=main_gas
        )
        comparisons.append((main_gas, initial_diameter, observation, model_prediction))
    return type, comparisons
