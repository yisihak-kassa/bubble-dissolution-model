"""
This module provides functionality to analyze sample data from an Excel file
containing information about exetainers and samples. It calculates various
parameters such as bubble volume, bubble diameter, adjusted oxygen content,
travel distance, and capture depth.
Functions:
- load_data(file_name): Load data from an Excel file into DataFrames.
- calculate_parameters(samples_df, exetainers_df): Perform calculations on
  the sample data to derive various parameters.
- main(): Main function to load data, process it, and save the output to a CSV file.
Constants:
- g: Gravitational constant (m/s^2).
- R: Universal gas constant (L*bar/(mol*K)).
"""

# Imports
import pandas as pd
from data_utils import (
    calculate_diameter,
    calculate_depth_adjusted_volume,
    adjust_oxygen_content,
    calculate_water_density,
    calculate_travel_distance,
    calculate_capture_depth,
)

g = 9.81  # m/s^2
R = 0.0831446261815324  # "L*bar/(mol*K)


# Define functions
def load_data(file_name):
    """
    Load data from an Excel file into DataFrames.

    Parameters:
    - file_name (str): The name of the Excel file to read.

    Returns:
    - Tuple of DataFrames: (samples_df, exetainers_df)
    """
    samples_df = pd.read_excel(file_name, sheet_name="samples")
    exetainers_df = pd.read_excel(file_name, sheet_name="exetainers")
    exetainers_df.set_index("label", inplace=True)
    return samples_df, exetainers_df


def calculate_parameters(samples_df, exetainers_df):
    """
    Perform calculations on the sample data to derive parameters such as
    bubble volume, bubble diameter, adjusted oxygen content, travel distance,
    and capture depth.

    Parameters:
    - samples_df (DataFrame): DataFrame containing sample data.
    - exetainers_df (DataFrame): DataFrame containing exetainer data.

    Returns:
    - DataFrame: A DataFrame containing calculated parameters.
    """
    analyzed_data_df = pd.DataFrame()
    # no = samples_df["no"]
    exetainer = samples_df["exetainer"]
    sampling_depth = samples_df["sampling_depth"].fillna(5)
    atm_pressure = samples_df["atm_pressure"].fillna(1013.25)
    oxygen_content = samples_df["oxygen_content"]
    oxygen_test_temperature = samples_df["test_T"].fillna(23.5)

    full_mass = (
        samples_df["exetainer"]
        .apply(lambda x: exetainers_df.loc[x, "full_mass"])
        .fillna(24.97)
    )

    water_density = samples_df["water_T"].fillna(19).apply(calculate_water_density)

    gas_volume = (full_mass - samples_df["sample_mass"]) / water_density

    bubble_count = (
        samples_df["count"]
        / samples_df["counting_duration"]
        * samples_df["sampling_duration"]
    )

    bubble_volume = gas_volume / bubble_count

    # Depth adjusted bubble volume in ml
    adjusted_bubble_volume = pd.Series(
        map(
            calculate_depth_adjusted_volume,
            bubble_volume,
            sampling_depth,
            atm_pressure,
            water_density,
        )
    )

    # in mm
    bubble_diameter = 10 * adjusted_bubble_volume.apply(calculate_diameter)

    adjusted_oxygen_content = pd.Series(
        map(
            adjust_oxygen_content,
            oxygen_content,
            gas_volume,
            oxygen_test_temperature,
        )
    )

    # Store calculated parameters in analysis_df
    analyzed_data_df["group"] = samples_df["group"]
    analyzed_data_df["exetainer"] = exetainer
    analyzed_data_df["gas"] = samples_df["gas"]
    analyzed_data_df["oxygen_content"] = oxygen_content.apply(lambda x: f"{x:.2f}")
    analyzed_data_df["adjusted_oxygen_content"] = adjusted_oxygen_content.apply(
        lambda x: f"{x:.2f}"
    )
    analyzed_data_df["bubble_diameter"] = bubble_diameter.apply(lambda x: f"{x:.2f}")
    analyzed_data_df["travel_distance"] = sampling_depth.apply(
        lambda x: f"{(calculate_travel_distance(x)):.2f}"
    )
    analyzed_data_df["capture_depth"] = sampling_depth.apply(
        lambda x: float(calculate_capture_depth(x))
    )

    P = (
        atm_pressure / 1000
        + water_density * g * analyzed_data_df["capture_depth"] * 1e-2
    )
    V = adjusted_bubble_volume * 1e-3

    T = samples_df["water_T"].fillna(19) + 273.15
    n = P * V / (R * T)
    analyzed_data_df["n"] = n
    return analyzed_data_df


def main() -> None:

    input_file = "data/grouped_data.xlsx"

    # Load and process the data
    samples_df, exetainers_df = load_data(input_file)
    analyzed_data_df = calculate_parameters(samples_df, exetainers_df)

    # Save the output to a CSV file
    output_file = "data/grouped_analyzed.csv"
    analyzed_data_df.to_csv(output_file, index=False, header=True)
    print(f"\nOutput written to {output_file}\n")


if __name__ == "__main__":
    main()
