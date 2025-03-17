"""
This module provides functionality to plot comparisons of the bubble dissolution model
against observed data. It includes custom legend handlers and a function to generate
comparison plots with error calculations.
Classes:
    HandlerRect: Custom handler to replace line markers with rectangular patches in legend.
Functions:
    plot_comparisons(data):
        Plots comparisons of observed and model-predicted data for bubble dissolution.
        Args:
            data (tuple): A tuple containing the type of comparison and a list of comparisons.
                          Each comparison is a tuple of (main_gas, diameter, observation, model_prediction).
        Returns:
            fig (matplotlib.figure.Figure): The generated plot figure.
            errors (list): List of relative errors for each comparison.
"""

from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import interp1d
from matplotlib.patches import Rectangle
from matplotlib.legend_handler import HandlerPatch

from compare import compare_sizes

# Y-axis lables depending on the type of comparison being made ("size", "oxygen", "n")
y_axis_labels = {"size": "Size (mm)", "oxygen": "Oxygen Content (%)", "n": "n (μmol)"}


# Custom handler to replace line markers with rectangular patches in legend
class HandlerRect(HandlerPatch):
    def create_artists(
        self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans
    ):
        p = Rectangle(
            [xdescent, ydescent],
            width,
            height,
            facecolor=orig_handle.get_facecolor(),
        )
        return [p]


def plot_comparisons(data):
    # Color-blind friendly palette
    color_palette = [
        "#808080",  # Gray
        "#a559aa",  # Purple
        "#298c8c",  # Teal
        "#f0c571",  # Gold
        "#a1a1a1",  # Gray
        "#082a54",  # Dark Blue
    ]

    fig, ax = plt.subplots(figsize=(9 / 2.54, 7 / 2.54), dpi=600)
    plt.rcParams["font.family"] = "Arial"

    for spine in ax.spines.values():
        spine.set_linewidth(0.6)

    # Prepare legend elements
    handles, labels = [], []
    type, comparisons = data

    errors = []

    for comparison in comparisons:
        main_gas, diameter, observation, model_prediction = comparison
        color = color_palette.pop(0)

        ax.scatter(
            observation[0],
            observation[1],
            marker="o",
            s=8,
            label=f"{diameter} mm",
            color=color,
        )
        ax.plot(
            model_prediction[0],
            model_prediction[1],
            linestyle="--",
            color=color,
            linewidth=1.2,
        )

        interpolator = interp1d(
            model_prediction[0],
            model_prediction[1],
            kind="cubic",
            fill_value="extrapolate",
        )
        modeled_content_at_observed_distances = interpolator(observation[0])
        observed_content = observation[1]

        relative_error = np.sqrt(
            np.mean((observed_content - modeled_content_at_observed_distances) ** 2)
        )

        errors.append(relative_error)
        print(f"{main_gas.capitalize()} {diameter} mm: {relative_error:.2f}")

        handles.append(
            Rectangle((0, 0), 0.5, 0.5, facecolor=color, edgecolor="black", lw=0.5)
        )
        labels.append(f"{diameter} {relative_error:.2f}")

    ax.tick_params(
        axis="both",
        which="both",
        direction="out",
        length=5,
        width=0.4,
        colors="black",
        top=False,
        bottom=True,
        left=True,
        right=False,
        labelsize=8,
    )

    if main_gas == "oxygen":
        ax.tick_params(axis="y", which="minor", length=2.5, width=0.25, labelsize=8)
        ax.yaxis.set_minor_locator(plt.MultipleLocator(1))

    ax.grid(
        visible=True,
        which="major",
        linestyle="-",
        linewidth=0.5,
        color="gray",
        alpha=0.5,
    )

    legend1 = ax.legend(
        handles,
        labels,
        title="$d_e$ | RMSE",
        title_fontsize=8,
        fontsize=8,
        loc="upper left",
        ncol=1,
        bbox_to_anchor=(1, 1),
        frameon=False,
        handletextpad=0.5,
        handler_map={Rectangle: HandlerRect()},
    )
    ax.add_artist(legend1)
    handles2 = [
        Line2D([0], [0], color="black", lw=0, marker="o", ms=2, label="Observed"),
        Line2D([0], [0], color="black", lw=0.8, linestyle="--", label="Predicted"),
    ]

    ax.legend(
        handles=handles2,
        loc="upper left",
        bbox_to_anchor=(0.1, 0.8),
        title="",
        title_fontsize=8,
        fontsize=8,
        frameon=True,
    )

    ax.yaxis.get_offset_text().set_visible(False)

    ax.set_title(
        f"{main_gas.capitalize()}",
        fontweight="bold",
        color="#444444",
        fontsize=10,
    )

    ax.set_xlabel("Travel Distance (m)", fontsize=10, fontweight="bold")
    ax.set_ylabel(f"{y_axis_labels[type]}", fontsize=10, fontweight="bold")

    plt.tight_layout()
    fig.subplots_adjust(right=0.75)
    return fig, errors


exp_gases = ["methane", "nitrogen", "argon", "helium"]
figs = []
errors = []
for gas in exp_gases:
    figure, error = plot_comparisons(compare_sizes("oxygen", gas))
    figs.append(figure)
    errors += error

errors = np.array(errors)
print(f"mean RMSE: {np.mean(errors):.2f}")
print(f"std RMSE: {np.std(errors):.2f}")
