import sys
import matplotlib.pyplot as plt
import datetime as dt
import numpy as np


def read_file_to_array(file_path):
    with open(file_path, "r") as file:
        # Read lines into an array, removing newline characters
        lines = [line.rstrip("\n") for line in file]
    return lines


def plot_bars(frequencies, bin_edges):
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    plt.bar(
        bin_centers, frequencies, width=bin_widths, align="center", edgecolor="black"
    )


if len(sys.argv) != 4:
    print(
        "Usage: python3 plot_bars.py [hgmt_debug_output.txt] [x_axis_name] [output location]"
    )
    sys.exit()
lines = read_file_to_array(sys.argv[1])[3:]
frequencies = np.array([float(x.split(": ")[1]) for x in lines])
bin_edges = np.array(
    [float(x.split("-")[0]) for x in lines]
    + [float(lines[-1].split("-")[1].split(": ")[0])]
)
plot_bars(frequencies, bin_edges)
plt.xlabel(sys.argv[2])
plt.ylabel("Frequency")
current_date = dt.datetime.now().strftime("%Y-%m-%d")
plt.text(
    0.98,
    0.98,
    f"{current_date}",
    transform=plt.gca().transAxes,
    fontsize=10,
    verticalalignment="top",
    horizontalalignment="right",
)
plt.savefig(sys.argv[3])
plt.show()
