import struct
import matplotlib.pyplot as plt
import numpy as np
import datetime as dt
import sys


def read_doubles_from_binary_file(filename):
    with open(filename, "rb") as f:
        data = f.read()
    return [d[0] for d in struct.iter_unpack("d", data)]


def plot_histogram(doubles, xmax):
    print("plotting histogram with " + str(len(doubles)) + " data points")
    counts, bin_edges = np.histogram(doubles, bins=50, range=(0, xmax))
    # Compute bin centers
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_widths = bin_edges[1:] - bin_edges[:-1]
    # normalize the data
    normalized = [float(counts[i]) / bin_widths[i] for i in range(len(counts))]
    normalized = counts / bin_widths
    normalized /= len(doubles)
    # Plot histogram as line graph using matplotlib
    plt.plot(bin_centers, normalized)


# replace 'file.bin' with your binary file
if len(sys.argv) != 6:
    print(str(len(sys.argv)) + "parameters, expected 4, usage:")
    print(
        "python3 plot_histogram [data_loc] [x_axis_label] [output_location] [x_axis_max] [y_axis_max]"
    )
    sys.exit()
doubles = read_doubles_from_binary_file(sys.argv[1])
plt.xlabel(sys.argv[2])
plt.ylabel("Frequency")
plt.xlim(0, float(sys.argv[4]))
plt.ylim(0, float(sys.argv[5]))
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
plt.gcf().canvas.get_default_filename = lambda: sys.argv[3]
plot_histogram(doubles, float(sys.argv[4]))
font = {"family": "normal", "weight": "bold", "size": 22}

plt.rc("font", **font)
plt.savefig(sys.argv[3])
plt.show()
