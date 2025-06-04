import sys
import numpy as np
import matplotlib.pyplot as plt


def render_z_slice(file_path, z_index):
    with open(file_path, "rb") as f:
        x_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        y_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        z_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        x_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        y_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        z_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        pixels = np.fromfile(f, dtype=np.float64, count=x_res * y_res * z_res)
    slice_2d = pixels.reshape((x_res, y_res, z_res), order="C")[:, :, z_index]
    plt.imshow(
        slice_2d.T,
        origin="lower",
        extent=[-x_len / 2, x_len / 2, -y_len / 2, y_len / 2],
        cmap="grey",
    )
    plt.colorbar(label="Intensity")
    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python render_image.py [image.pixels] [z-index]")
        sys.exit(1)
    output_dir = sys.argv[1]
    file_path = f"{sys.argv[1]}"
    render_z_slice(file_path, int(sys.argv[2]))
