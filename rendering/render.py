import sys
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python render_image.py [image.pixels]")
        sys.exit(1)
    file_path = sys.argv[1]
    with open(file_path, "rb") as f:
        x_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        y_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        z_len = np.fromfile(f, dtype=np.float64, count=1)[0]
        iterations = np.fromfile(f, dtype=np.int32, count=1)[0]
        x_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        y_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        z_res = np.fromfile(f, dtype=np.int32, count=1)[0]
        images = np.fromfile(
            f, dtype=np.float64, count=iterations * x_res * y_res * z_res
        )
    images = images.reshape((iterations, x_res, y_res, z_res))

    idx_a = 0
    idx_b = 0

    fig, ax = plt.subplots()
    plt.subplots_adjust(bottom=0.25)
    img = ax.imshow(
        images[idx_a, :, :, idx_b],
        cmap="gray",
        extent=[-x_len / 2, x_len / 2, -y_len / 2, y_len / 2],
        origin="lower",
    )
    cbar = fig.colorbar(img, ax=ax, label="Intensity")

    ax_a = plt.axes([0.15, 0.1, 0.65, 0.03])
    ax_b = plt.axes([0.15, 0.05, 0.65, 0.03])

    slider_a = Slider(ax_a, "iterations", 0, iterations - 1, valinit=idx_a, valstep=1)
    slider_b = Slider(ax_b, "z index", 0, z_res - 1, valinit=idx_b, valstep=1)

    def update(val):
        idx_a = int(slider_a.val)
        idx_b = int(slider_b.val)
        img.set_data(images[idx_a, :, :, idx_b])
        fig.canvas.draw_idle()

    slider_a.on_changed(update)
    slider_b.on_changed(update)

    plt.xlabel("X")
    plt.ylabel("Y")
    plt.show()
