import numpy as np
import sys
from vispy import scene
from vispy.scene.visuals import Mesh, Line, Text
from vispy.color import Color
from vispy.scene import visuals
from vispy.scene import Markers
import vispy.app
from dataclasses import dataclass
import struct
import os

vispy.app.use_app("pyqt6")

# Parameters
detector_length = 200  # cm
detector_thickness = 2.54  # cm
detector_inner_radii = np.array([45] * 12) + 5 * np.array(range(12))  # MUST BE SORTED
argument = 0

# Correctly determine the script and data directories
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd()  # Fallback for interactive environments
data_dir = os.path.join(script_dir, "../data/visualization.data")
print("Reading data from: " + data_dir)


@dataclass
class vec3d:
    x: float
    y: float
    z: float


@dataclass
class event:
    position: vec3d
    energy: float
    detected: bool


if len(sys.argv) > 1:
    argument = int(sys.argv[1])
    print(f"Visualizing event: {argument}")
else:
    print("No argument provided, visualizing event: 0")
print("Use j and k to scroll through annihilations.\n\n")


def read_events(f):
    n = struct.unpack("i", f.read(4))[0]
    events = []
    for _ in range(n):
        ex, ey, ez, energy, detected = struct.unpack("dddd?", f.read(33))
        e = event(vec3d(ex, ey, ez), float(energy), detected)
        events.append(e)
    return n, events


def read_file(file):
    annihils = []
    with open(file, "rb") as f:
        while True:
            vec_bytes = f.read(24)  # 3 doubles = 24 bytes
            if not vec_bytes:
                break
            vx, vy, vz = struct.unpack("ddd", vec_bytes)
            center = vec3d(float(vx), float(vy), float(vz))

            n1, events1 = read_events(f)
            n2, events2 = read_events(f)

            annihils.append((center, events1, events2))
    return annihils


def create_cylinder_mesh(radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)
        vertices.append([x, y, -length / 2])
        vertices.append([x, y, length / 2])

    bottom_center_idx = len(vertices)
    vertices.append([0, 0, -length / 2])
    top_center_idx = len(vertices)
    vertices.append([0, 0, length / 2])

    for i in range(num_segments):
        i0 = 2 * i
        i1 = 2 * ((i + 1) % num_segments)
        faces.extend(
            [
                [i0, i1, i1 + 1],
                [i0, i1 + 1, i0 + 1],
                [bottom_center_idx, i0, i1],
                [top_center_idx, i1 + 1, i0 + 1],
            ]
        )
        edges.extend([[i0 + 1, i1 + 1], [i0, i1], [i0, i0 + 1]])

    return np.array(vertices), np.array(faces), np.array(edges)


def draw_cylinder(radius, length, num_segments):
    vertices, faces, edges = create_cylinder_mesh(radius, length, num_segments)
    cylinder_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightblue"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(cylinder_mesh)
    return [wireframe, cylinder_mesh]


def create_tube_mesh(inner_radius, outer_radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x_inner, y_inner = inner_radius * np.cos(theta), inner_radius * np.sin(theta)
        x_outer, y_outer = outer_radius * np.cos(theta), outer_radius * np.sin(theta)
        vertices.extend(
            [
                [x_inner, y_inner, -length / 2],
                [x_outer, y_outer, -length / 2],
                [x_inner, y_inner, length / 2],
                [x_outer, y_outer, length / 2],
            ]
        )

    for i in range(num_segments):
        i0 = 4 * i
        i1 = 4 * ((i + 1) % num_segments)
        faces.extend(
            [
                [i0, i1, i1 + 2],
                [i0, i0 + 2, i1 + 2],
                [i0 + 1, i1 + 1, i1 + 3],
                [i0 + 1, i0 + 3, i1 + 3],
                [i0, i1, i0 + 1],
                [i1, i1 + 1, i0 + 1],
                [i0 + 2, i1 + 2, i0 + 3],
                [i1 + 2, i1 + 3, i0 + 3],
            ]
        )
        edges.extend(
            [
                [i0, i0 + 2],
                [i0 + 3, i0 + 1],
                [i0, i1],
                [i0 + 1, i1 + 1],
                [i0 + 2, i1 + 2],
                [i0 + 3, i1 + 3],
            ]
        )

    return np.array(vertices), np.array(faces), np.array(edges)


def draw_tube(inner_radius, outer_radius, length, num_segments):
    vertices, faces, edges = create_tube_mesh(
        inner_radius, outer_radius, length, num_segments
    )
    tube_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightyellow"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(tube_mesh)
    return [wireframe, tube_mesh]


# --- Vispy setup ---
canvas = scene.SceneCanvas(keys="interactive", bgcolor="white")
view = canvas.central_widget.add_view()
view.camera = scene.TurntableCamera(elevation=30, azimuth=30, fov=0, scale_factor=300)
# Initialize visual elements that will be updated
path = Line(pos=np.array([[0, 0, 0]]), color="red", width=2, method="gl")
labels = Text(
    text=[""],
    pos=[[0, 0, 0]],
    color="black",
    font_size=15,
    anchor_x="center",
    anchor_y="bottom",
)
markers = Markers(pos=np.array([[0, 0, 0]]), size=8, edge_color=None)

# Set GL state for transparency
for item in [path, markers, labels]:
    item.order = 1
    item.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
view.add(path)
view.add(labels)
view.add(markers)


def draw_path(events, pathid):
    if not events:
        return np.empty((0, 3)), np.empty((0, 4)), np.empty(0, dtype=str)

    points = np.array([[e.position.x, e.position.y, e.position.z] for e in events])
    colors = np.array([Color("green" if e.detected else "grey") for e in events])
    labels = np.array([f"{pathid}{i + 1}" for i in range(len(events))])
    for i in range(len(events)):
        print(labels[i] + ": " + str(events[i].energy))
    return points, colors, labels


def draw_annihilation(origin, events1, events2):
    center_pos = np.array([[origin.x, origin.y, origin.z]])

    points1, colors1, labels1 = draw_path(events1, "a")
    print("")
    points2, colors2, labels2 = draw_path(events2, "b")

    full_path_pos = np.concatenate([np.flip(points1, axis=0), center_pos, points2])
    path.set_data(pos=full_path_pos, color="purple")

    all_points_pos = np.concatenate([points1, center_pos, points2])
    all_labels_text = np.concatenate([labels1, [""], labels2])
    center_color = np.array([Color("red")])
    all_marker_colors = np.concatenate([colors1, center_color, colors2])

    labels.pos = all_points_pos
    labels.text = all_labels_text
    markers.set_data(pos=all_points_pos, face_color=all_marker_colors)


def redraw():
    origin, events1, events2 = annihilations[current_index]
    draw_annihilation(origin, events1, events2)


annihilations = read_file(data_dir)
current_index = argument if argument < len(annihilations) else 0

for inner_radius in detector_inner_radii:
    draw_tube(inner_radius, inner_radius + detector_thickness, detector_length, 60)
draw_cylinder(10.6, 4, 30)
redraw()


@canvas.events.key_press.connect
def on_key_press(event):
    global current_index
    if event.key == "j":
        current_index = (current_index - 1) % len(annihilations)
    elif event.key == "k":
        current_index = (current_index + 1) % len(annihilations)
    else:
        return

    print(f"\n\nVisualizing event: {current_index}\n")
    redraw()


canvas.show()

if __name__ == "__main__":
    if sys.flags.interactive != 1:
        vispy.app.run()
