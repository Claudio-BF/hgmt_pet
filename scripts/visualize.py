import numpy as np
import sys
from vispy import scene
from vispy.scene.visuals import Mesh, Line, Text
from vispy.color import Color
from vispy.scene import STTransform, visuals
from vispy.scene import Markers
import vispy.app
from dataclasses import dataclass
from vispy.visuals.transforms import MatrixTransform, ChainTransform
import struct
import os

vispy.app.use_app("pyqt6")

# Parameters
detector_length = 200  # cm
detector_thickness = 2.54  # cm
detector_inner_radii = np.array([45] * 12) + 5 * np.array(range(12))  # MUST BE SORTED
fly_camera = False
argument = 0

# Correctly determine the script and data directories
try:
    script_dir = os.path.dirname(os.path.abspath(__file__))
except NameError:
    script_dir = os.getcwd()  # Fallback for interactive environments
data_dir = os.path.join(script_dir, "../data/visualization.data")
print("Reading data from: " + data_dir)

if len(sys.argv) > 1:
    argument = int(sys.argv[1])
    print(f"Visualizing annihilation: {argument}")
else:
    print("No argument provided, visualizing annihilation: 0")
print("Use '[' and ']' keys to scroll through annihilations.\n\n")


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


@dataclass
class LOR:
    center: vec3d
    transform: tuple[float, ...]
    hit1: vec3d
    hit2: vec3d


def read_events(f):
    n_bytes = f.read(4)
    n = struct.unpack("i", n_bytes)[0]

    events = []
    for _ in range(n):
        ex, ey, ez, energy, detected = struct.unpack("dddd?", f.read(33))
        e = event(vec3d(ex, ey, ez), float(energy), detected)
        events.append(e)
    return events


def read_file(file_path):
    """Reads the entire binary file containing annihilation histories."""
    annihilations = []
    with open(file_path, "rb") as f:
        while True:
            # Read annihilation center (3 doubles = 24 bytes)
            center_bytes = f.read(24)
            if not center_bytes:
                break  # End of file

            vx, vy, vz = struct.unpack("ddd", center_bytes)
            center = vec3d(float(vx), float(vy), float(vz))

            # Read the event paths for both photons
            events1 = read_events(f)
            events2 = read_events(f)

            # Check if a Line of Response (LOR) was generated (1 byte bool)
            made_lor_byte = f.read(1)
            made_lor = struct.unpack("?", made_lor_byte)[0]

            lor_data = None
            if made_lor:
                # Read LOR center (3 doubles = 24 bytes)
                lor_c_x, lor_c_y, lor_c_z = struct.unpack("ddd", f.read(24))
                lor_center = vec3d(lor_c_x, lor_c_y, lor_c_z)

                # Read Cholesky transform matrix (6 doubles = 48 bytes)
                transform_matrix = struct.unpack("dddddd", f.read(48))

                # Read hit positions (3 doubles each = 24 * 2 = 48 bytes)
                h1_x, h1_y, h1_z = struct.unpack("ddd", f.read(24))
                h2_x, h2_y, h2_z = struct.unpack("ddd", f.read(24))
                hit1_pos = vec3d(h1_x, h1_y, h1_z)
                hit2_pos = vec3d(h2_x, h2_y, h2_z)

                lor_data = LOR(lor_center, transform_matrix, hit1_pos, hit2_pos)

            annihilations.append((center, events1, events2, lor_data))

    return annihilations


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
if fly_camera:
    view.camera = scene.FlyCamera(fov=60, scale_factor=1000)
else:
    view.camera = scene.TurntableCamera(
        elevation=30, azimuth=30, fov=0, scale_factor=300
    )
# Initialize visual elements that will be updated
path = Line(pos=np.empty((0, 3)), color="red", width=2, method="gl")
lor_line = Line(pos=np.empty((0, 3)), color="cyan", width=2, method="gl")
labels = Text(
    text=[""],
    pos=[[0, 0, 0]],
    color="black",
    font_size=10,
    anchor_x="center",
    anchor_y="bottom",
)
markers = Markers(pos=np.empty((0, 3)), size=6, edge_color=None)
view.add(path)
view.add(lor_line)
view.add(labels)
lor = scene.visuals.Sphere(radius=1, method="latitude", parent=view.scene, color="blue")
view.add(markers)
# Set GL state for transparency
for item in [path, markers, labels, lor_line, lor]:
    item.order = 1
    item.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )


def draw_lor(lor_data):
    if not lor_data:
        selected_hits, colors = np.empty((0, 3)), np.empty((0, 4))
        # Hide the LOR sphere if there is no data
        lor.visible = False
    else:
        # Ensure the sphere is visible when there is data
        lor.visible = True
        selected_hits = np.array(
            [[h.x, h.y, h.z] for h in [lor_data.hit1, lor_data.hit2]]
        )
        colors = np.array([Color("blue").rgba] * 2)

        # 1. Create a 4x4 matrix for the linear transformation part
        #    (scaling, rotation, shearing) from the Cholesky factors.
        linear_matrix = np.eye(4)
        L = lor_data.transform
        linear_matrix[0, 0] = L[0]
        linear_matrix[1, 0] = L[1]
        linear_matrix[1, 1] = L[2]
        linear_matrix[2, 0] = L[3]
        linear_matrix[2, 1] = L[4]
        linear_matrix[2, 2] = L[5]
        linear_tf = MatrixTransform(linear_matrix.T)
        translation_tf = STTransform(
            translate=(lor_data.center.x, lor_data.center.y, lor_data.center.z)
        )
        lor.transform = ChainTransform([translation_tf, linear_tf])

    lor_line.set_data(pos=selected_hits)
    return selected_hits, colors


def draw_path(events, pathid):
    if not events:
        return np.empty((0, 3)), np.empty((0, 4)), np.empty(0, dtype=str)

    points = np.array([[e.position.x, e.position.y, e.position.z] for e in events])
    colors = np.array([Color("green" if e.detected else "grey") for e in events])
    labels = np.array([f"{pathid}{i + 1}" for i in range(len(events))])
    for i in range(len(events)):
        print(
            (labels[i].upper() if events[i].detected else labels[i].lower())
            + ": "
            + str(events[i].energy)
        )
    return points, colors, labels


def draw_annihilation(origin, events1, events2, lor_data):
    center_pos = np.array([[origin.x, origin.y, origin.z]])

    points1, colors1, labels1 = draw_path(events1, "a")
    print("")
    points2, colors2, labels2 = draw_path(events2, "b")
    selected_hits, hit_colors = draw_lor(lor_data)

    full_path_pos = np.concatenate([np.flip(points1, axis=0), center_pos, points2])
    path.set_data(pos=full_path_pos, color="purple")
    lor_line.set_data(pos=selected_hits)

    all_points_pos = np.concatenate([points1, center_pos, points2])
    all_labels_text = np.concatenate([labels1, [""], labels2])
    labels.pos = all_points_pos
    labels.text = all_labels_text
    all_points_pos = np.concatenate([all_points_pos, selected_hits])
    center_color = np.array([Color("red")])
    all_marker_colors = np.concatenate([colors1, center_color, colors2, hit_colors])

    markers.set_data(
        pos=all_points_pos, face_color=all_marker_colors, size=6, edge_width=0
    )


def redraw():
    origin, events1, events2, lor_data = annihilations[current_index]
    draw_annihilation(origin, events1, events2, lor_data)


annihilations = read_file(data_dir)
current_index = argument if argument < len(annihilations) else 0

for inner_radius in detector_inner_radii:
    draw_tube(inner_radius, inner_radius + detector_thickness, detector_length, 60)
draw_cylinder(10.6, 4, 30)
if fly_camera:
    view.camera.scale_factor = 50
redraw()


@canvas.events.key_press.connect
def on_key_press(event):
    global current_index
    if event.key == "[" and current_index > 0:
        current_index -= 1
    elif event.key == "]" and current_index < len(annihilations) - 1:
        current_index += 1
    else:
        return

    print(f"\n\nVisualizing annihilation: {current_index}\n")
    redraw()


canvas.show()

if __name__ == "__main__":
    if sys.flags.interactive != 1:
        vispy.app.run()
