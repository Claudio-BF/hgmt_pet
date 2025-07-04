import numpy as np
import sys
from vispy import scene
from vispy.scene.visuals import Mesh, Line
from vispy.color import Color
from vispy.scene import visuals
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

script_dir = os.path.dirname(os.path.abspath(__file__))
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
        ex, ey, ez, energy = map(float, (ex, ey, ez, energy))
        e = event(vec3d(ex, ey, ez), energy, detected)
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

    bottom_center = len(vertices)
    top_center = bottom_center + 1
    vertices.append([0, 0, -length / 2])
    vertices.append([0, 0, length / 2])

    for i in range(num_segments):
        i0 = 2 * i
        i1 = 2 * ((i + 1) % num_segments)
        faces += [
            [i0, i1, i1 + 1],
            [i0, i1 + 1, i0 + 1],
            [bottom_center, i0, i1],
            [top_center, i1 + 1, i0 + 1],
        ]
        edges += [[i0 + 1, i1 + 1], [i0, i1], [i0, i0 + 1]]

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
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)
        vertices.append([x_inner, y_inner, -length / 2])
        vertices.append([x_outer, y_outer, -length / 2])
        vertices.append([x_inner, y_inner, length / 2])
        vertices.append([x_outer, y_outer, length / 2])

    for i in range(num_segments):
        i0 = 4 * i
        i1 = 4 * ((i + 1) % num_segments)
        faces += [
            [i0, i1, i1 + 2],
            [i0, i0 + 2, i1 + 2],
            [i0 + 1, i1 + 1, i1 + 3],
            [i0 + 1, i0 + 3, i1 + 3],
            [i0, i1, i0 + 1],
            [i1, i1 + 1, i0 + 1],
            [i0 + 2, i1 + 2, i0 + 3],
            [i1 + 2, i1 + 3, i0 + 3],
        ]
        edges += [
            [i0, i0 + 2],
            [i0 + 3, i0 + 1],
            [i0, i1],
            [i0 + 1, i1 + 1],
            [i0 + 2, i1 + 2],
            [i0 + 3, i1 + 3],
        ]

    return np.array(vertices), np.array(faces), np.array(edges)


# Create a canvas and view
canvas = scene.SceneCanvas(keys="interactive", bgcolor="white")
view = canvas.central_widget.add_view()
view.camera = scene.TurntableCamera(elevation=30, azimuth=30, distance=10000, fov=0.0)


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


def draw_sphere(center, color, radius=15):
    marker = scene.visuals.Markers(
        pos=np.array([center]), size=radius, face_color=color, edge_color=None
    )
    marker.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
    view.add(marker)
    return [marker]


def label(position, text):
    text_label = visuals.Text(
        text=text,
        pos=position,
        color="black",
        font_size=20,
        anchor_x="center",
        anchor_y="bottom",
    )
    view.add(text_label)
    return [text_label]


def draw_path(origin, events, pathid):
    visuals_created = []
    points = [[origin.x, origin.y, origin.z]] + [
        [e.position.x, e.position.y, e.position.z] for e in events
    ]
    path = Line(pos=points, color="red", width=1)
    path.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
    view.add(path)
    visuals_created.append(path)
    for i in range(len(events)):
        e = events[i]
        visuals_created += draw_sphere(
            [e.position.x, e.position.y, e.position.z],
            Color("green" if e.detected else "grey"),
        )
        visuals_created += label(points[i + 1], pathid + str(i + 1))
        new_pathid = pathid.upper() if e.detected else pathid
        print(new_pathid + str(i + 1) + ": " + str(e.energy))
    return visuals_created


def draw_annihilation(origin, events1, events2):
    visuals_created = []
    visuals_created += draw_path(origin, events1, "a")
    print("")
    visuals_created += draw_path(origin, events2, "b")
    visuals_created += draw_sphere([origin.x, origin.y, origin.z], Color("red"))
    return visuals_created


# --- Navigation logic starts here ---

annihilations = read_file(data_dir)
current_index = argument if argument < len(annihilations) else 0
visuals_list = []


def redraw():
    # Remove previous visuals
    for v in visuals_list:
        v.parent = None
    visuals_list.clear()

    # Draw new visuals and track them
    for inner_radius in detector_inner_radii:
        visuals_list.extend(
            draw_tube(
                inner_radius, inner_radius + detector_thickness, detector_length, 30
            )
        )
    visuals_list.extend(draw_cylinder(10.6, 4, 30))
    origin, events1, events2 = annihilations[current_index]
    visuals_list.extend(draw_annihilation(origin, events1, events2))


redraw()


@canvas.events.key_press.connect
def on_key_press(event):
    global current_index
    print("\n")
    if event.key == "j":
        current_index = (current_index - 1) % len(annihilations)
        print("Visualizing event: " + str(current_index) + "\n")
        redraw()
    elif event.key == "k":
        current_index = (current_index + 1) % len(annihilations)
        print("Visualizing event: " + str(current_index) + "\n")
        redraw()


canvas.show()

if __name__ == "__main__":
    canvas.app.run()
