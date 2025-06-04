import numpy as np
import six
import sys
from vispy import scene
from vispy.scene.visuals import Mesh, Line
from vispy.color import Color
from vispy.scene import visuals
import vispy.app

vispy.app.use_app("pyqt6")

# Parameters
detector_length = 200  # cm
detector_thickness = 2.54  # cm
detector_inner_radii = np.array([45] * 12) + 5 * np.array(range(12))  # MUST BE SORTED
argument = 0

if len(sys.argv) > 1:
    # The first argument is at index 1
    argument = int(sys.argv[1])
    print(f"visualizing event {argument}")
else:
    print("no argument provided, visualizing event 0")


def read_file(filename):
    with open(filename, "r") as file:
        file_content = file.read()
    groups = file_content.strip().split("\n\n\n")[argument].split("\n\n")

    # Initialize the result lists
    center = tuple(map(float, groups[0].split()))
    locs1, energies1, detected1 = map(
        list,
        zip(
            *[
                (
                    tuple(map(float, item.split()[:3])),
                    float(item.split()[3]),
                    int(item.split()[4]),
                )
                for item in groups[1].strip().split("\n")
            ]
        ),
    )
    locs2, energies2, detected2 = map(
        list,
        zip(
            *[
                (
                    tuple(map(float, item.split()[:3])),
                    float(item.split()[3]),
                    int(item.split()[4]),
                )
                for item in groups[2].strip().split("\n")
            ]
        ),
    )
    return center, locs1, locs2, energies1, energies2, detected1, detected2


def create_cylinder_mesh(radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    # Generate vertices for the cylinder
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x = radius * np.cos(theta)
        y = radius * np.sin(theta)

        # Bottom vertex
        vertices.append([x, y, -length / 2])
        # Top vertex
        vertices.append([x, y, length / 2])

    # Generate faces for the bottom and top caps
    bottom_center = len(vertices)
    top_center = bottom_center + 1
    vertices.append([0, 0, -length / 2])  # Bottom center
    vertices.append([0, 0, length / 2])  # Top center
    # Generate faces for the sides
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
    # Create the cylinder mesh
    vertices, faces, edges = create_cylinder_mesh(radius, length, num_segments)

    # Create the solid cylinder mesh visual
    cylinder_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightblue"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(cylinder_mesh)


# Create tube mesh
def create_tube_mesh(inner_radius, outer_radius, length, num_segments):
    vertices = []
    faces = []
    edges = []

    # Generate vertices for inner and outer surfaces
    for i in range(num_segments):
        theta = 2 * np.pi * i / num_segments
        x_inner = inner_radius * np.cos(theta)
        y_inner = inner_radius * np.sin(theta)
        x_outer = outer_radius * np.cos(theta)
        y_outer = outer_radius * np.sin(theta)

        # Bottom vertices
        vertices.append([x_inner, y_inner, -length / 2])
        vertices.append([x_outer, y_outer, -length / 2])
        # Top vertices
        vertices.append([x_inner, y_inner, length / 2])
        vertices.append([x_outer, y_outer, length / 2])

    # Generate faces for the sides
    for i in range(num_segments):
        i0 = 4 * i
        i1 = 4 * ((i + 1) % num_segments)
        # Faces
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
        # Edges
        edges += [
            [i0, i0 + 2],
            [i0 + 3, i0 + 1],
            [i0, i1],
            [i0 + 1, i1 + 1],
            [i0 + 2, i1 + 2],
            [i0 + 3, i1 + 3],
        ]

    return np.array(vertices), np.array(faces), np.array(edges)


# Create a canvas
# Create a canvas and view
canvas = scene.SceneCanvas(keys="interactive", bgcolor="white")
view = canvas.central_widget.add_view()
# Set up the camera
view.camera = scene.TurntableCamera(elevation=30, azimuth=30, distance=10000, fov=0.0)
# Add text to the overlay


def draw_tube(inner_radius, outer_radius, length, num_segments):
    # Create the tube mesh
    vertices, faces, edges = create_tube_mesh(
        inner_radius, outer_radius, length, num_segments
    )

    # Create the solid tube mesh visual
    tube_mesh = Mesh(vertices=vertices, faces=faces, color=Color("lightyellow"))
    wireframe = visuals.Line(
        pos=vertices, connect=edges, color=Color("black"), method="gl"
    )
    view.add(wireframe)
    view.add(tube_mesh)


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


def label(position, text):
    text_label = visuals.Text(
        text=text,  # Text to display
        pos=position,  # Position of the text (same as the marker)
        color="black",  # Text color
        font_size=20,  # Font size
        anchor_x="center",  # Center the text horizontally
        anchor_y="bottom",  # Position the text above the marker
    )
    view.add(text_label)


def draw_path(points, detected, pathid, energies):
    path = Line(pos=points, color="red", width=1)
    path.set_gl_state(
        "translucent",
        depth_test=False,
        blend=True,
        blend_func=("src_alpha", "one_minus_src_alpha"),
    )
    for i in range(1, len(points)):
        draw_sphere(points[i], Color("green" if detected[i] else "grey"))
        label(points[i], pathid + str(i))
        new_pathid = pathid.upper() if detected[i] else pathid
        print(new_pathid + str(i) + ": " + str(energies[i - 1]))
    view.add(path)


def draw_annihilation(locs1, energies1, detected1, locs2, energies2, detected2, origin):
    draw_path([origin] + locs1, [0] + detected1, "a", energies1)
    print("\n")
    draw_path([origin] + locs2, [0] + detected2, "b", energies2)
    # marker = scene.visuals.Markers(pos=np.array([origin]), size=10, face_color="red")
    draw_sphere(origin, Color("red"))


for inner_radius in detector_inner_radii:
    draw_tube(inner_radius, inner_radius + detector_thickness, detector_length, 30)
draw_cylinder(10.6, 4, 30)
origin, locs1, locs2, energies1, energies2, detected1, detected2 = read_file(
    "../data/visualization.data"
)
draw_annihilation(locs1, energies1, detected1, locs2, energies2, detected2, origin)
text = ""

canvas.show()
# Run the application
if __name__ == "__main__":
    canvas.app.run()
# Show the canvas
