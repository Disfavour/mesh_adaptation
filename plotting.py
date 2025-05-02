import gmsh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection
from matplotlib.patches import PathPatch, Circle
from matplotlib.path import Path
import matplotlib as mpl


def create_base_entities_patches():
    A = (0.0000, 3.6000)
    B = (-3.1177, -1.8000)
    C = (3.1177, -1.8000)
    triangle = np.array([A, B, C, A])

    # Круг (многоугольник по окружности)
    circle_center = (0, 0)
    radius = 1
    theta = np.linspace(0, 2*np.pi, 100)
    circle = np.stack((circle_center[0] + radius * np.cos(theta),
                    circle_center[1] + radius * np.sin(theta)), axis=1)

    # Создаем Path: одна внешняя граница и одна дырка
    vertices = np.concatenate([triangle, circle[::-1]])  # обратный обход круга
    codes = ([Path.MOVETO] + [Path.LINETO]*2 + [Path.CLOSEPOLY] +  # треугольник
            [Path.MOVETO] + [Path.LINETO]*(len(circle)-2) + [Path.CLOSEPOLY])  # круг

    path = Path(vertices, codes)
    patch = PathPatch(path, facecolor=(1, 0, 0, 0.3), edgecolor='none')

    circle_patch = Circle(circle_center, radius=radius, facecolor=(0, 1, 0, 0.3), edgecolor='none')

    return patch, circle_patch


def plot_base_entities(mesh_name, fname):
    gmsh.initialize()
    gmsh.open(mesh_name)
    gmsh.model.mesh.renumber_nodes()
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    A = (0.0000, 3.6000)
    B = (-3.1177, -1.8000)
    C = (3.1177, -1.8000)
    triangle = np.array([A, B, C, A])

    # Круг (многоугольник по окружности)
    circle_center = (0, 0)
    radius = 1
    theta = np.linspace(0, 2*np.pi, 100)
    circle = np.stack((circle_center[0] + radius * np.cos(theta),
                    circle_center[1] + radius * np.sin(theta)), axis=1)

    # Создаем Path: одна внешняя граница и одна дырка
    vertices = np.concatenate([triangle, circle[::-1]])  # обратный обход круга
    codes = ([Path.MOVETO] + [Path.LINETO]*2 + [Path.CLOSEPOLY] +  # треугольник
            [Path.MOVETO] + [Path.LINETO]*(len(circle)-2) + [Path.CLOSEPOLY])  # круг

    path = Path(vertices, codes)
    patch = PathPatch(path, facecolor=(1, 0, 0, 1), edgecolor='none')

    circle_patch = Circle(circle_center, radius=radius, facecolor=(0, 1, 0, 1), edgecolor='none')

    fig, ax = plt.subplots(figsize=np.array([xmax - xmin, ymax - ymin]))

    ax.add_patch(patch)
    ax.add_patch(circle_patch)

    #plt.margins(x=0, y=0)
    ax.axis('scaled')
    ax.axis(False)
    fig.tight_layout(pad=0)

    fig.savefig(fname, transparent=True)
    plt.close()


def plot_algorithm_steps(mesh_name, fname):
    gmsh.initialize()
    gmsh.open(mesh_name)
    gmsh.model.mesh.renumber_nodes()
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)
    patch, circle_patch = create_base_entities_patches()

    fig, ax = plt.subplots(figsize=np.array([xmax - xmin, ymax - ymin]))

    ax.add_patch(patch)
    ax.add_patch(circle_patch)

    ax.triplot(triangulation, 'o-b')

    #plt.margins(x=0, y=0)
    ax.axis('scaled')
    ax.axis(False)
    fig.tight_layout(pad=0)

    fig.savefig(fname, transparent=True)
    plt.close()


def plot_algorithm_node_affiliation(mesh_name, fname):
    gmsh.initialize()
    gmsh.open(mesh_name)
    gmsh.model.mesh.renumber_nodes()
    plane_surfaces = gmsh.model.get_entities(2)
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)

    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]

    triangle_nodes_for_surfaces = []
    node_coords_for_surfaces = []
    for plane_surface in plane_surfaces:
        triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type, plane_surface[1])
        uniq = np.unique(triangle_nodes)
        mapping = {old: new for new, old in zip(range(uniq.size), uniq)}
        triangle_nodes_for_surfaces.append(np.array([mapping[node] for node in triangle_nodes]).reshape(-1, 3))
        node_coords_for_surfaces.append(node_coords[np.isin(node_tags, uniq, assume_unique=True)])

    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulations = []
    for triangle_nodes_surf, node_coords in zip(triangle_nodes_for_surfaces, node_coords_for_surfaces):
        triangulations.append(Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes_surf))
    patch, circle_patch = create_base_entities_patches()

    fig, ax = plt.subplots(figsize=np.array([xmax - xmin, ymax - ymin]))

    ax.add_patch(patch)
    ax.add_patch(circle_patch)

    for tri, c in zip(triangulations, ('g', 'r')):
        ax.triplot(tri, f'o-{c}')

    #plt.margins(x=0, y=0)
    ax.axis('scaled')
    ax.axis(False)
    fig.tight_layout(pad=0)

    fig.savefig(fname, transparent=True)
    plt.close()


def plot_refining(mesh_name, fname):
    gmsh.initialize()
    gmsh.open(mesh_name)
    gmsh.model.mesh.renumber_nodes()
    plane_surfaces = gmsh.model.get_entities(2)
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)

    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]

    triangle_nodes_for_surfaces = []
    node_coords_for_surfaces = []
    for plane_surface in plane_surfaces:
        triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type, plane_surface[1])
        uniq = np.unique(triangle_nodes)
        mapping = {old: new for new, old in zip(range(uniq.size), uniq)}
        triangle_nodes_for_surfaces.append(np.array([mapping[node] for node in triangle_nodes]).reshape(-1, 3))
        node_coords_for_surfaces.append(node_coords[np.isin(node_tags, uniq, assume_unique=True)])

    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulations = []
    for triangle_nodes_surf, node_coords in zip(triangle_nodes_for_surfaces, node_coords_for_surfaces):
        triangulations.append(Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes_surf))
    patch, circle_patch = create_base_entities_patches()

    fig, ax = plt.subplots(figsize=np.array([xmax - xmin, ymax - ymin * 1.3]))

    #ax.add_patch(patch)
    #ax.add_patch(circle_patch)
    weight_functions = [lambda x: (1/2 * (-np.linalg.norm(x) + 1) - 1) ** 4, lambda x: (1/2 * (np.linalg.norm(x) - 1) - 1) ** 4]

    weight_for_nodes = []
    for node_coords, function in zip(node_coords_for_surfaces, weight_functions):
        weight_for_nodes.append([function(x) for x in node_coords])

    #weight_for_nodes = np.array(weight_for_nodes)
    vmin, vmax = min(weight_for_nodes[0] + weight_for_nodes[1]), max(weight_for_nodes[0] + weight_for_nodes[1])
    #print(vmin, vmax)

    some_plot = None
    for tri, weights in zip(triangulations, weight_for_nodes):
        some_plot = ax.tripcolor(tri, weights, vmin=vmin, vmax=vmax, cmap='viridis', shading='gouraud')

    for tri, c in zip(triangulations, ('g', 'r')):
        ax.triplot(tri, f'o-{c}')
    
    fig.colorbar(some_plot, location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    ax.axis('scaled')
    ax.axis(False)
    fig.tight_layout(pad=0)

    fig.savefig(fname, transparent=True)
    plt.close()


def plot_triangle_mesh(fname, mesh_fname, linewidth=1.5, markersize=6):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    gmsh.model.mesh.renumber_nodes()
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, ymax - ymin]))
    plt.triplot(triangulation, 'o-b', linewidth=linewidth, markersize=markersize)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()


def plot_function(fname, mesh_fname, f):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    gmsh.model.mesh.renumber_nodes()
    for i in range(3):
        gmsh.model.mesh.refine()
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.15]))
    
    some_plot = plt.tripcolor(triangulation, [f(x) for x in node_coords], shading='gouraud', norm=mpl.colors.LogNorm())

    plt.colorbar(some_plot, location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()

def plot_function2(fname, mesh_fname, f):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    gmsh.model.mesh.renumber_nodes()
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.15]))
    
    weights = [f(x) for x in node_coords]
    print(min(weights), max(weights))
    some_plot = plt.tripcolor(triangulation, weights, shading='gouraud')#, norm=mpl.colors.LogNorm())

    plt.colorbar(some_plot, location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()


def plot_function_contour(fname, mesh_fname, f):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    gmsh.model.mesh.renumber_nodes()
    #for i in range(3):
    #    gmsh.model.mesh.refine()
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.15]))
    
    #some_plot = plt.tripcolor(triangulation, [f(x) for x in node_coords], shading='gouraud', norm=mpl.colors.LogNorm())

    plt.tricontourf(triangulation, [f(x) for x in node_coords], levels=50)

    #plt.colorbar(some_plot, location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()


def plot_triangle_mesh_quality(quality, fname, mesh_fname):
    gmsh.initialize()
    
    gmsh.open(mesh_fname)

    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)

    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]

    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)

    triangle_qualities = gmsh.model.mesh.get_element_qualities(triangle_tags, qualityName=quality)

    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.2]) * 1.5)
    plt.tripcolor(triangulation, triangle_qualities)
    plt.triplot(triangulation, 'o-k')
    plt.colorbar(location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()


def test():
    x = np.linspace(-3, 3, 100)
    y = np.linspace(-3, 3, 100)
    X, Y = np.meshgrid(x, y)
    Z = np.sin(np.sqrt(X**2 + Y**2))  # пример 2D-функции

    plt.contourf(X, Y, Z, levels=50, cmap='viridis')
    plt.colorbar(label='z = f(x, y)')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.axis('equal')
    plt.show()





if __name__ == '__main__':

    # article_dir = os.path.join('images', 'article_1_black')
    # os.makedirs(article_dir, exist_ok=True)
    #plot_triangle_mesh('images/example_1.pdf', 'meshes/example_1.msh')

    #plot_base_entities('meshes/example_2_step_1.msh', 'images/base_entities.pdf')

    #plot_algorithm_steps('meshes/example_2_step_1.msh', 'images/algorithm_step_1.pdf')
    #plot_algorithm_steps('meshes/example_2_step_2.msh', 'images/algorithm_step_2.pdf')
    #plot_algorithm_steps('meshes/example_2_step_3.msh', 'images/algorithm_step_3.pdf')
    #plot_algorithm_steps('meshes/example_2.msh', 'images/algorithm_result.pdf')

    #plot_algorithm_node_affiliation('meshes/example_2_step_3.msh', 'images/algorithm_step_3.pdf')
    #plot_algorithm_node_affiliation('meshes/example_2.msh', 'images/algorithm_result.pdf')

    #plot_refining('meshes/example_2_refine.msh', 'images/algorithm_refine.pdf')

    qualities = [
    'minDetJac',    # the adaptively computed minimal Jacobian determinant
    'maxDetJac',    # the adaptively computed maximal Jacobian determinant
    'minSJ',        # sampled minimal scaled jacobien
    'minSICN',      # sampled minimal signed inverted condition number
    'minSIGE',      # sampled signed inverted gradient error
    'gamma',        # ratio of the inscribed to circumcribed sphere radius
    'innerRadius',
    'outerRadius',
    'minIsotropy',  # minimum isotropy measure
    'angleShape',   # angle shape measure
    'minEdge',      # minimum straight edge length
    'maxEdge',      # maximum straight edge length
    'volume',
    ]
    
    #plot_triangle_mesh_quality('gamma', 'images/test_gamma.pdf', 'meshes/example_1.msh')
    #plot_triangle_mesh_quality('volume', 'images/test_volume.pdf', 'meshes/example_1.msh')

    i = 11
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(x) - 3) ** 2 / 0.5 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 12
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(x) - 1.5) ** 2 / 0.25 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 13
    # f = lambda x: max(
    #     1.1 - 1 * np.exp(-(np.abs(x[0]) - 0) ** 16 / 3000 ** 2),
    #     1.1 - 1 * np.exp(-(np.abs(x[1]) - 0) ** 16 / 100 ** 2))
    f = lambda x: (1.1 - 1 * np.exp(-(np.abs(x[0]) - 0) ** 16 / 3000 ** 2)) * (1.1 - 1 * np.exp(-(np.abs(x[1]) - 0) ** 16 / 100 ** 2))
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 14
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(x - np.array((3, -0.5))) - 0) ** 2 / 0.5 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 15
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(x - np.array((-1, 0.5))) - 0) ** 2 / 0.5 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 16
    (x1, y1), (x2, y2) = (-1, -2), (0.7, 2)
    distance = lambda x: abs((y2 - y1)*x[0] - (x2 - x1)*x[1] + x2*y1 - y2*x1) / np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(distance(x)) - 0) ** 2 / 0.5 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)

    i = 17
    (x1, y1), (x2, y2) = (-1, -2), (0.7, 2)
    distance = lambda x: abs((y2 - y1)*x[0] - (x2 - x1)*x[1] + x2*y1 - y2*x1) / np.sqrt((y2 - y1)**2 + (x2 - x1)**2)
    f = lambda x: 0.1 + 1 * np.exp(-(np.linalg.norm(distance(x)) - 0) ** 2 / 0.25 ** 2)
    plot_triangle_mesh(f'images/example_{i}_step_3.pdf', f'meshes/example_{i}_step_3.msh', linewidth=1, markersize=0)
    plot_triangle_mesh(f'images/example_{i}.pdf', f'meshes/example_{i}.msh', linewidth=1, markersize=0)
    plot_function(f'images/example_{i}_function.pdf', f'meshes/example_{i}.msh', f)