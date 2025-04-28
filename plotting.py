import gmsh
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.tri import Triangulation
from matplotlib.collections import PolyCollection


def plot_triangle_mesh(fname, mesh_fname):
    gmsh.initialize()
    
    gmsh.open(mesh_fname)

    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)

    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)[:, :2]

    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)

    gmsh.finalize()

    triangulation = Triangulation(node_coords[:, 0], node_coords[:, 1], triangle_nodes - 1)

    plt.figure(figsize=np.array([xmax - xmin, ymax - ymin]) * 1.5)
    plt.triplot(triangulation, 'o-k')

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


if __name__ == '__main__':

    # article_dir = os.path.join('images', 'article_1_black')
    # os.makedirs(article_dir, exist_ok=True)
    plot_triangle_mesh('images/example_1.pdf', 'meshes/example_1.msh')

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
    
    plot_triangle_mesh_quality('gamma', 'images/test_gamma.pdf', 'meshes/example_1.msh')
    plot_triangle_mesh_quality('volume', 'images/test_volume.pdf', 'meshes/example_1.msh')