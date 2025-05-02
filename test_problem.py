from mpi4py import MPI
import dolfinx
import ufl
import dolfinx.fem.petsc
import basix.ufl
import numpy as np
import pyvista
from petsc4py import PETSc
import gmsh
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import matplotlib as mpl
import adaptation_v2
import plotting


def solve_test_problem(mesh_fname):
    domain, cell_markers, facet_markers = dolfinx.io.gmshio.read_from_msh(mesh_fname, MPI.COMM_WORLD, gdim=2)

    Q = dolfinx.fem.functionspace(domain, ("DG", 0))
    k = dolfinx.fem.Function(Q)

    k1, k2 = 1, 10

    inner_circle_tag = 1
    inner_circle_cells = cell_markers.find(inner_circle_tag)
    k.x.array[inner_circle_cells] = np.full_like(inner_circle_cells, k1, dtype=dolfinx.default_scalar_type)

    outer_circle_tag = 2
    outer_circle_cells = cell_markers.find(outer_circle_tag)
    k.x.array[outer_circle_cells] = np.full_like(outer_circle_cells, k2, dtype=dolfinx.default_scalar_type)

    #show_markers(domain, cell_markers)

    V = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    # f = - ufl.div(k * ufl.grad(u_ufl))
    # a = ufl.dot(k * ufl.grad(u), ufl.grad(v)) * ufl.dx
    # L = f * v * ufl.dx + ufl.dot(k*ufl.grad(u_ufl), n) * v * ufl.ds

    Neumann_condition = 0
    Dirichlet_condition = 0
    f = 0

    Neumann_condition = dolfinx.fem.Constant(domain, PETSc.ScalarType(Neumann_condition))
    #Dirichlet_condition = dolfinx.fem.Constant(domain, PETSc.ScalarType(Dirichlet_condition))
    f = dolfinx.fem.Constant(domain, PETSc.ScalarType(1))

    a = ufl.dot(k * ufl.grad(u), ufl.grad(v)) * ufl.dx
    L = f * v * ufl.dx  + Neumann_condition * v * ufl.ds

    def boundary_Dirichlet(x):
        x = x.T
        #return np.isclose(np.linalg.norm(x, axis=1), 3)
        return np.logical_and(np.isclose(np.linalg.norm(x, axis=1), 3), x[:, 1] > x[:, 0] - 3 + 1e-3)
    
    u_d = dolfinx.fem.Function(V)
    u_d.interpolate(lambda x: np.full((x.shape[1],), Dirichlet_condition))
    boundary_dofs = dolfinx.fem.locate_dofs_geometrical(V, boundary_Dirichlet)
    bc = dolfinx.fem.dirichletbc(u_d, boundary_dofs)

    problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    u = problem.solve()

    #show_function(u, V)

    # ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u)))   ufl.exp(ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u))) ** 0.8)
    weight_function_ufl = 0 + ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u))) ** 2
    weight_function_expresion = dolfinx.fem.Expression(weight_function_ufl, V.element.interpolation_points())
    weight_function = dolfinx.fem.Function(V)
    weight_function.interpolate(weight_function_expresion)

    points = domain.geometry.x[:, :2]
    triangles = domain.topology.connectivity(2, 0).array.reshape(-1, 3)
    values = weight_function.x.array

    #show_function(weight_function, V)
    
    return points, triangles, values


def show_function(u, V):
    p2 = pyvista.Plotter(window_size=[800, 800])
    grid_uh = pyvista.UnstructuredGrid(*dolfinx.plot.vtk_mesh(V))
    grid_uh.point_data["u"] = u.x.array.real
    grid_uh.set_active_scalars("u")
    p2.add_mesh(grid_uh, show_edges=True)
    p2.show()
    #p2.close()
    

def show_markers(domain, cell_markers):
    tdim = domain.topology.dim
    domain.topology.create_connectivity(tdim, tdim)
    topology, cell_types, x = dolfinx.plot.vtk_mesh(domain, tdim)
    grid = pyvista.UnstructuredGrid(topology, cell_types, x)
    num_local_cells = domain.topology.index_map(tdim).size_local
    grid.cell_data["Marker"] = cell_markers.values[cell_markers.indices < num_local_cells]
    grid.set_active_scalars("Marker")

    p = pyvista.Plotter(window_size=[800, 800])
    p.add_mesh(grid, show_edges=True)
    p.show()
    #p.close()


def plot_function(fname, mesh_fname, triangulation, values):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.15]))
    
    plt.tripcolor(triangulation, values, shading='gouraud', norm=mpl.colors.LogNorm())

    plt.colorbar(location='bottom', shrink=0.9, fraction=0.05, pad=0)

    #plt.margins(x=0, y=0)
    plt.axis('scaled')
    plt.axis(False)
    plt.tight_layout(pad=0)

    plt.savefig(fname, transparent=True)
    #plt.show()
    plt.close()


def iterative_adaptation():
    n = 40
    niter = 100

    mesh_fname = 'meshes/test_problem.msh'

    fname = 'meshes/iterative_adaptation'

    cur_mesh = f'{fname}_initial.msh'
    adaptation_v2.adaptation(adaptation_v2.example, [lambda x: 1, lambda x: 1], n, niter, f'{fname}_initial.msh', run_gmsh=False)

    plotting.plot_triangle_mesh(f'images/iterative_adaptation_initial_mesh.pdf', cur_mesh, linewidth=1, markersize=0)

    for i in range(1, 5 + 1):
        previous_mesh = cur_mesh
        cur_mesh = f'{fname}_step_{i}.msh'

        points, triangles, values = solve_test_problem(previous_mesh)

        triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
        #interpolator = tri.CubicTriInterpolator(triang, values)
        interpolator = tri.LinearTriInterpolator(triang, values)
        
        plot_function(f'images/iterative_adaptation_step_{i}_function.pdf', previous_mesh, triang, values)

        wight_function = lambda x: interpolator(*x[:2])
        adaptation_v2.adaptation(adaptation_v2.example, [wight_function, wight_function], n, niter, cur_mesh, run_gmsh=False)
        plotting.plot_triangle_mesh(f'images/iterative_adaptation_step_{i}_mesh.pdf', cur_mesh, linewidth=1, markersize=0)


if __name__ == '__main__':
    
    iterative_adaptation()
    exit()

    n = 40
    niter = 100

    mesh_fname = 'meshes/test_problem.msh'

    adaptation_v2.adaptation(adaptation_v2.example, [lambda x: 1, lambda x: 1], n, niter, mesh_fname, run_gmsh=False)

    points, triangles, values = solve_test_problem(mesh_fname)

    triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
    #interpolator = tri.CubicTriInterpolator(triang, values)
    interpolator = tri.LinearTriInterpolator(triang, values)

    plotting.plot_triangle_mesh('images/test_mesh.pdf', mesh_fname, linewidth=1, markersize=0)
    plot_function('images/test_function.pdf', mesh_fname, triang, values)

    #exit()

    wight_function = lambda x: interpolator(*x[:2])
    result_mesh_fname = 'meshes/test_problem_res.msh'
    adaptation_v2.adaptation(adaptation_v2.example, [wight_function, wight_function], n, niter, result_mesh_fname, run_gmsh=True)
    plotting.plot_triangle_mesh('images/test_mesh_res.pdf', result_mesh_fname, linewidth=1, markersize=0)

    

    #from plotting import plot_function2, plot_function_contour, plot_function2, plot_triangle_mesh
    #plot_triangle_mesh('images/test3.pdf', 'meshes/test_problem1.msh', linewidth=1, markersize=0)
    #plot_function2('images/test1.pdf', 'meshes/test_problem.msh', weight_function)

    # plot_function_contour('images/test.pdf', 'meshes/test_problem.msh', weight_function)

    # import matplotlib.pyplot as plt
    # import matplotlib.tri as tri
    # triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
    # asd = plt.tripcolor(triang, values,shading='gouraud', cmap='viridis')
    # plt.colorbar(asd, label="u")
    # plt.xlabel("x")
    # plt.ylabel("y")
    # plt.title("Решение задачи")
    # plt.show()
