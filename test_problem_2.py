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


def solve_test_problem(mesh_fname, eps=0.1):
    domain, cell_markers, facet_markers = dolfinx.io.gmshio.read_from_msh(mesh_fname, MPI.COMM_WORLD, gdim=2)

    V = dolfinx.fem.functionspace(domain, ("Lagrange", 1))
    u = ufl.TrialFunction(V)
    v = ufl.TestFunction(V)

    Dirichlet_condition = 0
    f = 1
    #eps = 1e-3 ** 0.5

    f = dolfinx.fem.Constant(domain, PETSc.ScalarType(f))

    a = eps * ufl.dot(ufl.grad(u), ufl.grad(v)) * ufl.dx + u*v * ufl.dx
    L = f * v * ufl.dx

    tdim = domain.topology.dim
    fdim = tdim - 1
    domain.topology.create_connectivity(fdim, tdim)
    boundary_facets = dolfinx.mesh.exterior_facet_indices(domain.topology)
    boundary_dofs = dolfinx.fem.locate_dofs_topological(V, fdim, boundary_facets)
    
    u_d = dolfinx.fem.Function(V)
    u_d.interpolate(lambda x: np.full((x.shape[1],), Dirichlet_condition))

    bc = dolfinx.fem.dirichletbc(u_d, boundary_dofs)

    problem = dolfinx.fem.petsc.LinearProblem(a, L, bcs=[bc], petsc_options={"ksp_type": "preonly", "pc_type": "lu"})
    u = problem.solve()

    #show_function(u, V)

    # ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u)))   ufl.exp(ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u))) ** 0.8)
    weight_function_ufl = 0.1 + ufl.sqrt(ufl.dot(ufl.grad(u), ufl.grad(u))) ** 1
    weight_function_expresion = dolfinx.fem.Expression(weight_function_ufl, V.element.interpolation_points())
    weight_function = dolfinx.fem.Function(V)
    weight_function.interpolate(weight_function_expresion)

    points = domain.geometry.x[:, :2]
    triangles = domain.topology.connectivity(2, 0).array.reshape(-1, 3)
    values_weight = weight_function.x.array
    values_solution = u.x.array

    #show_function(weight_function, V)
    
    return points, triangles, values_weight, values_solution


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


def plot_function2(fname, mesh_fname, triangulation, values):
    gmsh.initialize()
    gmsh.open(mesh_fname)
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.getBoundingBox(-1, -1)
    gmsh.finalize()

    plt.figure(figsize=np.array([xmax - xmin, (ymax - ymin) * 1.15]))
    
    plt.tripcolor(triangulation, values, shading='gouraud')

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

    for j, eps in enumerate((0.1, 0.01, 0.001,), 1):
        fname_mesh = f'meshes/iterative_adaptation_2_{j}'
        fname_images = f'images/iterative_adaptation_2_{j}'

        cur_mesh = f'{fname_mesh}_initial.msh'
        adaptation_v2.adaptation(adaptation_v2.example, [lambda x: 1, lambda x: 1], n, niter, cur_mesh, run_gmsh=False)

        plotting.plot_triangle_mesh(f'{fname_images}_initial_mesh.pdf', cur_mesh, linewidth=1, markersize=0)

        for i in range(1, 5 + 1):
            previous_mesh = cur_mesh
            cur_mesh = f'{fname_mesh}_step_{i}.msh'

            points, triangles, values_weight, values_solution = solve_test_problem(previous_mesh, eps)

            triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
            interpolator = tri.LinearTriInterpolator(triang, values_weight)
            
            plot_function(f'{fname_images}_step_{i}_function.pdf', previous_mesh, triang, values_weight)
            plot_function2(f'{fname_images}_step_{i}_solution.pdf', previous_mesh, triang, values_solution)

            wight_function = lambda x: interpolator(*x[:2])
            adaptation_v2.adaptation(adaptation_v2.example, [wight_function, wight_function], n, niter, cur_mesh, run_gmsh=False)
            plotting.plot_triangle_mesh(f'{fname_images}_step_{i}_mesh.pdf', cur_mesh, linewidth=1, markersize=0)


def iterative_adaptation2():
    n = 30
    niter = 100

    for j, eps in enumerate((0.001,), 4):
        fname_mesh = f'meshes/iterative_adaptation_2_{j}'
        fname_images = f'images/iterative_adaptation_2_{j}'

        cur_mesh = f'{fname_mesh}_initial.msh'
        adaptation_v2.adaptation(adaptation_v2.example, [lambda x: 1, lambda x: 1], n, niter, cur_mesh, run_gmsh=False)

        plotting.plot_triangle_mesh(f'{fname_images}_initial_mesh.pdf', cur_mesh, linewidth=1, markersize=0)

        for i in range(1, 5 + 1):
            previous_mesh = cur_mesh
            cur_mesh = f'{fname_mesh}_step_{i}.msh'

            points, triangles, values_weight, values_solution = solve_test_problem(previous_mesh, eps)

            triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
            interpolator = tri.LinearTriInterpolator(triang, values_weight)
            
            plot_function(f'{fname_images}_step_{i}_function.pdf', previous_mesh, triang, values_weight)
            plot_function2(f'{fname_images}_step_{i}_solution.pdf', previous_mesh, triang, values_solution)

            wight_function = lambda x: interpolator(*x[:2])
            adaptation_v2.adaptation(adaptation_v2.example, [wight_function, wight_function], n, niter, cur_mesh, run_gmsh=False)
            plotting.plot_triangle_mesh(f'{fname_images}_step_{i}_mesh.pdf', cur_mesh, linewidth=1, markersize=0)


if __name__ == '__main__':
    
    iterative_adaptation2()
    exit()

    n = 40
    niter = 100

    mesh_fname = 'meshes/test_problem.msh'
    #mesh_fname = 'meshes/iterative_adaptation_2_step_5.msh'

    #adaptation_v2.adaptation(adaptation_v2.example, [lambda x: 1, lambda x: 1], n, niter, mesh_fname, run_gmsh=False)

    points, triangles, values_weight, values_solution = solve_test_problem(mesh_fname)

    triang = tri.Triangulation(points[:, 0], points[:, 1], triangles)
    #interpolator = tri.CubicTriInterpolator(triang, values)
    interpolator = tri.LinearTriInterpolator(triang, values_weight)

    plotting.plot_triangle_mesh('images/test_mesh.pdf', mesh_fname, linewidth=1, markersize=0)
    plot_function('images/test_function.pdf', mesh_fname, triang, values_weight)

    exit()

    wight_function = lambda x: interpolator(*x[:2])
    result_mesh_fname = 'meshes/test_problem_res.msh'
    adaptation_v2.adaptation(adaptation_v2.example, [wight_function, wight_function], n, niter, result_mesh_fname, run_gmsh=False)
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
