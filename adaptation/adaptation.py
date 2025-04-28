import gmsh
import numpy as np


def adaptation(define_entities_and_weight_functions, fname=None, fnames_for_steps=[None, None, None], run_gmsh=False):
    gmsh.initialize()
        
    curves, curves_of_surfaces, fixed_points, weight_functions, initial_triangles_in_row, number_of_iterations = define_entities_and_weight_functions()

    curve_loops = [gmsh.model.occ.add_curve_loop([curve]) for curve in curves]

    surface_to_curve = {}
    for surface_curves in curves_of_surfaces:
        surface_tag = gmsh.model.occ.add_plane_surface([curve_loops[curves.index(curve)] for curve in surface_curves])
        surface_to_curve[surface_tag] = surface_curves

    gmsh.model.occ.synchronize()

    physical_groups = [gmsh.model.add_physical_group(2, [plane_surface]) for plane_surface in surface_to_curve]

    plane_surfaces = gmsh.model.get_entities(2)
    base_surface = plane_surfaces[0]

    # algorithm step 1
    create_initial_uniform_mesh(base_surface, initial_triangles_in_row)

    if fnames_for_steps[0] is not None:
        gmsh.write(fnames_for_steps[0])

    # algorithm step 2
    fixed_nodes = fix_specified_boundary_nodes(fixed_points)
    relocate_nearest_nodes_to_boundaries(surface_to_curve)

    # node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    # node_to_coords = {node: coords for node, coords in zip(node_tags, node_coords.reshape(-1, 3))}
    # for n in [102, 203, 304, 405]:
    #     n = np.uint64(n)
    #     coord = node_to_coords[n]
    #     inside = gmsh.model.is_inside(1, 5, coord)
    #     print(n, coord, inside)
    # gmsh.fltk.run()
    # exit()

    if fnames_for_steps[1] is not None:
        gmsh.write(fnames_for_steps[1])
    
    # algorithm step 3
    curve_to_nodes, plane_surface_to_nodes = classify_nodes()
    map_nodes_and_elements_to_entities(base_surface, curve_to_nodes, plane_surface_to_nodes)

    if fnames_for_steps[2] is not None:
        gmsh.write(fnames_for_steps[2])
    
    #gmsh.fltk.run()
    #exit()
    
    boundary_curves = gmsh.model.get_entities(1)
    nodes_of_plane_surfaces = [gmsh.model.mesh.get_nodes(*plane_surface, returnParametricCoord=False)[0] for plane_surface in plane_surfaces]
    nodes_of_boundary_curves = [np.setdiff1d(gmsh.model.mesh.get_nodes(*curve, returnParametricCoord=False)[0], fixed_nodes, assume_unique=True) for curve in boundary_curves]

    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_to_coords = {node: coords for node, coords in zip(node_tags, node_coords.reshape(-1, 3))}
    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_to_nodes = {triangle: nodes for triangle, nodes in zip(triangle_tags, triangle_nodes.reshape(-1, 3))}
    node_to_triangles = reverse_dict(triangle_to_nodes)

    # algorithm step 4 (loop iteration)
    for i in range(number_of_iterations):
        for plane_surface, inner_nodes, weight_function in zip(plane_surfaces, nodes_of_plane_surfaces, weight_functions):
            for inner_node in inner_nodes:
                adjacent_triangles = node_to_triangles[inner_node]
                triangle_areas = gmsh.model.mesh.get_element_qualities(adjacent_triangles, qualityName='volume')
                coords_of_triangle_nodes = [[node_to_coords[node] for node in triangle_to_nodes[adjacent_triangle]] for adjacent_triangle in adjacent_triangles]
                barycenters = np.sum(coords_of_triangle_nodes, axis=1) / 3
                weighted_areas = np.array([weight_function(barycenter) * triangle_area for barycenter, triangle_area in zip(barycenters, triangle_areas)])
                barycenter_area_product = barycenters * weighted_areas.reshape(-1, 1)
                new_coord = barycenter_area_product.sum(axis=0) / weighted_areas.sum()
                gmsh.model.mesh.set_node(inner_node, new_coord, [])
                node_to_coords[inner_node] = new_coord

        for boundary_curve, inner_nodes in zip(boundary_curves, nodes_of_boundary_curves):
            for inner_node in inner_nodes:
                adjacent_triangles = node_to_triangles[inner_node]
                triangle_areas = gmsh.model.mesh.get_element_qualities(adjacent_triangles, qualityName='volume')
                coords_of_triangle_nodes = [[node_to_coords[node] for node in triangle_to_nodes[adjacent_triangle]] for adjacent_triangle in adjacent_triangles]
                barycenters = np.sum(coords_of_triangle_nodes, axis=1) / 3
                barycenter_area_product = barycenters * triangle_areas.reshape(-1, 1)
                new_coord = barycenter_area_product.sum(axis=0) / triangle_areas.sum()
                new_coord = gmsh.model.get_closest_point(*boundary_curve, new_coord)[0]
                gmsh.model.mesh.set_node(inner_node, new_coord, [])
                node_to_coords[inner_node] = new_coord

    gmsh.model.mesh.renumber_nodes()

    #gmsh.option.setNumber("Mesh.MshFileVersion", 2)

    if fname is not None:
        gmsh.write(fname)

    if run_gmsh:
        gmsh.fltk.run()

    gmsh.finalize()


def reverse_dict(d):
    res = {}
    for k, vs in d.items():
        for v in vs:
            if v not in res:
                res[v] = []
            res[v].append(k)
    return res


def create_initial_uniform_mesh(plane_surface, triangles_in_row):
    xmin, ymin, zmin, xmax, ymax, zmax = gmsh.model.get_bounding_box(-1, -1)
    triangle_edge_lenght = (xmax - xmin) / triangles_in_row
    triangle_height = np.sqrt(triangle_edge_lenght ** 2 - (triangle_edge_lenght / 2) ** 2)
    xmin -= triangle_edge_lenght
    xmax += triangle_edge_lenght
    ymin -= triangle_height
    ymax += triangle_height
    #xmin, ymin, zmin, xmax, ymax, zmax = map(lambda x: 1.5 * x, gmsh.model.get_bounding_box(-1, -1))

    triangle_edge_lenght = (xmax - xmin) / triangles_in_row
    triangle_height = np.sqrt(triangle_edge_lenght ** 2 - (triangle_edge_lenght / 2) ** 2)
    triangles_in_column = int(np.ceil((ymax - ymin) / triangle_height))
    ymax = ymin + triangles_in_column * triangle_height

    x, y = np.meshgrid(np.linspace(xmin, xmax, triangles_in_row + 1), np.linspace(ymin, ymax, triangles_in_column + 1), indexing='xy')
    x[1::2] += triangle_edge_lenght / 2
    node_tags = range(1, (triangles_in_column + 1) * (triangles_in_row + 1) + 1)
    gmsh.model.mesh.add_nodes(*plane_surface, node_tags, [coord for xyz in zip(x.flat, y.flat, np.zeros(x.size)) for coord in xyz])

    node_tags = np.array(node_tags).reshape(x.shape)
    triangles_nodes = []
    for row1, row2 in zip(node_tags[::2], node_tags[1::2]):
        for n1, n2, n3, n4 in zip(row1[:-1], np.roll(row1, -1), row2, np.roll(row2, -1)):
            triangles_nodes.extend((n1, n2, n3, n2, n4, n3))
    for row1, row2 in zip(node_tags[1::2], node_tags[2::2]):
        for n1, n2, n3, n4 in zip(row1[:-1], np.roll(row1, -1), row2, np.roll(row2, -1)):
            triangles_nodes.extend((n1, n4, n3, n1, n2, n4))

    gmsh.model.mesh.add_elements_by_type(plane_surface[1], gmsh.model.mesh.get_element_type("Triangle", 1), [], triangles_nodes)


def fix_specified_boundary_nodes(fixed_points):
    one = np.uint64(1)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)

    fixed_nodes = []
    for point_tag in fixed_points:
        coord = gmsh.model.getValue(0, point_tag, [])
        node_tag = node_tags[np.argmin(np.linalg.norm(node_coords - coord, axis=1))]
        gmsh.model.mesh.set_node(node_tag, coord, [])
        node_coords[node_tag - one] = coord
        fixed_nodes.append(node_tag)

    return fixed_nodes


def relocate_nearest_nodes_to_boundaries(surface_to_curve):
    one = np.uint64(1)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)

    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(gmsh.model.mesh.get_element_type("Triangle", 1))
    triangle_nodes = triangle_nodes.reshape(-1, 3)

    plane_surfaces = gmsh.model.get_entities(2)

    for nodes_of_triangle in triangle_nodes:

        node_plane_surfaces = []
        for node_of_triangle in nodes_of_triangle:
            node_plane_surfaces.append([])
            for plane_surface in plane_surfaces:
                if gmsh.model.is_inside(*plane_surface, node_coords[node_of_triangle - one]):
                    node_plane_surfaces[-1].append(plane_surface)

        # if all nodes in one surface
        if set(node_plane_surfaces[0]).intersection(*node_plane_surfaces[1:]) or not any(node_plane_surfaces):
            continue
        
        # i, j in one surface and k in another
        for i, j in zip(range(3), np.roll(range(3), -1)):
            if set(node_plane_surfaces[i]).intersection(node_plane_surfaces[j]) or not (node_plane_surfaces[i] or node_plane_surfaces[j]):
                break
        k = set(range(3)).difference([i, j]).pop()

        # only max 2 surfaces can be
        surfaces = tuple(set(node_plane_surfaces[i]).union(node_plane_surfaces[k]))

        boundary_curve = None
        if len(surfaces) == 2:
            boundary_curve = set(surface_to_curve[surfaces[0][1]]).intersection(surface_to_curve[surfaces[1][1]]).pop()
        elif len(surfaces) == 1:
            # one surface and outer space
            curves = surface_to_curve[surfaces[0][1]]
            node_coord = node_coords[nodes_of_triangle[k] - one]
            closest_coords_on_curves = np.array([gmsh.model.get_closest_point(1, curve, node_coord)[0] for curve in curves])
            closest_curve = curves[np.argmin(np.linalg.norm(closest_coords_on_curves - node_coord, axis=1))]
            boundary_curve = closest_curve
        else:
            raise Exception('more than 2 surfaces')
        
        coords_of_triangle_nodes = [node_coords[node_of_triangle - one] for node_of_triangle in nodes_of_triangle]
        closest_coords = [gmsh.model.get_closest_point(1, boundary_curve, coords)[0] for coords in coords_of_triangle_nodes]
        distances = np.linalg.norm(np.array(closest_coords) - np.array(coords_of_triangle_nodes), axis=1)
        move_distance_of_nodes_ij = distances[[i, j]].max()
        move_distance_of_node_k = distances[k]
        if move_distance_of_nodes_ij < move_distance_of_node_k:
            for index in [i, j]:
                gmsh.model.mesh.set_node(nodes_of_triangle[index], closest_coords[index], [])
                node_coords[nodes_of_triangle[index] - one] = closest_coords[index]
        else:
            gmsh.model.mesh.set_node(nodes_of_triangle[k], closest_coords[k], [])
            node_coords[nodes_of_triangle[k] - one] = closest_coords[k]


def classify_nodes():
    one = np.uint64(1)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)
    plane_surfaces = gmsh.model.get_entities(2)
    curves = gmsh.model.get_entities(1)
    curve_to_nodes = {curve: [] for curve in curves}
    plane_surface_to_nodes = {plane_surface: [] for plane_surface in plane_surfaces}
    not_in_curve = True
    for node in node_tags:
        for curve in curves:
            if gmsh.model.is_inside(*curve, node_coords[node - one]):
                curve_to_nodes[curve].append(node)
                not_in_curve = False
                break
        if not_in_curve:
            for plane_surface in plane_surfaces:
                if gmsh.model.is_inside(*plane_surface, node_coords[node - one]):
                    plane_surface_to_nodes[plane_surface].append(node)
                    break
        not_in_curve = True

    return curve_to_nodes, plane_surface_to_nodes


def map_nodes_and_elements_to_entities(base_surface, curve_to_nodes, plane_surface_to_nodes):
    one = np.uint64(1)
    node_tags, node_coords, _ = gmsh.model.mesh.get_nodes()
    node_coords = node_coords.reshape(-1, 3)

    triangle_type = gmsh.model.mesh.get_element_type("Triangle", 1)
    triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type)
    triangle_nodes = triangle_nodes.reshape(-1, 3)

    gmsh.model.mesh.remove_elements(*base_surface)
    gmsh.model.mesh.reclassify_nodes()  # delete nodes

    for curve, nodes_of_curve in curve_to_nodes.items():
        gmsh.model.mesh.add_nodes(*curve, nodes_of_curve, np.array([node_coords[node_tag - one] for node_tag in nodes_of_curve]).flat)
    
    for plane_surface, nodes_of_plane_surface in plane_surface_to_nodes.items():
        gmsh.model.mesh.add_nodes(*plane_surface, nodes_of_plane_surface, np.array([node_coords[node_tag - one] for node_tag in nodes_of_plane_surface]).flat)

    plane_surfaces = gmsh.model.get_entities(2)
    plane_surface_to_triangle = {plane_surface: [] for plane_surface in plane_surfaces}

    for plane_surface, nodes_of_plane_surface in plane_surface_to_nodes.items():
        nodes_of_triangles_in_surface = triangle_nodes[np.any(np.isin(triangle_nodes, nodes_of_plane_surface), axis=1)]
        gmsh.model.mesh.add_elements_by_type(plane_surface[1], triangle_type, [], nodes_of_triangles_in_surface.flat)
    
    # 3 nodes of a triangle on boundary (example square, triangle)
    for curve, nodes_of_curve in curve_to_nodes.items():
        nodes_of_triangles_in_curve = triangle_nodes[np.all(np.isin(triangle_nodes, nodes_of_curve), axis=1)]
        if nodes_of_triangles_in_curve.size > 0:
            gmsh.model.mesh.add_elements_by_type(plane_surface[1], triangle_type, [], nodes_of_triangles_in_curve.flat)

    # verification
    for curve, nodes_of_curve in curve_to_nodes.items():
        boundary_nodes = gmsh.model.mesh.get_nodes(*curve)[0]
        if set(nodes_of_curve) != set(boundary_nodes):
            raise Exception('curve nodes problem')
    
    for plane_surface, nodes_of_plane_surface in plane_surface_to_nodes.items():
        boundary_nodes = gmsh.model.mesh.get_nodes(*plane_surface)[0]
        if set(nodes_of_plane_surface) != set(boundary_nodes):
            raise Exception('plane_surface nodes problem')

    for plane_surface in plane_surfaces:
        inner_nodes = gmsh.model.mesh.get_nodes(*plane_surface)[0]
        for boundary_curve in gmsh.model.get_boundary([plane_surface], oriented=False):
            boundary_nodes = gmsh.model.mesh.get_nodes(*boundary_curve, True)[0]
            if set(inner_nodes).intersection(boundary_nodes):
                raise Exception('boundary nodes in inner')

    if gmsh.model.mesh.get_duplicate_nodes().size > 0:
        raise Exception('duplicate nodes')
    
    for plane_surface in plane_surfaces:
        triangle_tags, triangle_nodes = gmsh.model.mesh.get_elements_by_type(triangle_type, plane_surface[1])
        areas = gmsh.model.mesh.get_element_qualities(triangle_tags, qualityName='volume')
        if areas.size - np.count_nonzero(areas) > 0:
            raise Exception('triangles with zero area')


def example_1():
    z = 0

    circle = gmsh.model.occ.add_circle(0, 0, 0, 1)
    ellipse = gmsh.model.occ.add_ellipse(0, 0, 0, 2.5, 1.5)
    circle2 = gmsh.model.occ.add_circle(0, 0, 0, 3)

    flower_radius = 4
    N = 100
    points = []
    for i in range(N):
        fi = 2 * np.pi * i / N
        r = flower_radius + 0.5 * np.cos(8 * fi)
        x = -r * np.cos(fi)
        y = -r * np.sin(fi)
        points.append(gmsh.model.occ.addPoint(x, y, z))
    flower = gmsh.model.occ.add_bspline(points + [points[0]])

    # Любые фигуры, задаваемые по точкам, задаются через bspline и точки должны быть заданы против часовой стрелки
    # gmsh.model.occ.add_rectangle создает сразу plane_surface, а нужно curve
    rectangle_points = []
    for i in range(-5, 6):
        rectangle_points.append(gmsh.model.occ.add_point(5, i, z))
    for i in range(4, -6, -1):
        rectangle_points.append(gmsh.model.occ.add_point(i, 5, z))
    for i in range(4, -6, -1):
        rectangle_points.append(gmsh.model.occ.add_point(-5, i, z))
    for i in range(-4, 5):
        rectangle_points.append(gmsh.model.occ.add_point(i, -5, z))
    rectangle = gmsh.model.occ.add_bspline(rectangle_points + [rectangle_points[0]], degree=1)

    curves = [circle, ellipse, circle2, flower, rectangle]

    curves_of_surfaces = [
        [circle],
        [ellipse, circle],
        [circle2, ellipse],
        [flower, circle2],
        [rectangle, flower],
    ]

    # Угловые точки, которые хотим сохранить в сетке (например, точки в углах квадрата)
    fixed_points = [rectangle_points[0], rectangle_points[10], rectangle_points[20], rectangle_points[30]]

    # linear        1 - 1/distance * (np.linalg.norm(x) - inner_radius)
    # quadratic     (1/distance * (np.linalg.norm(x) - inner_radius) - 1) ** 2
    # np.exp(10*(-(np.linalg.norm(x) - inner_radius)))
    weight_functions = [lambda x: 1, lambda x: np.exp(5*(-(np.linalg.norm(x) - 1))), lambda x: 1, lambda x: np.exp(5*(-(np.linalg.norm(x) - 3))), lambda x: 1]
    #weight_functions = [lambda x: 1, lambda x: 1, lambda x: 1, lambda x: 1, lambda x: 1]

    initial_triangles_in_row = 100
    number_of_iterations = 1
    
    return curves, curves_of_surfaces, fixed_points, weight_functions, initial_triangles_in_row, number_of_iterations


if __name__ == '__main__':
    fname = 'meshes/example_1.msh'
    fnames_for_steps = ['meshes/example_1_step_1.msh', 'meshes/example_1_step_2.msh', 'meshes/example_1_step_3.msh']
    adaptation(example_1, fname, fnames_for_steps, run_gmsh=True)
    #adaptation(example_1, run_gmsh=True)
