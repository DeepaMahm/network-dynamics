# -----------------------------------------------------------------------------
#
#  Gmsh Python extended tutorial 1
#
#  Geometry and mesh data
#  Reference:
# https://mail.google.com/mail/u/0/#search/mesh/KtbxLwhGPJBwdSHgPpGvkDnkxlXPmdPpdV
# https://gitlab.onelab.info/gmsh/gmsh/blob/gmsh_4_7_1/tutorial/python/x1.py
# https://forums.autodesk.com/t5/autocad-forum/adding-points-to-an-existing-geometry/m-p/9505285#M1019586
# https://gitlab.onelab.info/gmsh/gmsh/-/issues/1152 - coherence commands
# ------------------------------------------------------------------------------
# Note to self:
# While creating mesh from gmsh GUI using iges file exported from AutoCAD as input, make sure the units in iges file are
# are in micrometers and not in centimeters (option 9, 2HUM and not 10, 2HCM)

# STEPS:
# manual : trial0
# 1. create dxf file in cad
# 2. import dxf in comsol and create mesh by specifying min/max element size
# 3. export the mesh points from Mesh node in comsol to mphtxt
# 4. copy and save mesh points from mphtxt to txt file
# 5. In autocad, add points from step 4 to the geometry via Points option in autocad; save dxf
# 6. import dxf from step 5 into comsol and export the geometry from Geometry node in .step format (issue with iges)
# 7. import step in gmsh's GUI, remove duplicate nodes via Coherence -> save .geo; click mesh -> 1D to create elements
# 8. save .msh format
# 9. read .msh contents via meshio, generate graph
# 10. Note: directly exporting iges from cad via IGESEXPORT command doesn't help; not possible to generate elements


# Automated:
# 1. run extract_subgraph to retain just one inlet and one outlet in the user-defined graph
# 2. update the network topology in input.xlsx
# 3. run the untwiched version of the code to find out the right direction of flow (velocity computation) in MATLAB
# 4. update the edges with the right direction of flow in input.xlsx
# 5. run read_mesh.py to generate the coordinates of the mesh elements
# 6. update the pos_mesh coordinates in input.xlsx
# 7. run concentration simulations

# Note to self: Interpreter : python 3.6
# ------------------------------------------------------------------------------

import sys
from pathlib import Path

import gmsh
import pygmsh
import meshio
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

from vedo import *
# from vedo import io
from collections import OrderedDict
from matpancreas.settings_model import GRAPH_MESH_FILE
from matpancreas.utils.graph_utils_py import draw_graph3d_vedo, convert_graph_to_df, get_eucledian_lengths, \
    convert_coordinate_units, create_graph, get_eucledian_lengths_wrt_origin, remove_edge_attribute, get_graph
from matpancreas.utils.io_utils_py import write_output


def mesh_network(G):

    """
    loads a networkx graph and meshed the edges (/lines) to add intermediate nodes via gmsh library
    References:
    ----------
    https://stackoverflow.com/a/431868/8281509
    """
    # draw_graph3d_vedo([G], label_scale=5) #G=G, point_r=1, lw=1)
    H = G.copy()
    H = remove_edge_attribute(H, edge_attr=['l'])
    pos = nx.get_node_attributes(G, 'pos')

    # generate mesh
    gmsh.initialize()
    p = [gmsh.model.geo.addPoint(*point, tag=tag) for tag, point in pos.items()]
    gmsh.model.geo.synchronize()
    # print(gmsh.model.getEntities())

    for i, e in enumerate(sorted(G.edges())):
        t, h = e
        gmsh.model.geo.addLine(t, h, tag=i)
    gmsh.model.geo.synchronize()

    gmsh.option.setNumber("Mesh.MeshSizeMin", 9.5)
    gmsh.option.setNumber("Mesh.MeshSizeMax", 13.5)
    gmsh.model.mesh.generate(1)  # make a 1d mesh
    gmsh.write('test.msh')

    for i, e in enumerate(sorted(G.edges(data=True))):
        t, h, attr = e
        node_tags, node_coords = gmsh.model.mesh.getNodes(dim=1, tag=i)[0:2]
        # node_coords = np.reshape(node_coords, len(node_coords))
        nx.add_path(H, nodes_for_path=[t, *node_tags, h], d=attr['d'])
        nx.set_edge_attributes(G, {(t, h): len(node_tags)}, 'segment_i_node')
        # initial edge is removed when meshing creates interior nodes
        if len(node_tags) > 0:
            H.remove_edge(t, h)
    # ------------------------------------------------------------------------------------------------------------------
    # set positions of newly added mesh nodes and edge lengths
    # ------------------------------------------------------------------------------------------------------------------
    mesh_nodes, mesh_coords = gmsh.model.mesh.getNodes(dim=1)[0:2]
    mesh_coords = np.reshape(mesh_coords, (len(mesh_nodes), 3))
    mesh_pos = dict(zip(mesh_nodes, mesh_coords))

    nx.set_node_attributes(H, mesh_pos, 'pos')
    H = get_eucledian_lengths(H)
    gmsh.finalize()

    return G, H


if __name__ == '__main__':

    G = create_graph(attributes=True, actual_pos=False, directed=True)
    exit()
    # draw_graph3d_vedo([G], label_scale=3)
    # G.remove_edges_from(
    #     ebunch=[(289, 303), (291, 337), (339, 329),
    #             (14, 12), (229, 308), (168, 171),
    #             (323, 251), (144, 145), (68, 83),
    #             (32, 30), (172, 174), (314, 338),
    #             (178, 180), (88, 86), (159, 11)]
    # )

    G = get_eucledian_lengths(G)
    G, H = mesh_network(G=G)
    graph_df = convert_graph_to_df(H)
    write_output(d={'islet': graph_df}, file='data_model')
    # draw_graph3d_vedo([H], label_scale=3) #G=G, point_r=1, lw=1)
