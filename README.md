
README

# Geometry creation:
Includes the steps involved in conversion of the 3D vascular geometry into a weighted graph.

Geometry requirement for solving the transport physics: The format for specifying topology information of the graph is found in input.xlsx.
The specifics of the column fields provided in the excel sheet is detailed below.

*****************************************************

#### Step 1:   Original graph geomerty:
Represented as a graph formed by nodes and edges (or vessel segments). Nodes represent terminal ends or junctions at which
n-furcation of blood vessels occur. The edges are oriented in the flow direction and each is formed by a tail (t) and head (h) node. The
inputs required to form the weighted graph are:  
 
1D example: <br />
1--------2---------3----4  <br />
Here, 1,2,3, and 4 are the numbering of the nodes. (1,2), (2,3), (3,4) are the edges.

The following 11 fileds are the column headers present in the spreadsheet. These fileds are extracted from the skeletonized
version of the vascular geometry.<br />
* t - tail <br />
* h - head	 <br />
* r - radius of the vessel segment formed by (tail, head)	 <br />
* d - diameter of the vessel segment formed by (tail, head)	 <br />
* l - length of the vessel segment formed by (tail, head)	 <br />
* xpos -	x-coordinate of a node <br />
* ypos - 	y-coordinate of a node <br />
* zpos -  z-coordinate of a node	 <br />
* nodes - index/numbering of the graph nodes  <br /> 	
* hNode - inlet node at which fluid enters the domain (e.g. 1)  <br />	
* tNode - outlet node at which fluid leaves the domain (e.g. 4)	 <br />

*****************************************************

#### Step 2: Discretized graph domain:

1----5----2---6---7---3----4


* segment_i_nnode	index - no of intermediate(i) nodes between each head and tail nodes in the a given segment.  <br />
* nodes_mesh - index/numbering of nodes in the discretized domain  <br />
* xpos_mesh - x-coordinate of nodes in the discretized domain	 <br />
* ypos_mesh - y-coordinate of nodes in the discretized domain	 <br />
* zpos_mesh - z-coordinate of nodes in the discretized domain <br />
* *volume_ratio (computed internally, not user defined)*	 <br />

For carring out simulations with new vasculature, add a new workbook and populate the topology information corresponding to the 11 fields mentioned in Step 1.
To generate the data of the column fields mentioned in Step 2, use read_mesh.py for domain discretization in Gmsh.
1. Run extract_subgraph to retain just one inlet and one outlet in the user-defined graph
1. Update the network topology in input.xlsx 
1. Run read_mesh.py to generate the coordinates of the mesh elements by specifying the (minimum, maximum) bounds of element size. 
1. The output file with information fo the 
1. Update the pos_mesh coordinates in input.xlsx
1. The geomerty is ready; run static and dynamic simulations

*****************************************************


