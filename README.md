# simgraph

README

# Step 1: Geometry creation includes conversion of the 3D vascular geometry into a weighted graph.


Geometry requirement for solving the transport physics:
The format for specifying topology information of the graph is found in input.xlsx. The specifics of the column fields
provided in the excel sheet is detailed below.


# Original graph geomerty: Represented as a graph formed by nodes and edges (or vessel segments). Nodes represent terminal ends or junctions at which
n-furcation of blood vessels occur. The edges are oriented in the flow direction and each is formed by a tail (t) and head (h) node. The
inputs required to form the weighted graph are:  
 
1D example:
1--------2---------3----4
Here, 1,2,3, and 4 are the numbering of the nodes. (1,2), (2,3), (3,4) are the edges.

t - tail
h - head	
r - radius of the vessel segment formed by (tail, head)	
d - diameter of the vessel segment formed by (tail, head)	
l - length of the vessel segment formed by (tail, head)	
xpos -	x-coordinate of a node
ypos - 	y-coordinate of a node
zpos -  z-coordinate of a node	
nodes - index/numbering of the graph nodes 	
hNode - inlet node at which fluid enters the domain (e.g. 1)	
tNode - outlet node at which fluid leaves the domain (e.g. 4)	

# Discretized graph domain:

1----5----2---6---7---3----4


segment_i_nnode	index - no of intermediate(i) nodes between each head and tail nodes in the a given segment. 
nodes_mesh - index/numbering of nodes in the discretized domain
xpos_mesh - x-coordinate of nodes in the discretized domain	
ypos_mesh - y-coordinate of nodes in the discretized domain	
zpos_mesh - z-coordinate of nodes in the discretized domain
volume_ratio (computed internally, ignore)	

*****************************************************

