
# README

## Preprocessing: Preparing the input geometry

Geometry required for solving the transport physics involves conversion of the 3D vascular geometry into a weighted graph. The format for specifying topology information of the graph is found in [input.xlsx](https://github.com/DeepaMahm/simgraph/tree/main/input).
The specifics of the column fields provided in the excel sheet is detailed below.

*****************************************************

#### Step 1:   Original graph geomerty
Represented as a graph formed by nodes and edges (or vessel segments). Nodes represent terminal ends or junctions at which
n-furcation of blood vessels occur. The edges are oriented in the flow direction and each is formed by a tail (t) and head (h) node. The
inputs required to form the weighted graph are:  
 
Example: <br />
![A test image](/docs/images/test2.png)

<!-- <img src="https://github.com/DeepaMahm/simgraph/blob/main/test2.png" width="700"> -->

Here, 1, 2, 3, and 4 are the nodes labels. (1,2), (2,3), (3,4) are the edges.

The following 11 fileds are the column headers present in the spreadsheet. These values in these fileds are extracted from the skeletonized
version of the vascular geometry.<br />
* t - tail <br />
* h - head	 <br />
* r - radius of the vessel segment formed by (tail, head)	 <br />
* d - diameter of the vessel segment formed by (tail, head)	 <br />
* l - length of the vessel segment formed by (tail, head)	 <br />
* xpos -	x-coordinate of a node <br />
* ypos - 	y-coordinate of a node <br />
* zpos -  z-coordinate of a node	 <br />
* nodes(index) - index/numbering of the graph nodes  <br /> 	
* hNode - inlet node at which fluid enters the domain (e.g. 1)  <br />	
* tNode - outlet node at which fluid leaves the domain (e.g. 4)	 <br />

*****************************************************

#### Step 2: Discretizing graph domain

* nodes_mesh - index/numbering of nodes in the discretized domain  <br />
* xpos_mesh - x-coordinate of nodes in the discretized domain	 <br />
* ypos_mesh - y-coordinate of nodes in the discretized domain	 <br />
* zpos_mesh - z-coordinate of nodes in the discretized domain <br />
* segment_i_nnode	- no of intermediate(i) nodes between each head and tail nodes in the a given segment. e.g if 2 nodes
 (labelled 5 and 6) are added between (1, 2) segment_i_nnode=2. <br />

* *volume_ratio and other extras (computed internally, not user defined)*	 <br />

For carring out simulations with new vasculature, add a new workbook and populate the topology information corresponding to the 11 fields mentioned in Step 1.
To generate the data of the column fields mentioned in Step 2, use [read_mesh.py](https://github.com/DeepaMahm/simgraph/blob/main/preprocessing/read_mesh.py) for domain discretization in Gmsh. The output saved in [data_model.xlsx](https://github.com/DeepaMahm/simgraph/blob/main/preprocessing/data_model.xlsx) can be updated in 
the input spreadsheet prior to executing main.exe.

1. Run extract_subgraph to retain just one inlet and one outlet in the user-defined graph
1. Update the network topology in input.xlsx 
1. Run read_mesh.py to generate the coordinates of the mesh elements by specifying the (minimum, maximum) bounds of element size. 
1. The output file with information fo the 
1. Update the pos_mesh coordinates in input.xlsx
1. The geomerty is ready; run static and dynamic simulations

*****************************************************
## Simulation
*Generation of the metamodel, mathematical formulation and simulation is carried out in [main.exe](https://github.com/DeepaMahm/simgraph/blob/main/main/main.exe)*

#### Step 3: Setting up the computational domain for multiphysics simulation
Each node in the blood vessel that is modelled is an approximation of a volume element. The volume is computed considering a cylindrial geometry. The height of the cylinder(l) is equal to mesh size. The diameter of the cylinder (d) is equal to the diameter of the blood vessel that is modelled. The blood vessel volume element is surrounded by a monolayer of β-cell mass. The volum of this monolayer is approximated by the volume of a hollow cylinder also referred to as donut volume.  The idea is to treat 
each volume element like a continuous stirred tank reactor (CSTR). Convection , diffusion and reaction is modeled by writing mole balance equation for each volume element. The reaction nodes are highlighted in  green on the left image. Coupling between the donut element and the cylindrical volume element in the blood occurs at the green colored nodes. The bifurcation and trifurcation junctions are treated as spherical volume elements.
 
<!--  ![A test image](domain.svg) -->
<img src="https://github.com/DeepaMahm/simgraph/blob/main/docs/images/domain.svg" width="700">

![image](https://user-images.githubusercontent.com/29662579/128639317-29a7b18b-4b1b-433a-b042-53aadef1e4bc.png)

#### Step 4: Simulation of static flow fields and transient concentration dynamics 



## Postprocessing


Output files:

* File naming convention for mat files.*

e.g test1_default_bc1_v1_c2

{graph}_{mesh}_{boundary condition}_{no. of species in blood vessel}_{no. of species in cell}


Naming convection used: 
advection diffusion: {species}_{tool}_{physics}

species =  glc_ext/lac_ext/glc_cell/lac_cell
tool = simgraph/comsol
physics = ad/sink/' '/eq

e.g.
*  ad(advection diffusion) :  glc_ext_comosl_ad
*  sink( irreversible) : glc_ext_comsol_sink
*  ' ' (convection + diffusion + reaction): glc_ext_comsol
* 	eq (convection + diffusion + reaction - blood vessel and cell have equal volumes): glc_ext_comsol_eq 
* 	eqfull (convection + diffusion + reaction - blood vessel and cell have equal volumes): glc_ext_comsol_eqfull 
*  full (convection + diffusion + reaction - blood vessel and cell have unequal volumes): glc_ext_comsol_full 
