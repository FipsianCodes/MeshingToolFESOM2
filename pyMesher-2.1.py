#!/usr/bin/env python3

import src.mesher as mesher
import src.outputCDF as cdf

def main() -> None:
    
    # define path to PI mesh 2.1

    str_elem = "$PATH_TO_PIMESH/pi-mesh-2.1/pi-mesh-2.1/elem2d.out"
    str_node = "$PATH_TO_PIMESH/pi-mesh-2.1/pi-mesh-2.1/nod2d.out"
    str_aux  = "$PATH_TO_PIMESH/pi-mesh-2.1/aux3d.out"

    # read in files
    grid = mesher.convert_str_to_node(str_node)

    elements = mesher.convert_str_to_elem(str_elem)    

    [zlevel, bot_topo] = mesher.convert_str_to_aux(str_aux,grid)

    # compile PI mesh information 
    pimesh = mesher.MESH(grid,elements,bot_topo,zlevel)

    # if mesh needs rotation, change to true 
    pimesh.rotate_nodes(False)

    pimesh.refine_congruent(100,250)

    # linear interpolation of the bathymetry
    pimesh.refine_aux3d("linear")

    # rotation of the fined nodes 
    pimesh.rotate_refined_nodes(True)

    # write output of the refined mesh to working directory
    pimesh.write_refined_mesh("./")

    # additional option. needed for the gridfiles when conservative interpolation is required

    # pimesh.set_centroids()
    # pimesh_grid = cdf.DATASET(pimesh.nodes,pimesh.elements,pimesh.centroids)
    # pimesh_grid.write_gridfile("./pimesh_gridfile.txt")


    # additional option. Writes gridfiles of the refined mesh.

    # fineMesh = mesher.MESH(pimesh.nodes_refined,pimesh.elements_refined,pimesh.aux3d_refined,zlevel)
    # fineMesh.set_centroids()
    # fpimesh_grid = cdf.DATASET(fineMesh.nodes,fineMesh.elements,fineMesh.centroids)
    # fpimesh_grid.write_gridfile("./fpimesh_gridfile.txt")

    # not optimized. recommended for small test cases only
    #pimesh.show_plots()

if __name__ == "__main__":
    main()




