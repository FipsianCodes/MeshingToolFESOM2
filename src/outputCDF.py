

# Data and file creation for using the cdo operator suite

import numpy as np 
import sys
from dataclasses import dataclass, field 
import matplotlib.pyplot as plt 

import src.mesher as mesher 


@dataclass
class DATASET():

    nodes : mesher.NODES
    elems : mesher.ELEMS
    centr : mesher.NODES

    v_elem : np.ndarray = field(init=False)
    v_node : np.ndarray = field(init=False)

    X : np.ndarray = field(init=False)

    griddata : np.ndarray = field(init=False)

    def test_function(self,x,y):

        alpha = 10

        x0 = self.nodes.Nd[0].x
        y0 = self.nodes.Nd[0].y

        return np.exp( -pow(x-x0,2)/(alpha) - pow(y-y0,2)/(alpha) )


    def test_node_value(self):

        N = len(self.nodes.Nd)

        self.v_node = np.zeros((N,1))
        self.X = np.zeros((N,3))

        for i in range(N):
            x = self.nodes.Nd[i].x
            y = self.nodes.Nd[i].y


            self.v_node[i] = self.test_function(x,y) 
            self.X[i][0] = x
            self.X[i][1] = y
            self.X[i][2] = self.test_function(x,y)


    def test_elem_value(self):

        pass

    def save_node_values(self,fname):

        np.savetxt(fname,self.v_node)

    def save_text_array(self,fname):

        np.savetxt(fname,self.X)

    def write_griddata(self,fname):

        with open(fname,'w') as f:
            for i in range(len(self.nodes.Nd)):
                f.write(str(self.nodes.Nd[i].x))
                f.write("  ")
                f.write(str(self.nodes.Nd[i].y))
                f.write("\n")

    def write_gridfile(self,fname):

        with open(fname,'w') as f:
            f.write("gridtype = unstructured\n")
            f.write("gridsize = ")
            f.write(str(len(self.elems.IDX)))
            f.write("\n")
            f.write("nvertex = 3\n")
            
            f.write("xvals = ")
            for i in range(len(self.centr.Nd)):
                f.write(str(self.centr.Nd[i].x))
                f.write(" ")
           
            f.write("\n")
            f.write("xbounds = ")
            for i in range(len(self.elems.IDX)):
                id0 = self.elems.IDX[i].idx0
                id1 = self.elems.IDX[i].idx1
                id2 = self.elems.IDX[i].idx2


                f.write(str(self.nodes.Nd[id0].x))
                f.write(" ")
                f.write(str(self.nodes.Nd[id1].x))
                f.write(" ")
                f.write(str(self.nodes.Nd[id2].x))
                f.write("\n          ")

            f.write("\n")
            f.write("yvals = ")
            for i in range(len(self.centr.Nd)):
                f.write(str(self.centr.Nd[i].y))
                f.write(" ")

            f.write("\n")
            f.write("ybounds = ")
            for i in range(len(self.elems.IDX)):
                id0 = self.elems.IDX[i].idx0
                id1 = self.elems.IDX[i].idx1
                id2 = self.elems.IDX[i].idx2


                f.write(str(self.nodes.Nd[id0].y))
                f.write(" ")
                f.write(str(self.nodes.Nd[id1].y))
                f.write(" ")
                f.write(str(self.nodes.Nd[id2].y))
                f.write("\n          ")
