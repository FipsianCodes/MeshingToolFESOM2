
import numpy as np
import sys
from dataclasses import dataclass, field
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import math as m

plot_count = 0

@dataclass
class node():
    idx: int 
    x: float
    y: float
    bc: int

@dataclass
class NODES():
    Nd : node = field(init=False)
    Nd_init : bool = False

    def __post_init__(self) -> None:
        self.Nd = node(0,0,0,0)

    def add_node(self,node):
        if self.Nd_init:
            self.Nd = np.append(self.Nd,node) 
        else:
            self.Nd = node
            self.Nd_init = True

@dataclass
class elem(object):
    idx0: int  
    idx1: int  
    idx2: int  
    idx3: int  


@dataclass
class ELEMS():
    IDX: int = field(init=False)
    elem_init : bool = False

    def __post_init__(self) -> None:
        self.IDX = elem(0,0,0,0)

    def add_elem(self,elem):
        if self.elem_init:
            self.IDX = np.append(self.IDX,elem)
        else:
            self.IDX = elem
            self.elem_init = True

@dataclass
class edge():
    idx0: int 
    idx1: int

@dataclass
class EDGES():
    IDX: edge = field(init=False)
    idx_init : bool = False

    def __post_init__(self) -> None:
        self.IDX = edge(0,0)

    def add_edge(self,edge):
        if self.idx_init:
            self.IDX = np.append(self.IDX,edge)
        else:
            self.IDX = edge
            self.idx_init = True



def create_data_from_file(file_name):

    # check if input file exists
    try:
        data_check = open(file_name)
    except IOError:
        print(f"\n-->Input file '{data_file}' not accessable! \n")
        sys.exit()
    finally:
        data_check.close()

    # open file and save each content line as string in
    # a list.
    with open (file_name,"r") as query:
        data=[line.rstrip("\n") for line in query]
    return data

def convert_str_to_node(file_name):

    dat_str = create_data_from_file(file_name)
    n = len(dat_str)-1
    grid = NODES()

    dat_str.pop(0)
    for i in range(n):
        dat_split = dat_str[i].split()
        nd = node(int(dat_split[0])-1,float(dat_split[1]),float(dat_split[2]),int(dat_split[3]))
        grid.add_node(nd)
    return grid

def convert_str_to_edge(file_name):

    dat_str = create_data_from_file(file_name)
    n = len(dat_str)
    edges = EDGES()

    for i in range(n):
        dat_split = dat_str[i].split()
        ed = edge(int(dat_split[0])-1, int(dat_split[1])-1)
        edges.add_edge(ed)
    return edges

def convert_str_to_elem(file_name):

    dat_str = create_data_from_file(file_name)
    n = len(dat_str)-1
    elems = ELEMS()
    dat_str.pop(0)

    for i in range(n):
        dat_split = dat_str[i].split()
        el = elem(int(dat_split[0])-1,int(dat_split[1])-1,int(dat_split[2])-1,int(dat_split[0])-1)
        elems.add_elem(el)

    return elems

def convert_str_to_aux(file_name,nodes):
    dat_str = create_data_from_file(file_name)

    n  = len(dat_str)
    nz = int(dat_str[0])

    zlevel = np.zeros((nz,1))
    bottom = np.zeros((n-nz-1,3))

    for i in range(1,nz+1):
        zlevel[i-1] = float(dat_str[i]) 
       
    for i in range(nz+1,n):
        bottom[i-nz-1][0] = nodes.Nd[i-nz-1].x
        bottom[i-nz-1][1] = nodes.Nd[i-nz-1].y
        bottom[i-nz-1][2] = float(dat_str[i])

    return zlevel, bottom

#    np.savetxt("./bottom_topo.txt",bottom)

@dataclass
class MESH():
    nodes: NODES
    elements: ELEMS
    bottom: np.ndarray
    zlevel: np.ndarray

    edgdes: EDGES = field(init=False)

    nodes_refined: NODES = field(init=False)
    elements_refined: ELEMS = field(init=False)

    plot_count : int = 0

    boundary_nodes : NODES = field(init=False)
    internal_nodes : NODES = field(init=False)

    periodic_boundary : NODES = field(init=False)
    periodic_internal : NODES = field(init=False)
    
    aux3d_refined : np.ndarray = field(init=False)

    centroids: NODES = field(init=False)

    bottom_refined: np.ndarray = field(init=False)
    bott_ref: bool = False

    alpha : float = 50.0
    beta  : float = 15.0
    gamma : float = -90.0

    R : np.ndarray = field(init=False)

    def __post_init__(self) -> None:
        self.boundary_nodes = NODES()
        self.internal_nodes = NODES()

        self.periodic_boundary = NODES()
        self.periodic_internal = NODES()

        self.centroids = NODES()

        self.nodes_refined = NODES()
        self.elements_refined = ELEMS()

        rad = m.pi/180

        a = self.alpha  * rad
        b = self.beta * rad
        c = self.gamma * rad

        self.R = np.zeros((3,3))

        self.R[0,0] = m.cos(c) * m.cos(a) - m.sin(c) * m.cos(b) * m.sin(a)
        self.R[0,1] = m.cos(c) * m.sin(a) + m.sin(c) * m.cos(b) * m.cos(a)
        self.R[0,2] = m.sin(c) * m.sin(b)
        self.R[1,0] = -m.sin(c) * m.cos(a) - m.cos(c) * m.cos(b) * m.sin(a)
        self.R[1,1] = -m.sin(c) * m.sin(a) + m.cos(c) * m.cos(b) * m.cos(a)
        self.R[1,2] = m.cos(c) * m.sin(b)
        self.R[2,0] = m.sin(b) * m.sin(a)
        self.R[2,1] = -m.sin(b) * m.cos(a)
        self.R[2,2] = m.cos(b)

    def classify_nodes(self):

        for i in range(len(self.nodes.Nd)):
            if self.nodes.Nd[i].bc == 1:
                self.boundary_nodes.add_node(self.nodes.Nd[i])
            elif self.nodes.Nd[i].bc == 0:
                self.internal_nodes.add_node(self.nodes.Nd[i])
            elif self.nodes.Nd[i].bc == 5001:
                self.boundary_nodes.add_node(self.nodes.Nd[i])
                self.periodic_boundary.add_node(self.nodes.Nd[i])
            elif self.nodes.Nd[i].bc == 5000:
                self.internal_nodes.add_node(self.nodes.Nd[i])
                self.periodic_internal.add_node(self.nodes.Nd[i])

    def get_node_from_idx(self,index):

        return self.nodes.Nd[index]

    def set_centroids(self):

        N = len(self.elements.IDX)

        for i in range(N):
            id0 = self.elements.IDX[i].idx0 
            id1= self.elements.IDX[i].idx1 
            id2 = self.elements.IDX[i].idx2 

            nd0 = self.get_node_from_idx(id0)
            nd1 = self.get_node_from_idx(id1)
            nd2 = self.get_node_from_idx(id2)

            ax = nd0.x
            ay = nd0.y

            bx = nd1.x
            by = nd1.y

            cx = nd2.x
            cy = nd2.y

            F = 0.5*( (bx-ax)*(cy-ay) - (cx-ax)*(by-ay) )

            ltr = 50.0
            htr = 300.0
            mid = 160.0

            if F > 0.0:
                
                if ax < mid and bx > mid and cx > mid:
                    ax = ax+360.0

                elif ax < mid and bx < mid and cx > mid:
                    cx = cx-360.0
                
                elif ax < mid and bx > mid and cx < mid:
                    bx = bx-360.0
                
                elif ax > mid and bx > mid and cx < mid:
                    cx = cx+360.0
                
                elif ax > mid and bx < mid and cx < mid:
                    ax = ax-360.0

                elif ax > mid and bx < mid and cx > mid:
                    bx = bx+360.0

                xc = (ax+bx+cx)/3.0

                if xc > 360.0:
                    xc = xc-360.0
                elif xc < 0.0:
                    xc = xc+360.0

                yc = (nd0.y + nd1.y + nd2.y)/3.0
            else:
                xc = (nd0.x + nd1.x + nd2.x)/3.0
                yc = (nd0.y + nd1.y + nd2.y)/3.0

            ndc = node(i,xc,yc,0)

            self.centroids.add_node(ndc) 

    def get_centroid(self,id0,id1,id2):

        nd0 = self.get_node_from_idx(id0)
        nd1 = self.get_node_from_idx(id1)
        nd2 = self.get_node_from_idx(id2)

        xc = (nd0.x + nd1.x + nd2.x)/3.0
        yc = (nd0.y + nd1.y + nd2.y)/3.0

        ndc = node(-99,xc,yc,0) 

        return ndc


    def plot_elements(self,plot_option,switch,low_tresh,upp_tresh):

        global plot_count
        plt.figure(plot_count)
        n = len(self.elements.IDX)

        x = np.zeros((4,1))
        y = np.zeros((4,1))

        for i in range(n):

            a = self.get_node_from_idx(self.elements.IDX[i].idx0)
            x[0] = a.x
            y[0] = a.y

            a = self.get_node_from_idx(self.elements.IDX[i].idx1)
            x[1] = a.x
            y[1] = a.y

            a = self.get_node_from_idx(self.elements.IDX[i].idx2)
            x[2] = a.x
            y[2] = a.y

            a = self.get_node_from_idx(self.elements.IDX[i].idx3)
            x[3] = a.x
            y[3] = a.y

            if switch:

                if x[0] < low_tresh and x[1] < low_tresh and x[2] > upp_tresh:
                    x[2] = x[2] - 360.0
                elif x[0] < low_tresh and x[1] > upp_tresh and x[2] > upp_tresh:
                    x[1] = x[1] - 360.0
                    x[2] = x[2] - 360.0
                elif x[0] > upp_tresh and x[1] < low_tresh and x[2] < low_tresh:
                    x[0] = x[0] - 360.0
                    x[3] = x[3] - 360.0
                elif x[0] < low_tresh and x[1] > upp_tresh and x[2] < low_tresh:
                    x[1] = x[1] - 360.0
                elif x[0] > upp_tresh and x[1] < low_tresh and x[2] > upp_tresh:
                    x[0] = x[0] - 360.0
                    x[3] = x[3] - 360.0
                    x[2] = x[2] - 360.0
                elif x[0] > upp_tresh and x[1] > upp_tresh and x[2] < low_tresh:
                    x[0] = x[0] - 360.0
                    x[3] = x[3] - 360.0
                    x[1] = x[1] - 360.0
                else:
                    pass
            else:
                pass

            plt.plot(x,y,plot_option)

#        self.plot_count += 1
        plot_count += 1

    def plot_nodes(self,plot_option):
        plt.figure(self.plot_count)
        n = len(self.nodes.Nd)

        for i in range(n):
            plt.plot(self.nodes.Nd[i].x,self.nodes.Nd[i].y,plot_option)


        self.plot_count += 1

    def plot_centroids(self,plot_option):
        plt.figure(self.plot_count)
        n = len(self.centroids.Nd)

        for i in range(n):
            plt.plot(self.centroids.Nd[i].x,self.centroids.Nd[i].y,plot_option)
    
        self.plot_count += 1


    def plot_centroids_elements(self):
        plt.figure(self.plot_count)

        n = len(self.elements.IDX)

        x = np.zeros((4,1))
        y = np.zeros((4,1))

        for i in range(n):

            x[0] = self.nodes.Nd[self.elements.IDX[i].idx0].x
            y[0] = self.nodes.Nd[self.elements.IDX[i].idx0].y

            x[1] = self.nodes.Nd[self.elements.IDX[i].idx1].x
            y[1] = self.nodes.Nd[self.elements.IDX[i].idx1].y

            x[2] = self.nodes.Nd[self.elements.IDX[i].idx2].x
            y[2] = self.nodes.Nd[self.elements.IDX[i].idx2].y

            x[3] = self.nodes.Nd[self.elements.IDX[i].idx3].x
            y[3] = self.nodes.Nd[self.elements.IDX[i].idx3].y

            plt.plot(x,y,'b')

            plt.plot(self.centroids.Nd[i].x,self.centroids.Nd[i].y,'+r')

        self.plot_count += 1 

    def show_plots(self):
        plt.show()
   
    def get_mesh_subset(self,id_node):

        n = len(self.elements.IDX)

        elist = [] 
        nlist = []

        for i in range(n):
            if self.elements.IDX[i].idx0 == id_node:
                elist = np.append(elist,i)
                nlist = np.append(nlist,self.elements.IDX[i].idx0)
                nlist = np.append(nlist,self.elements.IDX[i].idx1)
                nlist = np.append(nlist,self.elements.IDX[i].idx2)
            elif self.elements.IDX[i].idx1 == id_node:
                elist = np.append(elist,i)
                nlist = np.append(nlist,self.elements.IDX[i].idx0)
                nlist = np.append(nlist,self.elements.IDX[i].idx1)
                nlist = np.append(nlist,self.elements.IDX[i].idx2)
            elif self.elements.IDX[i].idx2 == id_node:
                elist = np.append(elist,i)
                nlist = np.append(nlist,self.elements.IDX[i].idx0)
                nlist = np.append(nlist,self.elements.IDX[i].idx1)
                nlist = np.append(nlist,self.elements.IDX[i].idx2)
            else:
                pass

        nlist = np.unique(nlist)
        tlist = []

        ELEM_sub = ELEMS()
        NODE_sub = NODES()

        for i in range(len(elist)):
            ELEM_sub.add_elem(self.elements.IDX[int(elist[i])])

        for i in range(len(nlist)):
            NODE_sub.add_node(self.nodes.Nd[int(nlist[i])])
            tlist = np.append(tlist,i)

        for i in range(len(nlist)):
            NODE_sub.Nd[i].idx = int(tlist[i])

        for i in range(len(elist)):
            id0 = ELEM_sub.IDX[i].idx0
            id1 = ELEM_sub.IDX[i].idx1
            id2 = ELEM_sub.IDX[i].idx2
            id3 = ELEM_sub.IDX[i].idx3

            ID = [id0, id1, id2, id3]
            
            for l in range(len(nlist)):
               
                integer = nlist[l]
                a = [i for i,x in enumerate(ID) if x==integer] 

                if a:
                    for j in range(len(a)):
                        if a[j] == 0:
                            ELEM_sub.IDX[i].idx0 = int(tlist[l])
                        elif a[j] == 1:
                            ELEM_sub.IDX[i].idx1 = int(tlist[l])
                        elif a[j] == 2:
                            ELEM_sub.IDX[i].idx2 = int(tlist[l])
                        elif a[j] == 3:
                            ELEM_sub.IDX[i].idx3 = int(tlist[l])
                        else:
                            pass

                else:
                    pass

        return NODE_sub,ELEM_sub 
    
    def exist_node(self,x,y):

        incr = 0

        exist = False
        ii = -99

        for i in range(len(self.nodes.Nd),len(self.nodes_refined.Nd)):
            
            if self.nodes_refined.Nd[i].x == x and self.nodes_refined.Nd[i].y == y:
                exist = True
                ii   = self.nodes_refined.Nd[i].idx
                break
            else:
                pass

        return exist, ii

    def refine_congruent(self,low_tresh,upp_tresh):

        N = len(self.elements.IDX)
        
        for i in range(len(self.nodes.Nd)):
            self.nodes_refined.add_node(self.nodes.Nd[i])

        for i in range(N):
            ID = np.zeros((4,1),dtype=int)
            NC = NODES()
              
            ID0 = self.elements.IDX[i].idx0
            ID1 = self.elements.IDX[i].idx1
            ID2 = self.elements.IDX[i].idx2
            ID3 = self.elements.IDX[i].idx3

            # get coarse element nodes

            NC.add_node(self.get_node_from_idx(ID0))
            NC.add_node(self.get_node_from_idx(ID1))
            NC.add_node(self.get_node_from_idx(ID2))
            NC.add_node(self.get_node_from_idx(ID3))

            ax = NC.Nd[0].x
            ay = NC.Nd[0].y

            bx = NC.Nd[1].x
            by = NC.Nd[1].y

            cx = NC.Nd[2].x
            cy = NC.Nd[2].y

            F = 0.5*( (bx-ax)*(cy-ay) - (cx-ax)*(by-ay) )

            if F > 0: # counterclockwise triangles are reaching over the periodic boundary

                if NC.Nd[0].x < low_tresh and NC.Nd[1].x > upp_tresh:
                    xedge1 = (NC.Nd[0].x+(NC.Nd[1].x-360.0))/2.0
                    yedge1 = (NC.Nd[0].y+NC.Nd[1].y)/2.0
                elif NC.Nd[1].x < low_tresh and NC.Nd[0].x > upp_tresh:
                    xedge1 = ((NC.Nd[0].x-360.0)+NC.Nd[1].x)/2.0
                    yedge1 = (NC.Nd[0].y+NC.Nd[1].y)/2.0
                else:
                    xedge1 = (NC.Nd[0].x+NC.Nd[1].x)/2.0
                    yedge1 = (NC.Nd[0].y+NC.Nd[1].y)/2.0

                if xedge1 < 0.0:
                    xedge1 = xedge1 + 360.0
                else:
                    pass

                if NC.Nd[1].x < low_tresh and NC.Nd[2].x > upp_tresh: 
                    xedge2 = (NC.Nd[1].x+(NC.Nd[2].x-360.0))/2.0
                    yedge2 = (NC.Nd[1].y+NC.Nd[2].y)/2.0
                elif NC.Nd[2].x < low_tresh and NC.Nd[1].x > upp_tresh:
                    xedge2 = ((NC.Nd[1].x-360.0)+NC.Nd[2].x)/2.0
                    yedge2 = (NC.Nd[1].y+NC.Nd[2].y)/2.0
                else:
                    xedge2 = (NC.Nd[1].x+NC.Nd[2].x)/2.0
                    yedge2 = (NC.Nd[1].y+NC.Nd[2].y)/2.0

                if xedge2 < 0.0:
                    xedge2 = xedge2 + 360.0
                else:
                    pass

                if NC.Nd[2].x < low_tresh and NC.Nd[3].x > upp_tresh:
                    xedge3 = (NC.Nd[2].x+(NC.Nd[3].x-360.0))/2.0
                    yedge3 = (NC.Nd[2].y+NC.Nd[3].y)/2.0
                elif NC.Nd[3].x < low_tresh and NC.Nd[2].x > upp_tresh:
                    xedge3 = ((NC.Nd[2].x-360.0)+NC.Nd[3].x)/2.0
                    yedge3 = (NC.Nd[2].y+NC.Nd[3].y)/2.0
                else:
                    xedge3 = (NC.Nd[2].x+NC.Nd[3].x)/2.0
                    yedge3 = (NC.Nd[2].y+NC.Nd[3].y)/2.0

                if xedge3 < 0.0:
                    xedge3 = xedge3 + 360.0
                else:
                    pass

           

            else:
                xedge1 = (NC.Nd[0].x+NC.Nd[1].x)/2.0
                yedge1 = (NC.Nd[0].y+NC.Nd[1].y)/2.0
        
                xedge2 = (NC.Nd[1].x+NC.Nd[2].x)/2.0
                yedge2 = (NC.Nd[1].y+NC.Nd[2].y)/2.0

                xedge3 = (NC.Nd[2].x+NC.Nd[3].x)/2.0
                yedge3 = (NC.Nd[2].y+NC.Nd[3].y)/2.0
                

           
            vec_aux = np.zeros((1,3))
            # check the first edge for esixiting refinement

            check, ii = self.exist_node(xedge1,yedge1)

            if check:
                id0 = ii
            else:
                if NC.Nd[0].bc == 1 and NC.Nd[1].bc == 1:
                    nc0 = node(len(self.nodes_refined.Nd),xedge1,yedge1,1)
                elif (NC.Nd[0].bc == 5001 and NC.Nd[1].bc == 1) or (NC.Nd[0].bc == 1 and NC.Nd[1].bc == 5001):
                    nc0 = node(len(self.nodes_refined.Nd),xedge1,yedge1,1)
                elif (NC.Nd[0].bc == 5000 and NC.Nd[1].bc == 5000):
                    nc0 = node(len(self.nodes_refined.Nd),xedge1,yedge1,5000)
                elif (NC.Nd[0].bc == 5001 and NC.Nd[1].bc == 5000) or (NC.Nd[0].bc == 5000 and NC.Nd[1].bc == 5001):
                    nc0 = node(len(self.nodes_refined.Nd),xedge1,yedge1,5000)
                else:
                    nc0 = node(len(self.nodes_refined.Nd),xedge1,yedge1,0)
                id0 = len(self.nodes_refined.Nd)
                self.nodes_refined.add_node(nc0)

                vec_aux[0][0] = xedge1 
                vec_aux[0][1] = yedge1  
 
                v0 = self.bottom[ID0][2]
                v1 = self.bottom[ID1][2]
                
                if v1 > v0:
                    vec_aux[0][2] = v0
                elif v0 > v1:
                    vec_aux[0][2] = v1
                else:
                    vec_aux[0][2] = v0


#                print("V0 :: ",self.bottom[ID0][2],"   V1 :: ",self.bottom[ID1][2],"   vAux :: ",vec_aux[0][2])
#                vec_aux[0][2] = round((self.bottom[ID0][2] + self.bottom[ID1][2])/2.0)

                if self.bott_ref:
                    self.bottom_refined = np.append(self.bottom_refined,vec_aux,axis=0)
                else:
                    self.bott_ref = True
                    self.bottom_refined = vec_aux

          
            # check the second edge for esixiting refinement

            check, ii = self.exist_node(xedge2,yedge2)

            if check:
                id1 = ii
            else:
                if NC.Nd[1].bc == 1 and NC.Nd[2].bc == 1:
                    nc1 = node(len(self.nodes_refined.Nd),xedge2,yedge2,1)
                elif (NC.Nd[1].bc == 5001 and NC.Nd[2].bc == 1) or (NC.Nd[1].bc == 1 and NC.Nd[2].bc == 5001):
                    nc1 = node(len(self.nodes_refined.Nd),xedge2,yedge2,1)
                elif (NC.Nd[1].bc == 5000 and NC.Nd[2].bc == 5000):
                    nc1 = node(len(self.nodes_refined.Nd),xedge2,yedge2,5000)
                elif (NC.Nd[1].bc == 5001 and NC.Nd[2].bc == 5000) or (NC.Nd[1].bc == 5000 and NC.Nd[2].bc == 5001):
                    nc1 = node(len(self.nodes_refined.Nd),xedge2,yedge2,5000)
                else:
                    nc1 = node(len(self.nodes_refined.Nd),xedge2,yedge2,0)
                id1 = len(self.nodes_refined.Nd)
                self.nodes_refined.add_node(nc1)

                vec_aux[0][0] = xedge2 
                vec_aux[0][1] = yedge2  
                
                v1 = self.bottom[ID1][2]
                v2 = self.bottom[ID2][2]
                
                if v2 > v1:
                    vec_aux[0][2] = v1
                elif v1 > v2:
                    vec_aux[0][2] = v2
                else:
                    vec_aux[0][2] = v1


#                print("V1 :: ",self.bottom[ID1][2],"   V2 :: ",self.bottom[ID2][2],"   vAux :: ",vec_aux[0][2])
#                vec_aux[0][2] = round((self.bottom[ID1][2] + self.bottom[ID2][2])/2.0)
                
                self.bottom_refined = np.append(self.bottom_refined,vec_aux,axis=0)
                
           
            # check the third edge for esixiting refinement

            check, ii = self.exist_node(xedge3,yedge3)

            if check:
                id2 = ii
            else:
                if NC.Nd[2].bc == 1 and NC.Nd[3].bc == 1:
                    nc2 = node(len(self.nodes_refined.Nd),xedge3,yedge3,1)
                elif (NC.Nd[2].bc == 5001 and NC.Nd[3].bc == 1) or (NC.Nd[2].bc == 1 and NC.Nd[3].bc == 5001):
                    nc2 = node(len(self.nodes_refined.Nd),xedge3,yedge3,1)
                elif (NC.Nd[2].bc == 5000 and NC.Nd[3].bc == 5000):
                    nc2 = node(len(self.nodes_refined.Nd),xedge3,yedge3,5000)
                elif (NC.Nd[2].bc == 5001 and NC.Nd[3].bc == 5000) or (NC.Nd[2].bc == 5000 and NC.Nd[3].bc == 5001):
                    nc2 = node(len(self.nodes_refined.Nd),xedge3,yedge3,5000)
                else:
                    nc2 = node(len(self.nodes_refined.Nd),xedge3,yedge3,0)
                id2 = len(self.nodes_refined.Nd)
                self.nodes_refined.add_node(nc2)
                
                vec_aux[0][0] = xedge3 
                vec_aux[0][1] = yedge3  

                v2 = self.bottom[ID2][2]
                v3 = self.bottom[ID3][2]

                if v3 > v2:
                    vec_aux[0][2] = v2
                elif v2 > v3:
                    vec_aux[0][2] = v3
                else:
                    vec_aux[0][2] = v2


#                print("V2 :: ",self.bottom[ID2][2],"   V3 :: ",self.bottom[ID3][2],"   vAux :: ",vec_aux[0][2])
#                vec_aux[0][2] = round((self.bottom[ID2][2] + self.bottom[ID3][2])/2.0)
                self.bottom_refined = np.append(self.bottom_refined,vec_aux,axis=0)

        

            elem1 = elem(ID0,id0,id2,ID0)
            elem2 = elem(ID1,id1,id0,ID1)
            elem3 = elem(ID2,id2,id1,ID2)
            elem4 = elem(id0,id1,id2,id0)

            self.elements_refined.add_elem(elem1)
            self.elements_refined.add_elem(elem2)
            self.elements_refined.add_elem(elem3)
            self.elements_refined.add_elem(elem4)



    def refine_centroid(self):

        N = len(self.nodes.Nd)

        for i in range(N):
            self.nodes_refined.add_node(self.nodes.Nd[i])


        for i in range(len(self.elements.IDX)):
            i0 = self.elements.IDX[i].idx0
            i1 = self.elements.IDX[i].idx1
            i2 = self.elements.IDX[i].idx2
            i3 = self.elements.IDX[i].idx3

#            nc = self.get_centroid(i0,i1,i2) 
            nc = self.centroids.Nd[i]

            nc.idx = len(self.elements.IDX)+i+1

            self.nodes_refined.add_node(nc)

            tri = elem(i0,i1,nc.idx,i0)
            self.elements_refined.add_elem(tri)

            tri = elem(i1,i2,nc.idx,i1)
            self.elements_refined.add_elem(tri)

            tri = elem(i2,i3,nc.idx,i2)
            self.elements_refined.add_elem(tri)

    def return_refined_mesh(self):

        return self.nodes_refined, self.elements_refined
        
    def refine_aux3d(self,itype):

        if itype == "spline":

            from scipy.interpolate import CloughTocher2DInterpolator

            x = self.bottom[:,0]
            y = self.bottom[:,1]
            z = self.bottom[:,2]

            x = np.append(x,self.bottom[:,0]-360)
            y = np.append(y,self.bottom[:,1])
            z = np.append(z,self.bottom[:,2])

            x = np.append(x,self.bottom[:,0]+360)
            y = np.append(y,self.bottom[:,1])
            z = np.append(z,self.bottom[:,2])

            interpolator = CloughTocher2DInterpolator(np.array([x,y]).T,z)


            bottom_refined = np.zeros((len(self.nodes_refined.Nd),3))

            for i in range(0,len(bottom_refined)):
                bottom_refined[i,0] = self.nodes_refined.Nd[i].x
                bottom_refined[i,1] = self.nodes_refined.Nd[i].y
                bottom_refined[i,2] = interpolator(self.nodes_refined.Nd[i].x,self.nodes_refined.Nd[i].y)

#                print(bottom_refined[i,0],"  ",bottom_refined[i,1],"  ",bottom_refined[i,2])

            self.aux3d_refined = bottom_refined

        elif itype == "linear":

            self.aux3d_refined = np.append(self.bottom,self.bottom_refined,axis=0) 


        elif itype == "hottub":

            bottom_refined = np.zeros((len(self.nodes_refined.Nd),3))

            for i in range(0,len(bottom_refined)):

                bottom_refined[i,0] = self.nodes_refined.Nd[i].x
                bottom_refined[i,1] = self.nodes_refined.Nd[i].y
                if self.nodes_refined.Nd[i].bc == 1:
                    bottom_refined[i,2] = -30.0
                else:
                    bottom_refined[i,2] = -6000.0

#            self.aux3d_refined = np.append(self.bottom,self.bottom_refined,axis=0)
            self.aux3d_refined = bottom_refined

            np.savetxt("./auxfine.txt",self.aux3d_refined)

        else:
            print("Unknown interpolation type! No refinement possible ...\n")
            pass

#        np.savetxt("aux3d_refined",bottom_refined)

    def write_refined_mesh(self,fpath):

        fnode = fpath + 'nod2d.out'
        felem = fpath + 'elem2d.out'
        faux3 = fpath + 'aux3d.out'

        with open(fnode,'w') as f:
            f.write(str(len(self.nodes_refined.Nd)))
            f.write("\n")
           
            for i in range(len(self.nodes_refined.Nd)):
                f.write(str(i+1))
                f.write(" ")
                f.write(str(self.nodes_refined.Nd[i].x))
                f.write(" ")
                f.write(str(self.nodes_refined.Nd[i].y))
                f.write(" ")
                f.write(str(self.nodes_refined.Nd[i].bc))
                f.write("\n")


        with open(felem,'w') as f:
            f.write(str(len(self.elements_refined.IDX)))
            f.write("\n")
           
            for i in range(len(self.elements_refined.IDX)):
                f.write(str(self.elements_refined.IDX[i].idx0+1))
                f.write(" ")
                f.write(str(self.elements_refined.IDX[i].idx1+1))
                f.write(" ")
                f.write(str(self.elements_refined.IDX[i].idx2+1))
                f.write("\n")

        with open(faux3,'w') as f:
            f.write(str(len(self.zlevel)))
            f.write("\n")

            for i in range(len(self.zlevel)):
                f.write(str(self.zlevel[i][0]))
                f.write("\n")

            for i in range(len(self.aux3d_refined)):
                f.write("   ")
                f.write(str(self.aux3d_refined[i][2]))
                f.write("\n")


    def save_nodes_to_file(self,fname):

        with open(fname,'w') as f:
            for i in range(len(self.nodes.Nd)):
                f.write(str(self.nodes.Nd[i].idx))
                f.write(" ")
                f.write(str(self.nodes.Nd[i].x))
                f.write(" ")
                f.write(str(self.nodes.Nd[i].y))
                f.write(" ")
                f.write(str(self.nodes.Nd[i].bc))
                f.write("\n")

    def save_boundary_nodes_to_file(self,fname):

        with open(fname,'w') as f:
            for i in range(len(self.boundary_nodes.Nd)):
                f.write(str(self.boundary_nodes.Nd[i].idx))
                f.write(" ")
                f.write(str(self.boundary_nodes.Nd[i].x))
                f.write(" ")
                f.write(str(self.boundary_nodes.Nd[i].y))
                f.write(" ")
                f.write(str(self.boundary_nodes.Nd[i].bc))
                f.write("\n")

    def check_orientation(self):

        cw = 0
        ccw = 0

        N = len(self.elements.IDX)

        for i in range(N):

            ID = np.zeros((3,1),dtype=int)
            NC = NODES()

            ID0 = self.elements.IDX[i].idx0
            ID1 = self.elements.IDX[i].idx1
            ID2 = self.elements.IDX[i].idx2

            # get coarse element nodes

            NC.add_node(self.get_node_from_idx(ID0))
            NC.add_node(self.get_node_from_idx(ID1))
            NC.add_node(self.get_node_from_idx(ID2))

            ax = NC.Nd[0].x
            ay = NC.Nd[0].y

            bx = NC.Nd[1].x
            by = NC.Nd[1].y

            cx = NC.Nd[2].x
            cy = NC.Nd[2].y

            F = 0.5*( (bx-ax)*(cy-ay) - (cx-ax)*(by-ay) )

            if F > 0:
                ccw = ccw+1
            else:
                cw = cw+1

        print(N," elements with ",cw," clockwise and ",ccw," counterclockwise triangles\n")

    def rot_node(self,lon,lat,inv):

        rad = m.pi/180

        lat = lat * rad
        lon = lon * rad

        x = m.cos(lat) * m.cos(lon)
        y = m.cos(lat) * m.sin(lon)
        z = m.sin(lat)

        if inv:
            T = np.transpose(self.R)
        else:
            T = self.R

        xrot = T[0,0]*x + T[0,1]*y + T[0,2]*z
        yrot = T[1,0]*x + T[1,1]*y + T[1,2]*z
        zrot = T[2,0]*x + T[2,1]*y + T[2,2]*z

        latrot = m.asin(zrot) / rad
        lonrot = m.atan2(yrot,xrot) / rad


        return lonrot, latrot

    def rotate_nodes(self,inv):

        for i in range(len(self.nodes.Nd)):

            lat = self.nodes.Nd[i].y
            lon = self.nodes.Nd[i].x

            lonrot,latrot = self.rot_node(lon,lat,inv)
            
            if lonrot < 0:
                lonrot = 360.0+lonrot
            else:
                pass

            self.nodes.Nd[i].y = latrot
            self.nodes.Nd[i].x = lonrot

    def rotate_refined_nodes(self,inv):

        for i in range(len(self.nodes_refined.Nd)):

            lat = self.nodes_refined.Nd[i].y
            lon = self.nodes_refined.Nd[i].x

            lonrot,latrot = self.rot_node(lon,lat,inv)

            if lonrot < 0:
                lonrot = 360.0+lonrot
            else:
                pass

            self.nodes_refined.Nd[i].y = latrot
            self.nodes_refined.Nd[i].x = lonrot

   





        

