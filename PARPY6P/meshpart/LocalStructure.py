#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 23:10:58 2020

@author: kissami
"""
from mpi4py import MPI
from math import sqrt
import numpy as np
import meshio
from collections import OrderedDict 
import matplotlib.pyplot as plt
import matplotlib

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

def local_mesh():
    Cells = []
    Nodes = []
       
    file = open("mesh"+str(rank)+".txt","r")
    for line in file:
        #read elements
        if line == "Elements\n":
            continue   
        if (line == "EndElements\n"):
            continue
        if line == "Nodes\n":
            break
        Cells.append([int(x) for x in line.split()])   
    for line in file:
        #read Nodes
        if line == "Nodes\n":
            continue   
        if (line == "EndNodes\n"):
            continue
        if line == "HalosInt\n":
            break
        Nodes.append([float(x) for x in line.split()])
    size = 2    
    if size > 1:
        
        halos = []
        haloe = []
        gtol = OrderedDict()#[]
        ltog = OrderedDict()#[]
        neigh = []
        
        for line in file:
            #read halos
            if line == "HalosInt\n":
                continue   
            if (line == "EndHalosInt\n"):
                break
            halos.append(int(line))#[int(x) for x in line.split()])  
        
        for line in file:
            #read halos
            if line == "HalosExt\n":
                continue   
            if (line == "EndHalosExt\n"):
                break
            haloe.append([int(x) for x in line.split()])  
        
        cmpt = 0
        for line in file:
            #read Global cell to local
            if line == "GlobalCellToLocal\n":
                continue   
            if (line == "EndGlobalCellToLocal\n"):
                break
            gtol[int(line)] = cmpt
            cmpt +=1
            #gtol[int(line)]#[int(x) for x in line.split()])  
        
        cmpt = 0
        for line in file:
            #read Local Node To Global
            if line == "LocalNodeToGlobal\n":
                continue   
            if (line == "EndLocalNodeToGlobal\n"):
                break
            ltog[cmpt] = int(line)
            cmpt +=1
            #ltog.append(int(line))# for x in line.split()])  
        
        for line in file:
            #read LocalToGlobal
            if line == "Neigh\n":
                continue   
            if (line == "EndNeigh\n"):
                break
            neigh.append([int(x) for x in line.split()])#int(line))
        
     
    if size > 1:
        return Cells, Nodes, halos, haloe, gtol, ltog, neigh
    else:
        return Cells, Nodes


def create_Faces_information(Cells, Nodes, HaloExt, GlobalIndexNode):
    Faces = []
    cellF = []
    k = 0
    for i in range(len(Cells)):
        Faces.append([Cells[i][0], Cells[i][1]])
        Faces.append([Cells[i][1], Cells[i][2]])
        Faces.append([Cells[i][0], Cells[i][2]])
        cellF.append([Faces[k], Faces[k+1], Faces[k+2]])
        k = k+3
    
    Faces =   set( tuple(x) for x in Faces)
    Faces = list(Faces)
    
    Facesdict = OrderedDict()#{Faces[0]: 0}
    for i in range(len(Faces)):
        Facesdict[Faces[i]] = i#Faces[i]
    
    # create faceid of Cells
    cell_faceid = []
    for i in range(len(cellF)):
        cell_faceid.append([Facesdict.get(tuple(cellF[i][0])), Facesdict.get(tuple(cellF[i][1])), 
                            Facesdict.get(tuple(cellF[i][2]))])

    face_cellid = [(-1, -1)]*len(Faces)
    for i in range(len(Cells)):
        for j in range(3):
            if face_cellid[cell_faceid[i][j]] == (-1, -1):
                face_cellid[cell_faceid[i][j]]= (i, -1)
            if face_cellid[cell_faceid[i][j]][0] != i:
                face_cellid[cell_faceid[i][j]]= (face_cellid[cell_faceid[i][j]][0], i)
    
    haloFaces = []
    k = 0
    for i in range(len(HaloExt)):
        haloFaces.append([HaloExt[i][0], HaloExt[i][1]])
        haloFaces.append([HaloExt[i][1], HaloExt[i][2]])
        haloFaces.append([HaloExt[i][0], HaloExt[i][2]])
        k = k+3
        
    haloFacesdict = OrderedDict() #{tuple(haloFaces[0]) : 0}
    for i in range(len(haloFaces)):
        haloFacesdict[tuple(haloFaces[i])] = i
    
    for i in range(len(Faces)):
        if  haloFacesdict.get(tuple([GlobalIndexNode[Faces[i][0]], GlobalIndexNode[Faces[i][1]]])):
            face_cellid[i] = (face_cellid[i][0], -10)
            
    face_name = [0]*len(Faces)
    for i in range(len(Faces)):
        if (face_cellid[i][1] == -1):
            if (Nodes[Faces[i][0]][3] == Nodes[Faces[i][1]][3] ):
                face_name[i] = Nodes[Faces[i][0]][3]
            if ((Nodes[Faces[i][0]][3] == 1 and Nodes[Faces[i][1]][3] !=0) or 
                (Nodes[Faces[i][0]][3] != 0 and Nodes[Faces[i][1]][3] ==1)):
                face_name[i] = 1
            if ((Nodes[Faces[i][0]][3] == 2 and Nodes[Faces[i][1]][3] !=0) or 
                (Nodes[Faces[i][0]][3] != 0 and Nodes[Faces[i][1]][3] ==2)):
                face_name[i] = 2
        if (face_cellid[i][1] == -10):
            face_name[i] = 10
                
                
    return Faces, face_cellid, cell_faceid , face_name, haloFacesdict



def ExplicitScheme(w, Cells, Nodes, Cell_faceid, Faces, Face_cellid, Face_name, HaloFaces, HaloInt, Neigh,
                   GlobalIndexCell, GlobalIndexNode):
    
#    def ghost_value(w, face_cellid, face_name):
#        w_ghost = OrderedDict()#[-1]*len(face_cellid)
#        
#        for i in range(len(face_cellid)):
#            if face_name[i] == 1:
#                w_ghost[i] = 1;#w[face_cellid[i][0]]
#            else :
#                w_ghost[i] = w[face_cellid[i][0]]
#        return w_ghost
    
    
    def halo_value(w, Faces, haloFaces, face_cellid, halos, neigh, gtol, ltog):  
        from numpy import zeros
        w_halos = zeros(len(halos), 'd')#mmap.mmap(-1, len(halos))# array('B', [0]) * len(halos)
        #haloexts = zeros(len(halos), 'B')
        
        for i in range(len(halos)):
             w_halos[i] = w[gtol[halos[i]]]#, halos[i]]#(i+1)*(rank+1)#
        
    
        scount = [0  for i in range(size)]
        sdepl  = [0  for i in range(size)]
        rcount = [0  for i in range(size)]
        rdepl  = [0  for i in range(size)]
            
    
        for i in range(len(neigh[0])):
            scount[neigh[0][i]] = neigh[1][i]
        
        for i in range(size):
            if (i>0):
                sdepl[i] = sdepl[i-1] + scount[i-1]
        #obj	=	comm.alltoall(rank)	
        #comm.Alltoall(scount, rcount)
        
        rcount = comm.alltoall(scount)
        
        for i in range(size):
            if (i>0):
                rdepl[i] = rdepl[i-1] + rcount[i-1]
        
        taille = 0
        for i in range(size):
            taille += rcount[i]
                
        w_halor = zeros(taille, 'd')#mmap.mmap(-1, taille)#array('B', [0]) * taille  
        #haloextr = zeros(taille, 'B')
        
        s_msg = [w_halos, (scount, sdepl), MPI.DOUBLE_PRECISION]
        r_msg = [w_halor, (rcount, rdepl), MPI.DOUBLE_PRECISION]
      
        #s_msghalo = [haloexts, (scount, sdepl), MPI.BYTE]
        #r_msghalo = [haloextr, (rcount, rdepl), MPI.BYTE]    
        
        comm.Alltoallv(s_msg, r_msg)
        #comm.Alltoallv(s_msghalo, r_msghalo)
            
        w_halo = OrderedDict()#[-1]*len(Faces)
        for i in range(len(Faces)):
            if face_cellid[i][1] == -10:
                #print(haloFaces.get(tuple([ltog[Faces[i][0]], ltog[Faces[i][1]]])))
                w_halo[i] = w_halor[int(haloFaces.get(tuple([ltog[Faces[i][0]], ltog[Faces[i][1]]]))/3)]
               # w_halo[i] = w_halor[int(haloFaces.index([ltog[Faces[i][0]], ltog[Faces[i][1]]])/3)]
        return w_halo#w_halor, haloextr
    
    
    w_halo = halo_value(w, Faces, HaloFaces, Face_cellid, HaloInt, Neigh, GlobalIndexCell, GlobalIndexNode)
    #w_ghost = ghost_value(w, Face_cellid, Face_name)
    

#    def compute_flux(wl, wr, n):
#        c = 0
#        u = [0]*2
#        u[0] = 1
#        u[1] = 0.
#        q = np.dot(u,n) 
#        if (q >= 0):
#            c = wl
#        else:
#            c = wr   
#        flux = q*c
#        return flux
#    
#    def VecteurNormal(a,b,bary):
#
#        n = [None]*2
#        s = [None]*2
#        m = [None]*2
#        normal = [None]*2
#         
#       
#        n[0] = a[1] - b[1]
#        n[1] = b[0] - a[0];
#        
#        m[0] = 0.5 * (a[0] + b[0]);
#        m[1] = 0.5 * (a[1] + b[1]);
#        
#        s[0] = bary[0] - m[0] ;
#        s[1] = bary[1] - m[1] ;
#        
#        longueur = sqrt((a[0]-b[0])**2 + (a[1]-b[1])**2)
#    
#        n[0] = n[0]
#        n[1] = n[1]
#        
#        if ( (s[0] * n[0] + s[1] * n[1])> 0):
#            normal[0] = -1*n[0];
#            normal[1] = -1*n[1];
#        else:
#            normal[0] = n[0];
#            normal[1] = n[1];
#            
#        normal[0] = normal[0]/longueur
#        normal[1] = normal[1]/longueur 
#        
#        
#        return normal, longueur
#    
#    def barycentre(cells, nodes):
#        bary = [(0, 0)]*len(cells)
#          
#        for i in range(len(cells)):
#            s1 = cells[i][0]
#            s2 = cells[i][1]
#            s3 = cells[i][2]
#            
#            bary[i] = (1./3 * (nodes[s1][0] + nodes[s2][0] + nodes[s3][0]), 1./3 * 
#                (nodes[s1][1] + nodes[s2][1] + nodes[s3][1]))
#        return bary
#    
#    bary = barycentre(Cells, Nodes)
#    
#    
#    rezidus  = [0]*len(Cells)
#    for i in range(len(Faces)):
#        n, longueur = VecteurNormal(Nodes[Faces[i][0]], Nodes[Faces[i][1]], bary[Face_cellid[i][0]])
#        wl = w[Face_cellid[i][0]]
#        if (Face_cellid[i][1] != -1 and Face_cellid[i][1] != -10) :
#            wr = w[Face_cellid[i][1]]
#            flx = compute_flux(wl, wr, n)
#            rezidus[Face_cellid[i][0]] -= flx*longueur 
#            rezidus[Face_cellid[i][1]] += flx*longueur
#      
#        elif (Face_cellid[i][1] == -10):
#            wr = w_halo[i]
#            flx = compute_flux(wl, wr, n)
#            rezidus[Face_cellid[i][0]] -= flx*longueur 
#        elif(Face_cellid[i][1] == -1):
#            wr = w_ghost[i]
#            flx = compute_flux(wl, wr, n)
#            rezidus[Face_cellid[i][0]] -= flx*longueur
    rezidus = 0
    return rezidus

def save_paraview_results(w, n, dt, rank, Cells, Nodes):
    #paraview format
    elements = {"triangle": np.array(Cells)}
    points = [[-1, -1, -1] for i in range(len(Nodes))]
    for i in range(len(Nodes)):
        for j in range(3):
            points[i][j]  = Nodes[i][j]
        
    data = {"wcell" : np.array(w)}  
    data = {"wcell": data}
    maxw = np.zeros(1)
    maxw = max(w)
    integral_sum = np.zeros(1)
    
    if(n%10 == 0):
        comm.Reduce(maxw, integral_sum, MPI.MAX, 0)
        if rank == 0:
            print("saving paraview results, iteration number ", n)
            print("max w =", integral_sum, "dt =", dt)
            
        meshio.write_points_cells("results/visu"+str(rank)+"-"+str(n)+".vtu", points, elements, cell_data=data)

def showMeshPlot(nodes, elements):

    y = []#nodes[:][0]
    z = []#nodes[:][1]
    
    for i in range(len(nodes)):
        y.append((nodes[i][0]))
        z.append((nodes[i][1]))
        
    def quatplot(y,z, quatrangles, ax=None, **kwargs):

        if not ax: ax=plt.gca()
        yz = np.c_[y,z]
        verts= yz[quatrangles]
        pc = matplotlib.collections.PolyCollection(verts, **kwargs)
        ax.add_collection(pc)
        ax.autoscale()

    plt.figure()
    plt.gca().set_aspect('equal')
   
    quatplot(y,z, np.asarray(elements), ax=None, color="crimson", facecolor="None")
    if nodes:            
        plt.plot(y,z,  ls="", color="crimson")

    plt.title('This is the plot for: quad')
    plt.xlabel('Y Axis')
    plt.ylabel('Z Axis')


    plt.show()
    
#Cells, Nodes, HaloInt, HaloExt, GlobalIndexCell, GlobalIndexNode, Neigh = local_mesh()
#Faces, Face_cellid, Cell_faceid , Face_name, HaloFaces = create_Faces_information(Cells, \
#Nodes, HaloExt, GlobalIndexNode)

