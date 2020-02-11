#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 23:10:58 2020

@author: kissami
"""
from mpi4py import MPI
from collections import OrderedDict 
import numpy as np
import meshio

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


__all__ = ["generate_mesh"]


class Cells:
    nodeid = OrderedDict()
    faceid = OrderedDict()
    center = OrderedDict()
    volume = OrderedDict()
    globalIndex = OrderedDict()
    

class Nodes:
    vertex = OrderedDict()
    bound = OrderedDict()
    #pas encore créé
    cellid = OrderedDict()
    globalIndex = OrderedDict()

class Faces:
    nodeid = OrderedDict()
    cellid = OrderedDict()
    bound = OrderedDict()
    normal = OrderedDict()

class Halo:
    halosInt = OrderedDict()
    halosExt = OrderedDict()
    neigh = OrderedDict()
    faces = OrderedDict()

def VecteurNormal(a,b,bary):

    n = [None]*2
    s = [None]*2
    m = [None]*2
    normal = [None]*2
     
   
    n[0] = a[1] - b[1]
    n[1] = b[0] - a[0];
    
    m[0] = 0.5 * (a[0] + b[0]);
    m[1] = 0.5 * (a[1] + b[1]);
    
    s[0] = bary[0] - m[0] ;
    s[1] = bary[1] - m[1] ;
    
    if ( (s[0] * n[0] + s[1] * n[1])> 0):
        normal[0] = -1*n[0];
        normal[1] = -1*n[1];
    else:
        normal[0] = n[0];
        normal[1] = n[1];
        
    normal[0] = normal[0]#/longueur
    normal[1] = normal[1]#/longueur 
    
    
    return normal#, longueur
 
def create_local_mesh(file):
    
    #Lecture des cellules à partir du fichier mesh..txt
    k = 0
    #file = open("mesh"+str(rank)+".txt","r")
    for line in file:
        #read elements
        if line == "Elements\n":
            continue   
        if (line == "EndElements\n"):
            continue
        if line == "Nodes\n":
            break
        Cells.nodeid[k] = [int(x) for x in line.split()] 
        k +=1
    #Lecture des coordonnées des noeuds à partir du fichier mesh..txt
    k = 0
    for line in file:
        #read Nodes
        if line == "Nodes\n":
            continue   
        if (line == "EndNodes\n"):
            continue
        if line == "HalosInt\n":
            break
        Nodes.vertex[k] = [float(x) for x in line.split()]
        Nodes.bound[k] = int(Nodes.vertex[k][3])
        k +=1
        
    
    #calcul du barycentre      
    for i in range(len(Cells.nodeid)):
        s1 = Cells.nodeid[i][0]
        s2 = Cells.nodeid[i][1]
        s3 = Cells.nodeid[i][2]
        
        x1 = Nodes.vertex[s1][0]
        y1 = Nodes.vertex[s1][1]
        x2 = Nodes.vertex[s2][0]
        y2 = Nodes.vertex[s2][1]
        x3 = Nodes.vertex[s3][0]
        y3 = Nodes.vertex[s3][1]
        
        
        Cells.center[i] = (1./3 * (x1 + x2 + x3), 1./3*(y1 + y2 + y3))
        Cells.volume[i] = (1./2) * abs((x1-x2)*(y1-y3)-(x1-x3)*(y1-y2))

    #Création des faces  
    Cellule = Cells.nodeid
    cellF = []
    faces = []
    k = 0
    for i in range(len(Cellule)):
        faces.append([Cellule[i][0], Cellule[i][1]])
        faces.append([Cellule[i][1], Cellule[i][2]])
        faces.append([Cellule[i][0], Cellule[i][2]])
        cellF.append([faces[k], faces[k+1], faces[k+2]])
        k = k+3
    
    faces =   set( tuple(x) for x in faces)
    faces = list(faces)
    
    facesdict = OrderedDict()
    for i in range(len(faces)):
        facesdict[faces[i]] = i
        Faces.nodeid[i] = faces[i]
        Faces.cellid[i] = (-1,-1)
 
            
     
    #Création des 3 faces de chaque cellule            
    for i in range(len(cellF)):
        Cells.faceid[i]  = [facesdict.get(tuple(cellF[i][0])), facesdict.get(tuple(cellF[i][1])), 
                            facesdict.get(tuple(cellF[i][2]))]
        
    #Faces.cellid = [(-1, -1)]*len(faces)
    for i in range(len(Cells.faceid)):
        for j in range(3):
            if Faces.cellid[Cells.faceid[i][j]] == (-1, -1):
                Faces.cellid[Cells.faceid[i][j]]= (i, -1)
            if Faces.cellid[Cells.faceid[i][j]][0] != i:
                Faces.cellid[Cells.faceid[i][j]]= (Faces.cellid[Cells.faceid[i][j]][0], i) 
                
    #Faces aux bords (1,2,3,4), Faces à l'interieur 0    A VOIR !!!!!   
    for i in range(len(faces)):
        Faces.bound[i] = 0
        if (Faces.cellid[i][1] == -1 and Faces.cellid[i][1] != -10):
            if (Nodes.bound[faces[i][0]] == Nodes.bound[faces[i][1]] ):
                Faces.bound[i] = Nodes.bound[faces[i][0]]
            #if ((Nodes.bound[faces[i][0]] == 1 and Nodes.bound[faces[i][1]] !=0) or 
            #    (Nodes.bound[faces[i][0]] != 0 and Nodes.bound[faces[i][1]] ==1)):
            #    Faces.bound[i] = 1
#            if ((Nodes.bound[faces[i][0]] == 2 and Nodes.bound[faces[i][1]] !=0) or 
#                (Nodes.bound[faces[i][0]] != 0 and Nodes.bound[faces[i][1]] ==2)):
#                Faces.bound[i] = 2
        Faces.normal[i] = VecteurNormal(Nodes.vertex[faces[i][0]], Nodes.vertex[faces[i][1]], 
                    Cells.center[Faces.cellid[i][0]])
                
    

    return Cells, Nodes, Faces

def create_halo_structure(file):
     
    #txt_file = open(filename)
    k = 0
    for line in file:
        if "EndHalosInt" in line:
            break
        Halo.halosInt[k] = int(line)#.append(int(line))#for x in line.split()])
        k += 1
    k = 0
    for line in file:
        # process the line
        if "HalosExt" in line :
            continue
        if "GlobalCellToLocal" in line:
            break
        Halo.halosExt[k] = [int(x) for x in line.split()]
        k +=1
    
    cmpt = 0
    for line in file:
    #read Global cell to local
        if line == "GlobalCellToLocal\n":
            continue   
        if (line == "EndGlobalCellToLocal\n"):
            break
        Cells.globalIndex[int(line)] = cmpt
        cmpt +=1
    
    cmpt = 0
    for line in file:
        #read Local Node To Global
        if line == "LocalNodeToGlobal\n":
            continue   
        if (line == "EndLocalNodeToGlobal\n"):
            break
        Nodes.globalIndex[cmpt] = int(line)
        cmpt +=1
        #Nodes.GlobalIndex.append(int(line))# for x in line.split()])  
    k = 0
    for line in file:
        #read LocalToGlobal
        if line == "Neigh\n":
            continue   
        if (line == "EndNeigh\n"):
            break
        Halo.neigh[k] = [int(x) for x in line.split()]#int(line))
        k +=1
       
    #halofaces = []
    if size > 1:
        k = 1
        for i in range(len(Halo.halosExt)):
            Halo.faces[tuple([Halo.halosExt[i][0], Halo.halosExt[i][1]])] = k
            Halo.faces[tuple([Halo.halosExt[i][1], Halo.halosExt[i][2]])] = k + 1
            Halo.faces[tuple([Halo.halosExt[i][0], Halo.halosExt[i][2]])] = k + 2
            k = k+3
    
        for i in range(len(Faces.nodeid)):
            
            if  Halo.faces.get(tuple([Nodes.globalIndex[Faces.nodeid[i][0]],
                                      Nodes.globalIndex[Faces.nodeid[i][1]]])):
               # print(tuple([Nodes.globalIndex[Faces.nodeid[i][0]],
                #                      Nodes.globalIndex[Faces.nodeid[i][1]]]))
                Faces.cellid[i] = (Faces.cellid[i][0], -10)
                Faces.bound[i] = 10
                #print(Faces.bound[i], i, rank)

    return Halo


def ghost_value(w, Faces):
    w_ghost = OrderedDict()#[-1]*len(face_cellid)
    
    for i in range(len(Faces.cellid)):
        if Faces.bound[i] == 2:
            #print(-1*w[Faces.cellid[i][0]])
            w_ghost[i] = w[Faces.cellid[i][0]]
        else :
            w_ghost[i] = w[Faces.cellid[i][0]]
    return w_ghost

def halo_value(w, Cells, Faces, Halo):  
   
    w_halo = OrderedDict()
    if size > 1 :
        from numpy import zeros
        w_halosend = zeros(len(Halo.halosInt), 'd')
        #haloexts = zeros(len(halos), 'B')
       
        
        for i in range(len(Halo.halosInt)):
             w_halosend[i] = w[Cells.globalIndex[Halo.halosInt[i]]]#, halos[i]]#(i+1)*(rank+1)#
        
    
        scount = [0  for i in range(size)]
        sdepl  = [0  for i in range(size)]
        rcount = [0  for i in range(size)]
        rdepl  = [0  for i in range(size)]
            
    
        for i in range(len(Halo.neigh[0])):
            scount[Halo.neigh[0][i]] = Halo.neigh[1][i]
        
        for i in range(size):
            if (i>0):
                sdepl[i] = sdepl[i-1] + scount[i-1]
        
        rcount = comm.alltoall(scount)
        
        for i in range(size):
            if (i>0):
                rdepl[i] = rdepl[i-1] + rcount[i-1]
        
        taille = 0
        for i in range(size):
            taille += rcount[i]
                
        w_halorecv = zeros(taille, 'd')#mmap.mmap(-1, taille)#array('B', [0]) * taille  
       
        
        s_msg = [w_halosend, (scount, sdepl), MPI.DOUBLE_PRECISION]
        r_msg = [w_halorecv, (rcount, rdepl), MPI.DOUBLE_PRECISION]
      
        
        comm.Alltoallv(s_msg, r_msg)            
      
        for i in range(len(Faces.cellid)):
            if Faces.cellid[i][1] == -10:      
                w_halo[i] = w_halorecv[int((-1+Halo.faces.get(
                        tuple([Nodes.globalIndex[Faces.nodeid[i][0]], 
                               Nodes.globalIndex[Faces.nodeid[i][1]]])))/3)]

    return w_halo




def generate_mesh():
    
   
    filename = 'mesh'+str(rank)+'.txt'
    txt_file = open(filename)
    
    
    
    cells, nodes, faces = create_local_mesh(txt_file)
    halos = create_halo_structure(txt_file)
    
    txt_file.close()
    
    return cells, nodes, faces, halos
   

def compute_flux(wl, wr,u, n):
        c = 0
        q = np.dot(u,n) 
        if (q >= 0):
            c = wl
        else:
            c = wr   
        flux = q*c
        return flux

def ExplicitScheme(w, u, w_ghost, w_halo, Faces):
    
    rezidus  = [0]*len(w)
    for i in range(len(Faces.cellid)):
        wl = w[Faces.cellid[i][0]]
        n = Faces.normal[i]
        
        if (Faces.bound[i] == 0):
            wr = w[Faces.cellid[i][1]]
            flx = compute_flux(wl, wr,u[i], n)
            rezidus[Faces.cellid[i][0]] -= flx 
            rezidus[Faces.cellid[i][1]] += flx
      
        elif (Faces.bound[i] == 10):
            wr = w_halo[i]
            flx = compute_flux(wl, wr, u[i], n)
            rezidus[Faces.cellid[i][0]] -= flx 
        
        else:
            wr = w_ghost[i]
            flx = compute_flux(wl, wr,u[i],  n)
            rezidus[Faces.cellid[i][0]] -= flx
        
  
    return rezidus    


def save_paraview_results(w, n, dt, rank, cells, nodes):
    elements = {"triangle": np.array(list(cells))}
    points = [[-1, -1, -1] for i in range(len(nodes))]
    for i in range(len(nodes)):
        for j in range(3):
            points[i][j]  = nodes[i][j]
        
    data = {"wcell" : np.array(w)}  
    data = {"wcell": data}
    maxw = np.zeros(1)
    maxw = max(w)
    integral_sum = np.zeros(1)
    
    if(n%50 == 0):
        comm.Reduce(maxw, integral_sum, MPI.MAX, 0)
        if rank == 0:
            print("saving paraview results, iteration number ", n)
            print("max w =", integral_sum, "dt =", dt)
            
        meshio.write_points_cells("results/visu"+str(rank)+"-"+str(n)+".vtu", 
                                  points, elements, cell_data=data, file_format="vtu")
        
        if(rank == 0 ):
            with open("results/visu"+str(n)+".pvtu", "a") as text_file:
                text_file.write("<?xml version=\"1.0\"?>\n")
                text_file.write("<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\" byte_order=\"LittleEndian\">\n")
                text_file.write("<PUnstructuredGrid GhostLevel=\"0\">\n") 
                text_file.write("<PPoints>\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"binary\"/>\n")
                text_file.write("</PPoints>\n") 
                text_file.write("<PCells>\n") 
                text_file.write("<PDataArray type=\"Int64\" Name=\"connectivity\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Int64\" Name=\"offsets\" format=\"binary\"/>\n")
                text_file.write("<PDataArray type=\"Int64\" Name=\"types\" format=\"binary\"/>\n")
                text_file.write("</PCells>\n") 
                text_file.write("<PCellData Scalars=\"wcell\">\n")
                text_file.write("<PDataArray type=\"Float64\" Name=\"wcell\" format=\"binary\"/>\n")
                text_file.write("</PCellData>\n")
                for i in range(size):
                	name1="visu";
                	bu1 = [10];
                	bu1 = str(i)
                	name1+=bu1;
                	name1+="-"+str(n);
                	name1+=".vtu";
                	text_file.write("<Piece Source=\""+str(name1)+"\"/>\n")                    
                text_file.write("</PUnstructuredGrid>\n") 
                text_file.write("</VTKFile>")
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
