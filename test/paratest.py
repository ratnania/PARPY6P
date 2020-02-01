from PARPY6P import meshpart
import numpy as np
from mpi4py import MPI
import timeit

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()
#meshpart.LocalStructure.
#meshpart.LocalStructure.ExplicitScheme()


#
start = timeit.default_timer()
Cells, Nodes, HaloInt, HaloExt, GlobalIndexCell, GlobalIndexNode, Neigh = meshpart.LocalStructure.local_mesh()
Faces, Face_cellid, Cell_faceid , Face_name, HaloFaces = meshpart.LocalStructure.create_Faces_information(Cells, \
Nodes, HaloExt, GlobalIndexNode)
stop = timeit.default_timer()
print("Creating local structures", stop-start)


t = 0
Tfinal = 8
dt = 0
w=wn=np.zeros(len(Cells))

start = timeit.default_timer()
rezidus = meshpart.LocalStructure.ExplicitScheme(w, Cells, Nodes, Cell_faceid, Faces, Face_cellid, Face_name, HaloFaces, \
                                                      HaloInt, Neigh, GlobalIndexCell, GlobalIndexNode)
stop = timeit.default_timer()
print("Compute rezidus", stop-start)


#w=wn=np.zeros(len(Cells))
#
#n = 0    
#while(t<Tfinal):
#    dt = 0.01
#    t = t + dt
#    n += 1    
#    rezidus = meshpart.LocalStructure.ExplicitScheme(w, Cells, Nodes, Cell_faceid, Faces, Face_cellid, Face_name, HaloFaces, \
#                                                      HaloInt, Neigh, GlobalIndexCell, GlobalIndexNode)
#    
#
#    
#    for i in range(len(Cells)):
#        wn[i]= w[i] + (dt/0.017)*rezidus[i]
#        
#        
#    meshpart.LocalStructure.save_paraview_results(w, n, dt, rank, Cells, Nodes)
#
#    
#    w = wn
#    






