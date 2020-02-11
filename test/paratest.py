from PARPY6P import meshpart
from mpi4py import MPI
import numpy as np
import timeit

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


start = timeit.default_timer()

cells, nodes, faces, halos = meshpart.generate_mesh()

w=wn=np.zeros(len(cells.center))
vit = [0] * len(faces.nodeid)

for i in range(len(vit)):
    vit[i] = [2.,0]

nbelements = len(cells.center)
w=wn=np.zeros(nbelements)

for i in range(nbelements):
    if(cells.center[i][0] > 20 and cells.center[i][0] < 25 ):
        w[i] = 10
        

t = 0
Tfinal = 20
dt = 0

n = 0    
while(t<Tfinal):
    dt = 0.01
    t = t + dt
    n += 1
    w_halo = meshpart.halo_value(w, cells, faces, halos)
    w_ghost = meshpart.ghost_value(w, faces)      
    rezidus = meshpart.ExplicitScheme(w, vit, w_ghost, w_halo, faces)
    
    for i in range(nbelements):
        wn[i]= w[i] + dt * (rezidus[i]/cells.volume[i])
               
    meshpart.save_paraview_results(w, n, dt, rank, cells.nodeid.values(), nodes.vertex)

    
    w = wn
stop = timeit.default_timer()

if rank == 0 : print(stop - start)

