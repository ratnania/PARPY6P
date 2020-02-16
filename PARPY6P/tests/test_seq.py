from PARPY6P import meshpart
import os

# ... get the mesh directory
try:
    mesh_dir = os.environ['MESH_DIR']

except:
    base_dir = os.path.dirname(os.path.realpath(__file__))
    base_dir = os.path.join(base_dir, '..', '..')
    mesh_dir = os.path.join(base_dir, 'mesh')
# ...

def test_mesh_1():
    filename = os.path.join(mesh_dir, "mesh.msh")
    meshpart.MeshPartitioning.MeshPart(2, filename)
