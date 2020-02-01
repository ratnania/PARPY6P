Il faut d'abord executer le fichier seqtest.py pour créer les fichiers locaux.
La fonction MeshPart prend deux arguments le nombre de partitions souhaitées
et le fichier du maillage (.msh)
Exemple :
MeshPart(2, "mesh.msh") : crée deux fichiers mesh0.txt et mesh1.txt

Par la suite il faut faire :
mpirun -n 2 python partest.py
