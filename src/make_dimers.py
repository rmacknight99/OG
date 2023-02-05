import cctk
import numpy as np
import copy
import itertools
import os, sys

def center_molecule(molecule):
        """
        Moves a ``Molecule`` object's centroid to the origin
        """
        atoms = np.arange(0, molecule.num_atoms())
        molecule.translate_molecule(-molecule.geometry[atoms].mean(axis=0))
        return molecule

def spherical_random(radius=1):
    """
    Generates a random point on a sphere of radius ``radius``.
    """
    v = np.random.normal(size=3)
    v = v / np.linalg.norm(v)
    return v * radius

def form_complex(xyz_files, radius=6, title=''):

    mol1 = cctk.XYZFile.read_file(f"./data/output_data/{xyz_files[0]}").molecule
    mol2 = cctk.XYZFile.read_file(f"./data/output_data/{xyz_files[1]}").molecule
    mol1 = center_molecule(mol1)
    mol2 = center_molecule(mol2)
    complex_size = len(xyz_files)

#    title = '+'.join(xyz_files)
    num_structures = 10

    for i in range(num_structures):
        trans_v = spherical_random(radius=radius)
        for j in range(complex_size-1):
            x = copy.deepcopy(mol2)
            x.translate_molecule(trans_v)
            x.rotate_molecule(np.array([1, 0, 0]), np.random.random() * 360)
            x.rotate_molecule(np.array([0, 1, 0]), np.random.random() * 360)
            x.rotate_molecule(np.array([0, 0, 1]), np.random.random() * 360)
            atoms = np.hstack((mol1.atomic_numbers.T, x.atomic_numbers.T))
            geoms = np.vstack((mol1.geometry, x.geometry))
            mx = cctk.Molecule(atomic_numbers=atoms, geometry=geoms, charge=0, multiplicity=1)
            cctk.XYZFile.write_molecule_to_file(f"./data/test/{title}.xyz", mx, title=f"{title}")

opt_files = []
for i in os.listdir('./data/output_data/'):
  if 'opt' and 'xyz' in i:
    opt_files.append(i)
print(f"{len(opt_files)} total molecules")
combos = list(itertools.combinations(opt_files, 2))

select_combos = []
try: 
        search = sys.argv[1]
except:
        search = None

if search is not None:
        for i in combos:
                if search in i[0] or search in i[1]:
                        select_combos.append(i)
else:
        select_combos = combos
                        
print(f"{len(select_combos)} potential EDA complexes")
titles = []
for plex in select_combos:
        ID = ''
        for item in plex:
                for index, character in enumerate(item.split('_')):
                        if character == 'opt':
                                ID += item.split('_')[index-1]+'_'
        titles.append(ID[:-1])

for index, files in enumerate(select_combos):
  xyz_files = list(files)
  form_complex(xyz_files,radius=8,title=titles[index])
