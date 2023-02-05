import funsies as f
import argparse, os, sys
from funsies.types import Encoding

lst = sys.argv[0].split("/")[:-2]
root = '/'.join(lst)

parser = argparse.ArgumentParser(description="Run CREST conformer generation for EDA complexes")
parser.add_argument("COMPLEX", help="", type=str)
parser.add_argument("-c", "--charge", help="charge", type=int, default=0)
parser.add_argument("-s", "--spin", help="spin", type=int, default=1)

args = parser.parse_args()

COMPLEX_XYZ = args.COMPLEX
COMPLEX_ROOT = COMPLEX_XYZ.split("/")[-1][:-4]
CHARGE = args.charge
SPIN = args.spin
SPIN_PLUS_TWO = SPIN + 2

# Files 
xTB_output = COMPLEX_ROOT + '_' + 'xTB' + '_' + 'opt' + '.xyz'

with open(COMPLEX_XYZ, "r") as xyz_file:
    next(xyz_file)
    next(xyz_file)
    for line in xyz_file:
        linelist = line.split()
        atoms.append(linelist[0])

indices = []

for i,e in enumerate(atoms):
    if e == 'O':
        if [atoms[i-2], atoms[i-1], e] == ['O', 'S', 'O']:
            indices.append(i+1)
        elif e == 'C':
            if [atoms[i-2], atoms[i-1], e, atoms[i+1]] == ['F','F', 'C', 'F']:
                indices.append(i+1)
        elif e == 'S':
            if [atoms[i-1], e, atoms[i+1]] == ['F', 'S', 'C']:
                indices.append(i+1)

if indices == []:
    indices = [1, 2, 3]
    print(indices)

# load the xyz file into a variable to give to funsies                                                                                                
with open(COMPLEX_XYZ, "rb") as tmp:
    COMPLEX_XYZ = tmp.read()

with f.Fun():
    COMPLEX_XYZ = COMPLEX_XYZ
    # encode charge and spin                                                                                                                                 
    chrg = f.morph(lambda x: str(x).encode(), CHARGE, out=Encoding.blob)
    s0 = f.morph(lambda x: str(x).encode(), SPIN, out=Encoding.blob)
    s1 = f.morph(lambda x: str(x).encode(), SPIN_PLUS_TWO, out=Encoding.blob)

    # optimize initial cctk complex with xTB                                                                                                          
    xTB_opt = f.shell(
        "xtb complex.xyz --opt",
        inp={"complex.xyz": COMPLEX_XYZ, ".CHRG": chrg, ".s0": s0, ".s1": s1},
        out=["xtbopt.xyz"],
        env={"OMP_NUM_THREADS": str(os.cpu_count())},
    )

    # generate conformers using crest
    crest = f.shell(
        "crest input.xyz -gfn2//gfnff -niceprint -nci",
        inp={"input.xyz": xTB_opt.out["xtbopt.xyz"], ".CHRG": chrg, ".s0": s0, ".s1": s1},
        out=["crest_best.xyz", "crest_conformers.xyz", "crest.energies", "crest_rotamers.xyz", "cre_members"],
        env={"OMP_NUM_THREADS": str(os.cpu_count())},
    )

    # constrain atoms by index
    con_crest = f.shell(
        f"crest input.xyz --constrain {indices[0]},{indices[1]},{indices[2]}",
        inp={"input.xyz": xTB_opt.out["xtbopt.xyz"],
             ".CHRG": chrg, ".s0": s0, ".s1": s1},
        out=[".xcontrol.sample", "coord.ref"],
        env={"OMP_NUM_THREADS": NTASKS},
    )

    # cinp argument
    crest_topo = f.shell(
        f"crest input.xyz -gfn2//gfnff -niceprint -nci --cinp constraints.inp",
        inp={"input.xyz": xTB_opt.out["xtbopt.xyz"],
             "constraints.inp": con_crest.out[".xcontrol.sample"],
             ".CHRG": chrg, ".s0": s0, ".s1": s1,
             "coord.ref": con_crest.out["coord.ref"]},
        out=["crest_best.xyz", "crest_conformers.xyz", "crest.energies", "crest_rotamers.xyz", "cre_members"],
        env={"OMP_NUM_THREADS": NTASKS},
    )

    if COMPLEX_ROOT.endswith('_4')
        f.execute(crest_topo)
        f.wait_for(crest_topo)
    elif COMPLEX_ROOT.startswith('4_'):
        f.execute(crest_topo)
        f.wait_for(crest_topo)
    else:
        f.execute(crest)
        f.wait_for(crest)

    directory = f"{root}/crest_conformers/{COMPLEX_ROOT}/"
    try:
        os.mkdir(directory)
    except:
        pass
    for i in crest.out:
        f.takeout(crest.out[i], directory+i)
    
    
    
