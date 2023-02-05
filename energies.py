import pickle as pkl
import pandas as pd
import os, json, shutil, sys
import numpy as np

def get_code(string):
    return string.split("/")[-1].split("electric")[0][:-1]

def find_geometry(code, mult=True, level="PBE0"):
    found = False
    for i in os.listdir(f"data/crest_dimers/multiple_conformers_{level}/finalized/"):
        if i == str(code) + "_" + str(level) + ".xyz":
            found = True
            break
    full_path = f"data/crest_dimers/multiple_conformers_{level}/finalized/" + str(code) + "_" + str(level) + ".xyz"
    return found, full_path

def write_input(functional, basis_set, dispersion, cores=32, xyz_file="input.xyz"):

    header = f"! SP {functional} {basis_set} {dispersion} RIJCOSX\n"
    processors = f"%PAL NPROCS {cores} END\n"
    memory = f"%maxcore 2000\n"
    geometry = f"*xyzfile 0 1 {xyz_file}\n"
    with open("orca_opt.inp", "w") as input_file:
        input_file.write(header)
        input_file.write(processors)
        input_file.write(memory)
        input_file.write(geometry)

def gather_geometries(samples):

    master = pd.read_pickle(samples)
    comp_files = master["comp_file"].tolist()
    geometries = []
    for list_of_files in comp_files:
        conf_geometries = []
        first_file = list_of_files[0]
        if "mult" in first_file:
            for f in list_of_files:
                code = get_code(f)
                found, path_to_XYZ = find_geometry(code, mult=True)
                if found:
                    conf_geometries.append(path_to_XYZ)
                else:
                    print(f"Could not find XYZ file for {f}")
        if conf_geometries != []:
            geometries.append(conf_geometries)
    
    return geometries


def dft_sp(geometries, levels={}, dispersions=[], orca_path="0rca", benchmark=True, random=5, run_all=False):

    ROOT = os.getcwd()
    os.system("mkdir -p sp_dir/")
    functionals = levels["functionals"]
    basis_sets = levels["basis_sets"]
    dispersions = dispersions
    random_indices = np.random.choice(len(geometries), size=random, replace=False)
    locations = {}
    
    for n, conf in enumerate(geometries):
        if n % 25 == 0 or n == len(geometries):
            print(f"Working on {n}/{len(geometries)}...", flush=True)
        for conf_geom in conf:
            location = []
            code = conf_geom.split("/")[-1].split(".")[0][:-5]
            if not os.path.exists(f"sp_dir/{code}/"):
                for idx, func in enumerate(functionals):
                    func_name = func.split("(")[0]
                    code = conf_geom.split("/")[-1].split(".")[0][:-5]
                    os.system(f"mkdir -p sp_dir/{code}/{func_name}")
                    shutil.copy(conf_geom, f"sp_dir/{code}/{func_name}")
                    write_input(func, basis_sets[idx], dispersions[idx], cores=32, xyz_file=code+"_PBE0.xyz")
                    shutil.move("orca_opt.inp", f"sp_dir/{code}/{func_name}")
            if benchmark and n in random_indices:
                sub_dirs = []
                print(f"\n----- RUNNING BENCHMARK ON: {code} -----\n", flush=True)
                os.chdir(f"sp_dir/{code}/")
                for i in os.listdir("."):
                    sub_dirs.append(i)
                for d in sub_dirs:
                    print(f"\t {d}", flush=True)
                    os.chdir(d)
                    os.system(f"{orca_path} orca_opt.inp > orca_opt.out")
                    print(f"\t DONE", flush=True)
                    os.chdir(ROOT)
                    os.chdir(f"sp_dir/{code}/")

            os.chdir(ROOT)
            if run_all:
                sub_dirs = []
                print(f"\n----- RUNNING {functionals[0]} ON: {code} -----\n", flush=True)
                if not os.path.exists(f"sp_dir/{code}/{functionals[0]}/IN_PROGRESS") and not os.path.exists(f"sp_dir/{code}/{functionals[0]}/DONE"):
                    os.chdir(f"sp_dir/{code}/{functionals[0]}/")
                    with open(f"IN_PROGRESS", "w") as progress:
                        progress.write("IN PROGRESS")
                    os.system(f"{orca_path} orca_opt.inp > orca_opt.out")
                    with open(f"DONE", "w") as complete:
                        complete.write("DONE")
                    os.remove("IN_PROGRESS")
                    print(f"\t DONE", flush=True)
                else:
                    print(f"\t PREVIOUSLY COMPLETED or IN PROGRESS", flush=True)
                    
            os.chdir(ROOT)

            if os.path.exists(f"sp_dir/{code}/{functionals[0]}/orca_opt.out"):
                for idx, func in enumerate(functionals):
                    func_name = func.split("(")[0]
                    location.append(f"sp_dir/{code}/{func_name}/")
                locations[code] = location
                
    return locations

def get_energies(locations):

    E = {}
    for code, locs in locations.items():
        energies = [parse_sp(loc) for loc in locs]
#        energies = [e - energies[-1] for e in energies]
        E[code] = energies
        
    return E

def parse_sp(loc):


    energy = 0.0
    with open(loc+"orca_opt.out", "r") as orca_out:
        for line in orca_out:
            if "FINAL SINGLE POINT ENERGY" in line:
                words = line.split()
                for word in words:
                    try:
                        energy = float(word)
                    except:
                        pass
    return energy

def get_weights(E):

    codes = ["_".join(code.split("_")[:-1]) for code in E.keys()]
    ensemble_energies = {}
    for code in codes:
        energies = []
        for key in E.keys():
            if code in key:
                energies.append(E[key][0])
        if len(energies) == 13:
            print(code)
        ensemble_energies[code] = [(e - min(energies)) for e in energies]
        
    return ensemble_energies
    
geometries = gather_geometries("complexes_1168.pkl")
#levels = {"functionals": ["B3LYP", "PBE0", "M062X", "DLPNO-CCSD(T)"],
#          "basis_sets": ["DEF2-TZVPD def2/J", "DEF2-TZVPD def2/J", "DEF2-TZVPD", "cc-pVTZ cc-pVTZ/C"]}
levels = {"functionals": ["PBE0"],
          "basis_sets": ["DEF2-TZVPD def2/J"]}
#dispersions = ["D4", "D4", "D3zero", ""]
dispersions = ["D4"]
orca_path = "0rca"

locations = dft_sp(geometries,
                   levels=levels,
                   dispersions=dispersions,
                   orca_path=orca_path,
                   benchmark=False,
                   run_all=False)

energies = get_energies(locations)
ensemble_energies = get_weights(energies)

#df = pd.DataFrame.from_dict(energies, orient="index", columns=["B3LYP", "PBE0", "M062X", "DLPNO-CCSD"])
#print(df.min(axis=0))

