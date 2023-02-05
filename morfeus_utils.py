import os
from morfeus.conformer import ConformerEnsemble


def load_ensemble(code):

    ensemble = []
    for i in sorted(os.listdir("./outs/")):
        if i.startswith(code):
            ensemble.append(i[:-4])

    return sorted(ensemble)

def gen_crest(list_xyz_files):
    # generate crest conformer ensemble (crest_conformers.xyz, cre_members, crest.energies)
    # from a group of xyz files for the same molecule (list_xyz_files)

    # initialize lists
    cre_members = []
    energies = []
    xyz = []
    raw_files = []
    
    cre_members.append(len(list_xyz_files)) # cre_members starts with number of conformers

    col_0 = 1 # always 1 degeneracy ?
    col_1 = 0 # initialize cumulative sum
    col_2 = 0 # initialize cumulative sum
    for xyz_file in list_xyz_files:
        with open(f"./opts/{xyz_file}.xyz", "r") as f: # open xyz files
            lines = f.readlines()
            num_atoms = lines[0] # first line as number of atoms
        with open(f"./outs/{xyz_file}.out") as orca_out: # open orca output
            orca_lines = orca_out.readlines() 
            for orca_line in orca_lines:
                if "FINAL SINGLE POINT" in orca_line:
                    energy = float(orca_line.split()[-1]) # get energy for this conformer

        energies.append(energy) # append energy

        xyz_lines = [] # for coordinates
        for line in lines[2:]:
            xyz_lines.append(line) # append coordinated

        col_0 = 1
        col_1 += 1
        col_2 += 1

        cre_members_string = f"\t{col_0}\t{col_1}\t{col_2}"
        cre_members.append(cre_members_string) # append to cre_members
        
        xyz.append([str(num_atoms)] + [str(energy)+"\n"] + xyz_lines) # make xyz_file lines
        raw_files.append(xyz_file) # append to list of all files

    energies = [(e - min(energies))*627.5 for e in energies] # convert to kcal/mol from Hartree
    zipped = zip(energies, xyz, raw_files) # zipped list for sorting
    zipped = list(zipped)
    res = sorted(zipped, key = lambda x: x[0])
    energies, xyz, raw_files = [list(tup) for tup in zip(*res)] # unpack sorted by energy

    # write the files
    
    for index, xyz_f in enumerate(xyz):
        xyz_f[1] = str(energies[index])+"\n"
    try:
        os.remove(f"crest_conformers.xyz")
        os.remove(f"cre_members")
        os.remove(f"crest.energies")
    except:
        pass

    with open(f"crest_conformers.xyz", "w") as conformers:
        for xyz_f in xyz:
            for line in xyz_f:
                conformers.write(line)

    with open(f"cre_members", "w") as cre_memb:
        for line in cre_members:
            cre_memb.write(str(line)+"\n")

    with open(f"crest.energies", "w") as e:
        for index, energy in enumerate(energies):
            e.write(str(index+1)+"\t"+str(energy)+"\n")

    return energies

def get_boltz_weights(path, prune=True):

    ce = ConformerEnsemble.from_crest(path)
    if prune:
        ce.prune_energy()
    return ce.boltzmann_weights()

def prune(list_of_xyz_files, energies, ID, conf_track_dict, energy_thresholds=[3.0], RMSD_thresholds=[1.5,2.0,2.5]):
    # morfeus conformer pruning
    
    ce = ConformerEnsemble.from_crest(".")

    energy_file_dict = {}
    for i, filename in enumerate(list_of_xyz_files):
         energy_file_dict[round(energies[i], 3)] = filename

    # prune Energy
    energy_prune_files = {key: [] for key in energy_thresholds}
    for energy in energy_thresholds:
        print(f"\tEnergies more than {energy} kcal/mol away from lowest energy conformer will be removed ...")
        ce.prune_energy(threshold=energy)
        relative_energies = ce.get_relative_energies()
        for re in relative_energies:
            energy_prune_files[energy].append(energy_file_dict[round(re,3)])
        conf_track_dict[ID].append(len(energy_prune_files[energy]))

    # prune RMSD
    rmsd_prune_files = {key: [] for key in RMSD_thresholds}
    for rmsd in RMSD_thresholds:
        print(f"\tPruning RMSD, structures with pair-wise RMSD less than {rmsd} A will be considered geometrically degenerate")
        ce.prune_rmsd(thres=rmsd, method="obrms-iter")
        relative_energies = ce.get_relative_energies()
        for re in relative_energies:
            rmsd_prune_files[rmsd].append(energy_file_dict[round(re,3)])
        conf_track_dict[ID].append(len(rmsd_prune_files[rmsd]))

    return conf_track_dict, energy_prune_files, rmsd_prune_files