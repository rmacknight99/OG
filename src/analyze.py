import sys, os, shutil
from morfeus.conformer import ConformerEnsemble

# bottom of file contains a for loop running the functions defined in this script on
# the ./opts/ directory, which contains an ensemble of xyz files

def gen_crest(list_xyz_files, ID):
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
        with open(f"./opts/{xyz_file}", "r") as f: # open xyz files
            lines = f.readlines()
            num_atoms = lines[0] # first line as number of atoms
        with open(f"./outs/{xyz_file[:-4]}.out") as orca_out: # open orca output
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

    for xyz_f in xyz:
        with open(f"crest_conformers.xyz", "a") as conformers:
            for line in xyz_f:
                conformers.write(line)

    with open(f"cre_members", "w") as cre_memb:
        for line in cre_members:
            cre_memb.write(str(line)+"\n")

    with open(f"crest.energies", "w") as e:
        for index, energy in enumerate(energies):
            e.write(str(index+1)+"\t"+str(energy)+"\n")

    return energies

def prune(list_of_xyz_files, energies, ID, conf_track_dict):
    # morfeus conformer pruning
    
    ce = ConformerEnsemble.from_crest(".")

    energy_file_dict = {}
    for i, filename in enumerate(list_of_xyz_files):
         energy_file_dict[round(energies[i], 3)] = filename

    # prune energy
    ce.prune_energy()
    relative_energies = ce.get_relative_energies()
    energy_prune_files = []
    for re in relative_energies:
        energy_prune_files.append(energy_file_dict[round(re,3)])
    conf_track_dict[ID].append(len(energy_prune_files))

    # prune rmsd=1.5
    ce.prune_rmsd(thres=1.5, method="obrms-iter")
    relative_energies = ce.get_relative_energies()
    rmsd_15_prune_files = []
    for re in relative_energies:
        rmsd_15_prune_files.append(energy_file_dict[round(re,3)])
    conf_track_dict[ID].append(len(rmsd_15_prune_files))

    # prune rmsd=2.0
    ce.prune_rmsd(thres=2.0, method="obrms-iter")
    relative_energies = ce.get_relative_energies()
    rmsd_20_prune_files = []
    for re in relative_energies:
        rmsd_20_prune_files.append(energy_file_dict[round(re,3)])
    conf_track_dict[ID].append(len(rmsd_20_prune_files))

    # prune rmsd=2.5
    ce.prune_rmsd(thres=2.5, method="obrms-iter")
    relative_energies = ce.get_relative_energies()
    rmsd_25_prune_files = []
    for re in relative_energies:
        rmsd_25_prune_files.append(energy_file_dict[round(re,3)])
    conf_track_dict[ID].append(len(rmsd_25_prune_files))

    return conf_track_dict


conf_track_dict = {} # initialize empty dict for tracking conformers
ID_list = [] # initialize ID list for tracking
all_xyz_files = [] # for tracking
for xyz_file in os.listdir('./opts/'): # get IDs and filenames based on path
    ID = xyz_file[:-17].split("_")[0] + "_" + xyz_file[:-17].split("_")[1] # this will be without confID
    ID_list.append(ID) # essentially a counter for conformers 
    all_xyz_files.append(xyz_file.split("/")[-1]) # filename
ID_set = set(ID_list) # this is the set of all IDs in ./opts/

for ID in ID_set: # iterate through unique Ids in ./opts/

    conf_track_dict[ID] = [ID_list.count(ID)] # update tracking dict with the starting number of conformers
    list_of_xyz_files = [] # get list of conformer xyz files for this ID
    for f in all_xyz_files:
        if ID in f:
            list_of_xyz_files.append(f) # append if ID is in the filename

    energies = gen_crest(list_of_xyz_files, ID) # run gen_crest on list of xyz files
    try:
        conf_track_dict = prune(list_of_xyz_files, energies, ID, conf_track_dict) # get the updated tracking dict with prune
    except:
        pass

print(conf_track_dict) 
