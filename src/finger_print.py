from rdkit import Chem
from rdkit.Chem import AllChem
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('--mol', type=str,
			help = 'Mol file')

args = parser.parse_args()

mol = Chem.MolFromMolFile(args.mol)
#fp = Chem.RDKFingerprint(mol)

radius=3
nBits=1024

ECFP6 = AllChem.GetMorganFingerprintAsBitVect(mol,radius=radius, nBits=nBits) 
ecfp6_bits = list(ECFP6)
ecfp6_bits = [str(element) for element in ECFP6]
ecfp6_bits = ''.join(ecfp6_bits)


print(args.mol)         
with open('finger_prints.csv','a') as f:
	f.write(f"{args.mol[:-4]},{ecfp6_bits}\n")
