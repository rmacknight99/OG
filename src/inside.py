from __future__ import annotations
import pandas as pd
import argparse
import funsies as f
from funsies.types import Encoding


df = pd.read_csv("../data/raw_data/smiles.csv")
smiles = df["SMILES"].tolist()
id = df["ID"].tolist()
output = df["NAME"].tolist()

def not_empty(inp: bytes) -> bytes:
    """Raise an error if a file is empty."""
    if inp == b"":
        raise RuntimeError("input was empty")
    else:
        return inp

orca_templates = [
    """
! Opt PBE0 def2-SVP def2/J D4 RIJCOSX Pal4
%maxcore 1200
*xyz {{charge}} {{spin}}
{{structure}}
*
""",
    """
! Freq wB97X def2-TZVP def2/J RIJCOSX Pal4
%maxcore 1200
*xyzfile {{charge}} {{spin}} orca_input0.xyz
""",
    """
! wB97X def2-TZVP def2/J RIJCOSX Pal4
%maxcore 1200
%tddft
maxdim 5
nroots 20
end
*xyzfile {{charge}} {{spin}} orca_input0.xyz
""",
    """
! wB97X def2-TZVP def2/J RIJCOSX Pal4
%maxcore 1200
%tddft
maxdim 5
nroots 20
end
*xyzfile {{charge}} {{spin_plus_two}} orca_input0.xyz
""",
]

for p, e in enumerate(smiles):
    SMILES = smiles[p]
    ID = id[p]
    OUTPUT = output[p]
    SPIN = 1
    CHARGE = 0
    SPIN_PLUS_TWO = 3
    OUTPUTFILE0 = OUTPUT + "_" + str(ID) + "_" + "opt_S0.out"
    OUTPUTFILE1 = OUTPUT + "_" + str(ID) + "_" + "freq_S0.out"
    OUTPUTFILE2 = OUTPUT + "_" + str(ID) + "_" + "tddft_S0.out"
    OUTPUTFILE3 = OUTPUT + "_" + str(ID) + "_" + "tddft_S1.out"
    XYZFILE0 = OUTPUT + "_" + str(ID) + "_" + "opt_S0.xyz"

    with f.Fun():

        xyzs = f.shell(
        f"obabel -:'{SMILES}' --addhs --gen3d --ff uff --minimize -O mol1.xyz",
        out=["mol1.xyz"],
        )

        opt1 = f.shell(
        "xtb mol1.xyz --opt -P 4",
        inp={"mol1.xyz": xyzs.out["mol1.xyz"]},
        out=["xtbopt.xyz"],
        )		

        crest_conf = f.shell(
        "crest xtbopt.xyz --gfnff//gfn2 --nci",
        inp={"xtbopt.xyz": opt1.out["xtbopt.xyz"]},
        out=["crest_best.xyz"],
        )

        best = f.morph(
        not_empty,
        crest_conf.out["crest_best.xyz"],
        )

        orca_input0 = f.template(
        orca_templates[0],
        dict(
            charge=CHARGE,
            spin=SPIN,
            structure=f.utils.truncate(best, top=2),
        ),
        )

        orca_input1 = f.template(
        orca_templates[1],
        dict(
            charge=CHARGE,
            spin=SPIN
        ),
        )

        orca_input2 = f.template(
        orca_templates[2],
        dict(
            charge=CHARGE,
            spin=SPIN
        ),
        )

        orca_input3 = f.template(
        orca_templates[3],
        dict(
            charge=CHARGE,
            spin_plus_two=SPIN_PLUS_TWO
        ),
        )

        orca0 = f.shell(
            "/opt/orca/orca orca_input0.inp > orca_output0.out",
            inp={"orca_input0.inp": orca_input0},
            out=["orca_output0.out", "orca_input0.xyz"])
    
        orca1 = f.shell(
                "/opt/orca/orca orca_input1.inp > orca_output1.out",
                inp={"orca_input1.inp": orca_input1, "orca_input0.xyz": orca0.out["orca_input0.xyz"]},
                out=["orca_output1.out"]) #save hess; save xyz

        orca2 = f.shell(
                "/opt/orca/orca orca_input2.inp > orca_output2.out",
                inp={"orca_input2.inp": orca_input2, "orca_input0.xyz": orca0.out["orca_input0.xyz"]},
                out=["orca_output2.out"]) #save xyz
        
        orca3 = f.shell(
                "/opt/orca/orca orca_input3.inp > orca_output3.out",
                inp={"orca_input3.inp": orca_input3, "orca_input0.xyz": orca0.out["orca_input0.xyz"]},
                out=["orca_output3.out"]) #save xyz

        f.execute(orca1)
        f.execute(orca2)
        f.execute(orca3)
        f.wait_for(orca1)
        f.wait_for(orca2)
        f.wait_for(orca3)
        f.takeout(orca0.out["orca_output0.out"], OUTPUTFILE0)
        f.takeout(orca0.out["orca_input0.xyz"], XYZFILE0)
        f.takeout(orca1.out["orca_output1.out"], OUTPUTFILE1)
        #f.takeout(orca1.out["orca_output1.out"], "orca_output1.out")
        #f.takeout(orca1.out["orca_input1.hess"], HESSIANFILE1)
        f.takeout(orca2.out["orca_output2.out"], OUTPUTFILE2)
        f.takeout(orca3.out["orca_output3.out"], OUTPUTFILE3)

