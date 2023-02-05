"""Generate a single conformer from a SMILES string."""
from __future__ import annotations

import argparse, os, glob
import funsies as f
from funsies.types import Encoding

def not_empty(inp: bytes) -> bytes:
    """Raise an error if a file is empty."""
    if inp == b"":
        raise RuntimeError("input was empty")
    else:
        return inp


# Parameters
parser = argparse.ArgumentParser(
    description="Compute the spectrum of a molecule using the Vertical Gradient method."
)

parser.add_argument("-C", "--complex", help="complex", type=str)
parser.add_argument("-c", "--charge", help="charge", type=int, default=0)
parser.add_argument("-s", "--spin", help="spin", type=int, default=1)

args = parser.parse_args()

if len(glob.glob("*xyz")) == 1:
    COMPLEX_XYZ = glob.glob('*xyz')[0] # complex file name 
else:
    COMPLEX_XYZ = args.complex # complex file name
    
COMPLEX_ROOT = COMPLEX_XYZ[:-4] # ID

CHARGE = args.charge
SPIN = args.spin
SPIN_PLUS_TWO = SPIN + 2

ID = COMPLEX_ROOT # ID

# output file names for takeout

low_opt_out = ID + "_low_opt_S0.out"
low_opt_xyz= ID + "_low_opt_S0.xyz"
high_opt_out = ID + "_high_opt_S0.out"
high_opt_xyz = ID + "_high_opt_S0.xyz"

freq_out = ID + "_freq_S0.out"
hess = ID + "_freq_hessian_S0.hess"

tddft_s0 = ID + "_tddft_S0.out"
tddft_s1 = ID + "_tddft_S1.out"

# ORCA templates

orca_templates = [
 """
! Opt PBE def2-SVPD def2/J D4 RIJCOSX
%PAL NPROCS 32 END
%maxcore 2000
*xyzfile {{charge}} {{spin}} complex.xyz
 """,
    """
! Opt PBE0 def2-TZVPD def2/J D4 RIJCOSX 
%PAL NPROCS 32 END
%maxcore 2000
*xyzfile {{charge}} {{spin}} orca_low_opt.xyz
""",
    """
! Freq wB97X def2-TZVPD D4 def2/J RIJCOSX 
%PAL NPROCS 32 END
%maxcore 2000
*xyzfile {{charge}} {{spin}} orca_high_opt.xyz
""",
    """
! wB97X def2-TZVPD def2/J D4 RIJCOSX 
%PAL NPROCS 32 END
%maxcore 2000
%tddft
maxdim 5
nroots 20
end
*xyzfile {{charge}} {{spin}} orca_high_opt.xyz
""",
    """
! wB97X def2-TZVPD def2/J D4 RIJCOSX 
%PAL NPROCS 32 END
%maxcore 2000
%tddft
maxdim 5
nroots 20
end
*xyzfile {{charge}} {{spin_plus_two}} orca_high_opt.xyz
""",
]

# load the xyz file into a variable to give to funsies
with open(COMPLEX_XYZ, "rb") as tmp:
    COMPLEX_XYZ = tmp.read()


with f.Fun():
    COMPLEX_XYZ = COMPLEX_XYZ

    # Setup ORCA input
    orca_low_opt = f.template(
            orca_templates[0],
            dict(
            charge=CHARGE,
            spin=SPIN),
            )

    orca_high_opt = f.template(
            orca_templates[1],
            dict(
            charge=CHARGE,
            spin=SPIN),
            )

    orca_freq = f.template(
            orca_templates[2],
            dict(
            charge=CHARGE,
            spin=SPIN),
            )
    
    orca_tddft_s0 = f.template(
            orca_templates[3],
            dict(
            charge=CHARGE,
            spin=SPIN),
            )
    
    orca_tddft_s1 = f.template(
            orca_templates[4],
            dict(
            charge=CHARGE,
            spin_plus_two=SPIN_PLUS_TWO),
            )


    orca0 = f.shell(
            "/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca orca_low_opt.inp > orca_low_opt.out",
            inp={"orca_low_opt.inp": orca_low_opt, "complex.xyz": COMPLEX_XYZ},
            out=["orca_low_opt.out", "orca_low_opt.xyz"]
    )


    orca1 = f.shell(
            "/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca orca_high_opt.inp > orca_high_opt.out",
            inp={"orca_high_opt.inp": orca_high_opt, "orca_low_opt.xyz": orca0.out["orca_low_opt.xyz"]},
            out=["orca_high_opt.out", "orca_high_opt.xyz"]
    )

    orca2 = f.shell(
            "/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca orca_freq.inp > orca_freq.out",
            inp={"orca_freq.inp": orca_freq, "orca_high_opt.xyz": orca1.out["orca_high_opt.xyz"]},
            out=["orca_freq.out", "orca_freq.hess"]
    )

    orca3 = f.shell(
            "/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca orca_tddft_s0.inp > orca_tddft_s0.out",
            inp={"orca_tddft_s0.inp": orca_tddft_s0, "orca_high_opt.xyz": orca1.out["orca_high_opt.xyz"]},
            out=["orca_tddft_s0.out"]
    )

    orca4 = f.shell(
            "/opt/packages/orca/orca_5_0_1_linux_x86-64_shared_openmpi411/orca orca_tddft_s1.inp > orca_tddft_s1.out",
            inp={"orca_tddft_s1.inp": orca_tddft_s1, "orca_high_opt.xyz": orca1.out["orca_high_opt.xyz"]},
            out=["orca_tddft_s1.out"]
    )

    f.execute(orca1)
    f.execute(orca2)
    f.execute(orca3)
    f.execute(orca4)

    f.wait_for(orca1)
    f.wait_for(orca2)
    f.wait_for(orca3)
    f.wait_for(orca4)

    f.takeout(orca0.out["orca_low_opt.out"], low_opt_out)
    f.takeout(orca0.out["orca_low_opt.xyz"], low_opt_xyz)

    f.takeout(orca1.out["orca_high_opt.out"], high_opt_out)
    f.takeout(orca1.out["orca_high_opt.xyz"], high_opt_xyz)

    f.takeout(orca2.out["orca_freq.out"], freq_out)
    f.takeout(orca2.out["orca_freq.hess"], hess)

    f.takeout(orca3.out["orca_tddft_s0.out"], tddft_s0)

    f.takeout(orca4.out["orca_tddft_s1.out"], tddft_s1)
    
