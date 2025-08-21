# ruff: noqa: T201

import logging
from pathlib import Path
import os

import numpy as np
import stk
import stko
from metallicious import supramolecular_structure

def build_complex(path_to_complex, opt_dist):

    bb1 = stk.BuildingBlock(
            smiles='[Pd+2]',
            functional_groups=(
                stk.SingleAtom(stk.Pd(0, charge=2))
                for i in range(4)
            ),
            position_matrix=[[0, 0, 0]],
        )

    bb2 = stk.BuildingBlock(
                smiles='C1=CC=NC=C1',
                functional_groups=[
                    stk.SmartsFunctionalGroupFactory(
                        smarts='[#6]~[#7X2]~[#6]',
                        bonders=(1, ),
                        deleters=(),
                    ),
                ],
            )

    complex = stk.ConstructedMolecule(
                topology_graph=stk.metal_complex.SquarePlanar(
                    metals=bb1,
                    ligands=bb2,
                    optimizer=stk.Collapser(step_size=0.1, distance_threshold=opt_dist, scale_steps=True)
                ),
            )
    #print(complex._bonds)
    complex.write(path_to_complex)
    return complex

def main() -> None:
        """Run metallicious example 1"""
        cage = supramolecular_structure('example_1/ru_pd.pdb', metal_charges={'Ru': 2, 'Pd':2 }, topol='example_1/ru_pd.top', LJ_type='uff')
        cage.parametrize(out_coord='out.pdb', out_topol='out.top')
        
        print("example 1 ran successfully")

        """Run metallicious example 2"""
        cage = supramolecular_structure('example_2/ru_pd.xyz', metal_charges={'Ru':2, 'Pd':2 }, LJ_type='uff')
        cage.prepare_initial_topology()
        cage.parametrize(out_coord='out.pdb', out_topol='out.top')

        print("example 2 ran successfully")

        """Run the stk example."""
        script_dir = Path(__file__).resolve().parent  #We want the output to be relative to the location of the script and not the working directory

        output =script_dir/"simple_met_output"
        output.mkdir(exist_ok=True, parents=True)

        stk_mol = build_complex(output / 'PdComplex.xyz',2.5)

        os.chdir(Path(__file__).resolve().parent) #to prevent nested init_topol folders
        complex = supramolecular_structure(str(output / 'PdComplex.xyz'), metal_charges={'Pd': 2}, LJ_type='uff', truncation_scheme = "angle")#, covalent_cutoff = 2, rmsd_cutoff = 10)
        complex.stk_molecule = stk_mol  
        complex.parametrize(out_coord='out.pdb', out_topol='out.top', prepare_initial_topology=True)

        print("stk example ran successfully")

if __name__ == "__main__":
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
    )
    main()
