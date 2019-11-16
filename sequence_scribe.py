#!/usr/bin/env python3

from __future__ import print_function
import csv
from sys import argv
from rdkit import Chem
from rdkit.Chem import Draw



def main():
    # read in sequence of masses and monomer pool
    script, sequence_raw, monomer_list = argv

    parent_masses = []
    sequence = []

    monomer_pool = {}

    # open csvs and populate parent masses list
    with open(sequence_raw, 'rt') as csvfile:
        mass_data = csv.reader(csvfile)
        next(mass_data)
        for entry in mass_data:
            entry = float(entry[1])
            parent_masses.append(entry)

    # open csvs and populate monomer pool list
    with open(monomer_list, 'rt') as csvfile:
        monomers_data = csv.reader(csvfile)
        next(monomers_data)
        for entry in monomers_data:
            monomer_pool[float(entry[0])] = entry[1]

    # subtract the mass lost to the parent mass to get the individual monomer mass
    current_mass = 0
    for mass in parent_masses:
        monomer_mass = abs(abs(current_mass) - mass)
        for monomer in monomer_pool:
            if ((monomer_mass > monomer-1) and (monomer_mass < monomer+1)):
                sequence.append(monomer_pool[monomer])
            else:
                continue
        current_mass = abs(mass)

    # get the final mass and determine the final monomer that is fluorophored labelled
    final_fluorophore = parent_masses[len(parent_masses)-1]
    for monomer in monomer_pool:
        if ((final_fluorophore > monomer-2) and (final_fluorophore < monomer+2)):
            sequence.append(monomer_pool[monomer])
        else:
            continue

    molecule_list = []


    # convert smiles code to molecules
    for x in sequence:
        mol = Chem.MolFromSmiles(x)
        molecule_list.append(mol)

    # write molecules into a png
    molecules_per_row = len(molecule_list)/2
    img = Draw.MolsToGridImage(molecule_list, molsPerRow=4, subImgSize=(200, 200))
    img.save(sequence_raw[:len(sequence_raw) - 4] +".png")


main ()


