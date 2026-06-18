from functools import partial
from time import time
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdForceFieldHelpers, rdDistGeom
from rdkit.Chem.rdchem import RWMol, BondType, BondStereo
from rdkit.Chem.rdmolops import AddHs
import random

def smiles_to_xyz_MMFF(smi):
    xyz = None
    success = None

    # special treat for H2
    if smi == "[H][H]":
        success = True
        xyz = (("H","0.000000","0.000000","0.3"),("H","0.000000","0.000000","-0.3"))
        return success, {"xyz": xyz}
    else:                    
        try:
            mol = AddHs(Chem.MolFromSmiles(smi))  
            # ETDG
            p = rdDistGeom.ETKDG()
            p.randomSeed = 42
            p.numThreads = 1
            p.useBasicKnowledge = False
            rdDistGeom.EmbedMultipleConfs(mol,numConfs=40,params=p)
            molff = None        
            
            if rdForceFieldHelpers.MMFFHasAllMoleculeParams(mol):
                cand = [1 for _ in range(mol.GetNumConformers())]
                ffprops = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
                molff = partial(rdForceFieldHelpers.MMFFGetMoleculeForceField, pyMMFFMolProperties=ffprops)
            elif rdForceFieldHelpers.UFFHasAllMoleculeParams(mol):
                cand = [1 for _ in range(mol.GetNumConformers())]
                ffprops = None
                molff = partial(rdForceFieldHelpers.UFFGetMoleculeForceField)
            else:
                cand = []

        except:
            success = False 
            return success, {"xyz": xyz}

        # MMFF loop
        for i in range(4):
            forcetol = 0.1*(0.1**i)
            energytol = 0.001*(0.1**i)
            en_cand_all = []
            # MMFF
            for id_conf in range(len(cand)):
                if cand[id_conf]:
                    ff = molff(mol, confId=id_conf)
                    if ff.Minimize(maxIts=1000, forceTol=forcetol, energyTol=energytol):
                        cand[id_conf] = 0
                        continue
                    en = ff.CalcEnergy()
                    en_cand_all.append((en, id_conf))
            en_cand_all.sort()

            if sum(cand) == 0:
                success = False 
                return success, {"xyz": xyz}
            
            # if final loop determine conformer
            if i == 3:
                break

            # Delete similar or high energy comformers
            for j, state in enumerate(en_cand_all):
                en, id = state
                if j != 0 and en - en_cand_all[0][0] > 15.0:
                    cand[id] = 0
                    continue
                for k in range(j):
                    if getConformerMaxAtomDistance(mol, en_cand_all[k][1],id) < 0.5:
                        cand[id] = 0
                        break

            if sum(cand) == 1:
                break
        
        success = True
        xyz = []
        pos = mol.GetConformers()[en_cand_all[0][1]].GetPositions()
        random.seed(314)
        for k,atom in enumerate(mol.GetAtoms()):
            xyz.append((atom.GetSymbol(),f'{pos[k][0]+random.gauss(0,0.01): .6f}',f'{pos[k][1]+random.gauss(0,0.01): .6f}',f'{pos[k][2]+random.gauss(0,0.01): .6f}'))
        xyz = tuple(xyz)
        return success, {"xyz": xyz}
    
def getConformerMaxAtomDistance(mol, confId1, confId2, atomIds=None, prealigned=False):
    if not prealigned:
        if atomIds:
            AllChem.AlignMolConformers(mol, confIds=[confId1, confId2], atomIds=atomIds)
        else:
            AllChem.AlignMolConformers(mol, confIds=[confId1, confId2])
    conf1 = mol.GetConformer(id=confId1)
    conf2 = mol.GetConformer(id=confId2)
    dmax = 0
    for i in range(mol.GetNumAtoms()):
        dmax = max(dmax,conf1.GetAtomPosition(i).Distance(conf2.GetAtomPosition(i)))
    return dmax

def output_xyz(xyz, filename):
    with open(filename,"w",encoding="utf-8") as f:
        f.write(f"{len(xyz)}\n\n")
        for atom in xyz:
            f.write(f"{atom[0]} {atom[1]} {atom[2]} {atom[3]}\n")





    
