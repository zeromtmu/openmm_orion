#############################################################################
# Copyright (C) 2018 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np
import openeye.oechem as oechem
import mdtraj as md

def ExtractProtLigActsiteResNums( mol, fromLigCutoff=5.0):
    '''Extracts the protein and ligand from a single OEMol containing a protein-ligand
    complex plus other components. A list of protein residues within a cutoff distance
    from the ligand is also returned.
    Inputs:
        mol: The OEMol containing protein, ligand, and all other components
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
    Returns:
        protein: An OEMol containing only the protein
        ligand: An OEMol containing only the ligand
        actSiteResNums: A list of integers, one per residue number for protein residues
            with any atom within fromLigCutoff distance of the ligand.'''
    # perceive residue hierarchy of total system
    if not oechem.OEHasResidues(mol):
        oechem.OEPerceiveResidues(mol, oechem.OEPreserveResInfo_All)
    # split the total system into components
    ligand = oechem.OEMol()
    protein = oechem.OEMol()
    water = oechem.OEMol()
    other = oechem.OEMol()
    #sopts = oechem.OESplitMolComplexOptions('MOL')
    sopts = oechem.OESplitMolComplexOptions()
    oechem.OESplitMolComplex(ligand, protein, water, other, mol, sopts)
    # use residue-based distance cutoff (in angstroms) from ligand to define an active site residue
    nn = oechem.OENearestNbrs(protein, fromLigCutoff)
    actSiteResNums = set()
    for nbrs in nn.GetNbrs(ligand):
        residue = oechem.OEAtomGetResidue(nbrs.GetBgn())
        actSiteResNums.add(residue.GetResidueNumber())
    return protein, ligand, actSiteResNums


def ExtractAlignedProtLigTraj_hdf5( mol, traj_hdf5Filename, fromLigCutoff=5.0, skip=0):
    '''Extracts the aligned protein trajectory and aligned ligand trajectory from
    a MD trajectory of a larger system that includes other components (eg water).
    The passed in OEMol must have the topology that matches the trajectory, and its xyz
    coordinates are the reference for the alignment. The alignment is done on the
    alpha carbons (atom name CA) of the those residues specified by atom number
    in the passed in list actSiteResNums. Once the alignment is done, the protein
    and ligand trajectories are each placed into a separate OEMol, one conformer
    per trajectory frame.
    Inputs:
        mol: An OEMol giving the topology for the trajectory and the reference xyz
            coordinates for the alignment.
        trajDCDFilename: The filename of the hdf5-format MD trajectory.
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
        skip: number of frames to skip at the beginning of the trajectory.
    Outputs:
        protTraj: A multiconformer OEMol for the protein, one conformer per frame.
        ligTraj: A multiconformer OEMol for the ligand, one conformer per frame.'''
    # extract protein and ligand molecules from the larger multicomponent system
    # and identify residue numbers for residues within fromLigCutoff of the ligand.
    protein, ligand, actSiteResNums = ExtractProtLigActsiteResNums( mol, fromLigCutoff)
    # get the topology from 1st frame of the traj file
    topologyTraj = md.load_hdf5(traj_hdf5Filename, frame=1)
    # Make a list of the atom indices of the carbon-alphas of the active site residues;
    # assume residue numbering matches the mol
    actSiteCA = [atom.index for atom in topologyTraj.topology.atoms
                    if ((atom.residue.resSeq in actSiteResNums) and (atom.name == 'CA'))]
    # extract protein and ligand subset topology
    #   Note: the ligand must have residue name 'MOL' or 'LIG' (bad, should change)
    protligIdx = topologyTraj.topology.select('protein or resname == MOL or resname == LIG')
    protlig = topologyTraj.atom_slice( protligIdx)
    # Read the protein-ligand subsystem of the trajectory file
    trj_initial = md.load_hdf5(traj_hdf5Filename, atom_indices=protligIdx)
    if skip>0 and len(trj_initial)>skip:
        trj = trj_initial[skip:]
    else:
        trj = trj_initial
    # Fit the protein-ligand trajectory to the active site carbon-alphas of the reference
    protligAtoms = [ atom for atom in protlig.topology.atoms]
    inplace = True
    trjImaged = trj.image_molecules(inplace, [protligAtoms])
    trjImaged.superpose( protlig,0,actSiteCA)
    #Generate a multiconformer representation of the ligand trajectory
    ligIdx = protlig.topology.select('resname == MOL or resname == LIG')
    ligTraj = oechem.OEMol(ligand)
    ligTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in ligIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = ligTraj.NewConf( confxyz)
    # Generate a multiconformer representation of the protein trajectory
    protIdx = protlig.topology.select('protein')
    protTraj = oechem.OEMol(protein)
    protTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in protIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = protTraj.NewConf( confxyz)
    return protTraj, ligTraj


def ExtractAlignedProtLigTraj_dcd( mol, trajDCDFilename, fromLigCutoff=5.0, skip=0):
    '''Extracts the aligned protein trajectory and aligned ligand trajectory from
    a MD trajectory of a larger system that includes other components (eg water).
    The passed in OEMol must have the topology that matches the trajectory, and its xyz
    coordinates are the reference for the alignment. The alignment is done on the
    alpha carbons (atom name CA) of the those residues specified by atom number
    in the passed in list actSiteResNums. Once the alignment is done, the protein
    and ligand trajectories are each placed into a separate OEMol, one conformer
    per trajectory frame.
    Inputs:
        mol: An OEMol giving the topology for the trajectory and the reference xyz
            coordinates for the alignment.
        trajDCDFilename: The filename of the DCD-format MD trajectory.
        fromLigCutoff: The cutoff distance in angstroms to include protein residues
            close to the ligand.
        skip: number of frames to skip at the beginning of the trajectory.
    Outputs:
        protTraj: A multiconformer OEMol for the protein, one conformer per frame.
        ligTraj: A multiconformer OEMol for the ligand, one conformer per frame.'''
    # extract protein and ligand molecules from the larger multicomponent system
    # and identify residue numbers for residues within fromLigCutoff of the ligand.
    protein, ligand, actSiteResNums = ExtractProtLigActsiteResNums( mol, fromLigCutoff)
    # write out pdb file of mol so MDTraj can read it in to get the topology
    pdbFilename = 'ExtractAlignedProtLigTraj_forMDTrajTopology.pdb'
    with oechem.oemolostream(pdbFilename) as ofs:
        flavor = ofs.GetFlavor(oechem.OEFormat_PDB) ^ oechem.OEOFlavor_PDB_OrderAtoms
        ofs.SetFlavor(oechem.OEFormat_PDB, flavor)
        oechem.OEWriteConstMolecule( ofs, mol)
        ofs.close()
    # read the pdb file into an MDTraj topology
    refpdb = md.load_pdb(pdbFilename)
    # Make a list of the atom indices of the carbon-alphas of the active site residues
    actSiteCA = [atom.index for atom in refpdb.topology.atoms
                    if ((atom.residue.resSeq in actSiteResNums) and (atom.name == 'CA'))]
    # extract protein and ligand subset topology
    #   Note: the ligand must have residue name 'MOL' or 'LIG' (bad, should change)
    protligIdx = refpdb.topology.select('protein or resname == MOL or resname == LIG')
    protlig = refpdb.atom_slice( protligIdx)
    # Read the protein-ligand subsystem of the trajectory file
    trj_initial = md.load_dcd(trajDCDFilename, refpdb.topology, None, protligIdx)
    if skip>0 and len(trj_initial)>skip:
        trj = trj_initial[skip:]
    else:
        trj = trj_initial
    # Fit the protein-ligand trajectory to the active site carbon-alphas of the reference
    protligAtoms = [ atom for atom in protlig.topology.atoms]
    inplace = True
    trjImaged = trj.image_molecules(inplace, [protligAtoms])
    trjImaged.superpose( protlig,0,actSiteCA)
    #Generate a multiconformer representation of the ligand trajectory
    ligIdx = protlig.topology.select('resname == MOL or resname == LIG')
    ligTraj = oechem.OEMol(ligand)
    ligTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in ligIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = ligTraj.NewConf( confxyz)
    # Generate a multiconformer representation of the protein trajectory
    protIdx = protlig.topology.select('protein')
    protTraj = oechem.OEMol(protein)
    protTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in protIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = protTraj.NewConf( confxyz)
    return protTraj, ligTraj


