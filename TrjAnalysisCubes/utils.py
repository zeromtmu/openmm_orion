#############################################################################
# Copyright (C) 2018 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np
import openeye.oechem as oechem
import mdtraj as md
from datarecord import (Types, OEField, OERecord)


def GetCardinalOrderOfProteinResNums( mol):
    # make map of protein res nums to the residue cardinal order index
    resmap = {}
    currRes = -10000
    currIdx = -1
    for atom in mol.GetAtoms(oechem.OEIsBackboneAtom()):
        thisRes = oechem.OEAtomGetResidue(atom)
        resnum = thisRes.GetResidueNumber()
        if resnum!=currRes:
            currIdx += 1
            currRes = resnum
            resmap[currRes] = currIdx
    return resmap, currIdx


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
    alpha carbons (atom name CA) of the active site residues within fromLigCutoff
    from the ligand. Once the alignment is done, the protein and ligand trajectories
    are each placed into a separate OEMol, one conformer per trajectory frame.
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
    # get the topology from 1st frame of the traj file
    topologyTraj = md.load_hdf5(traj_hdf5Filename, frame=1)
    # Put the reference mol xyz into the 1-frame topologyTraj to use as a reference in the fit
    molXyz = oechem.OEDoubleArray( 3*mol.GetMaxAtomIdx())
    mol.GetCoords( molXyz)
    molXyzArr = np.array( molXyz)
    molXyzArr.shape = (-1,3)
    # convert from angstroms to nanometers and slice out the protein-ligand complex
    topologyTraj.xyz[0] = molXyzArr/10.0
    # extract protein and ligand molecules from the larger multicomponent system
    # and identify residue numbers for residues within fromLigCutoff of the ligand.
    protein, ligand, actSiteResNums = ExtractProtLigActsiteResNums( mol, fromLigCutoff)
    protResMap, numProtRes = GetCardinalOrderOfProteinResNums( protein)
    actSiteResIdxs = set()
    for resnum in actSiteResNums:
        actSiteResIdxs.add( protResMap[resnum])
    # extract protein atom indices: cannot trust mdtraj protein selection so
    # assume they are contiguous and starting the atom list and just get the same
    # number of atoms as in the OpenEye protein
    protOEIdx = np.array( [ atom.GetIdx() for atom in protein.GetAtoms()] )
    # extract ligand atom indices
    #   Note: the ligand must have residue name 'MOL' or 'LIG' (bad, should change)
    ligIdx = topologyTraj.topology.select('resname == MOL or resname == LIG')
    protligIdx = np.append( protOEIdx, ligIdx)
    #print( 'numAtoms prot, lig, protlig:', len(protOEIdx), len(ligIdx), len(protligIdx))
    #protligIdx = topologyTraj.topology.select('protein or resname == MOL or resname == LIG')
    # Read the protein-ligand subsystem of the trajectory file
    trj_initial = md.load_hdf5(traj_hdf5Filename, atom_indices=protligIdx)
    if skip>0 and len(trj_initial)>skip:
        trj = trj_initial[skip:]
    else:
        trj = trj_initial
    # Image the protein-ligand trajectory so the complex does not jump across box boundaries
    protlig = topologyTraj.atom_slice( protligIdx)
    protligAtoms = [ atom for atom in protlig.topology.atoms]
    inplace = True
    trjImaged = trj.image_molecules(inplace, [protligAtoms])
    # Make a list of the atom indices of the carbon-alphas of the active site residues;
    # assume residue numbering matches the mol
    actSiteCA = [atom.index for atom in topologyTraj.topology.atoms
                    if ((atom.residue.resSeq in actSiteResIdxs) and (atom.name == 'CA'))]
    # Fit the protein-ligand trajectory to the active site carbon-alphas of the reference
    trjImaged.superpose( protlig,0,actSiteCA)
    #Generate a multiconformer representation of the ligand trajectory
    ligTraj = oechem.OEMol(ligand)
    ligTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in ligIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = ligTraj.NewConf( confxyz)
    # Generate a multiconformer representation of the protein trajectory
    strNumProteinAtomsToSelect = 'index '+str(protOEIdx[0])+' to '+str(protOEIdx[-1])
    protIdx = protlig.topology.select( strNumProteinAtomsToSelect)
    protTraj = oechem.OEMol(protein)
    protTraj.DeleteConfs()
    for frame in trjImaged.xyz:
        xyzList = [ 10*frame[idx] for idx in protIdx ]
        confxyz = oechem.OEFloatArray( np.array(xyzList).ravel() )
        conf = protTraj.NewConf( confxyz)
    return protTraj, ligTraj


def RequestOEField( record, field, rType):
    if not record.has_value(OEField(field,rType)):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field))
        raise ValueError('The record does not have field {}'.format( field))
    else:
        return record.get_value(OEField(field,rType))

def RequestOEFieldType( record, field):
    if not record.has_value(field):
        #opt['Logger'].warn('Missing record field {}'.format( field))
        print( 'Missing record field {}'.format( field.get_name() ))
        raise ValueError('The record does not have field {}'.format( field.get_name() ))
    else:
        return record.get_value(field)

def MakeClusterInfoText(dataDict):
    # Generate text string about Clustering information
    #
    text = []
    nFrames = dataDict['nFrames']
    text.append('Cluster method {}\n'.format( dataDict['ClusterMethod']) )
    text.append('- Clustered {} frames\n'.format(nFrames) )
    text.append('- Used alpha={:.2f}\n'.format( dataDict['HDBSCAN_alpha']))
    #
    if dataDict['nClusters']<2:
        text.append('- Produced {} cluster\n'.format( dataDict['nClusters']))
    else:
        text.append('- Produced  {} clusters\n'.format( dataDict['nClusters']))
    nOutliers = dataDict['ClusterVec'].count(-1)
    text.append('    with {:4d} outliers\n\n'.format( nOutliers))
    #
    text.append(' Cluster Size Status\n')
    text.append(' ------- ---- ------\n')
    for i, count in enumerate(dataDict['ClusterCounts']):
        if nFrames/count>10:
            text.append('  {:2d}    {:4d}  minor\n'.format( i, count))
        else:
            text.append('  {:2d}    {:4d}  major\n'.format( i, count))
    #
    return text


