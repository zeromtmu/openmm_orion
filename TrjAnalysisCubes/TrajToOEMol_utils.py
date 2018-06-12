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


def AnalyseProteinLigandTrajectoryOEMols( ligTraj, protTraj):
    '''Analyse an MD trajectory which has been mapped as conformers on to a
    ligand OEMol and a protein OEMol.
      Input:
        ligTraj: the multiconformer trajectory OEMol for the ligand.
        protTraj: the multiconformer trajectory OEMol for the protein.
      Output:
        ligMedian: a single conformer OEMol for the median ligand.
        protForMedianLig: a single conformer OEMol for the lig-median protein.
        ligAverage: a single conformer OEMol for the average ligand.
        protForMedianLig: a single conformer OEMol for the average protein.
    The analyses are B-Factor-like rms fluctuations around the mean and
    dihedral fluctuations on ligand heavy-atom rotors. The results are returned
    on a single configuration OEMol for the ligand and protein. The single
    configuration corresponds to the median configuration for the ligand;
    the protein configuration comes from the same frame in the trajectory as
    the median ligand configuration. The B-Factors are set in the
    OEResidue data on the protein and ligand atoms. The ligand heavy
    atom rotor fluctuations are mapped on the ligand rotors as Generic Data
    with tags "TorIdxs", "TorAngles", "TorWrapped", "TorMean", and "TorStdev":
      TorIdxs: integer quadruples referring to the atoms that make the torsion
               around that rotor. The references are to atom Generic data "cIdx"
      TorAngles: a one-per-conformer list of the torsion angles in radians.
      TorWrapped: like TorAngles but wrapped around the circle so the torsion
               angle values are continuous around the mean value (see TorMean).
      TorMean: the mean value of the wrapped torsion angles.
      TorStdev: the standard deviation of the wrapped torsion angles. When this
               value is high it means the torsion fluctuates more.'''
    #Calculate B-Factors for the ligand
    ligAvg = ConfAverage( ligTraj)
    ligMedian = BFactorsFromTrajMolAndAverageMol( ligTraj, ligAvg)

    #Add rotor fluctuation statistics to rotors (as Generic Data on bonds)
    ligTrajWithRotorStats = RotorConfFluctuationStats( ligTraj)

    # Copy rotor fluctuation statistics to ligMedian
    # Note: this was not put in a helper function because ligTrajWithRotorStats and
    # ligMedian must have identical graphs and atom/bond ordering, and we are not
    # checking for it here.
    for fromAtom, toAtom in zip(ligTrajWithRotorStats.GetAtoms(), ligMedian.GetAtoms()):
        if fromAtom.HasData("cIdx"):
            toAtom.SetData("cIdx", fromAtom.GetData("cIdx"))
    tags = ["TorIdxs", "TorAngles", "TorWrapped", "TorMean", "TorStdev"]
    for fromBond, toBond in zip(ligTrajWithRotorStats.GetBonds(), ligMedian.GetBonds()):
        for tag in tags:
            if fromBond.HasData(tag):
                toBond.SetData(tag, fromBond.GetData(tag))

    # Calculate B-Factors for the protein
    protAvg = ConfAverage( protTraj)
    protMedian = BFactorsFromTrajMolAndAverageMol( protTraj, protAvg)

    #Extract protein snapshot corresponding to ligand median; copy B-Factors from protMedian
    if ligMedian.HasData("Traj_MedianConfIdx"):
        ligMedianConfIdx = ligMedian.GetData("Traj_MedianConfIdx")
        protForMedianLig = oechem.OEMol( protTraj.GetConf(oechem.OEHasConfIdx( ligMedianConfIdx)))
        if not oechem.OEHasResidues(protForMedianLig):
            oechem.OEPerceiveResidues(protForMedianLig, OEPreserveResInfo_All)
        for atomFrom, atomTo in zip( protMedian.GetAtoms(), protForMedianLig.GetAtoms() ):
            resFrom = oechem.OEAtomGetResidue(atomFrom)
            resTo = oechem.OEAtomGetResidue(atomTo)
            resTo.SetBFactor( resFrom.GetBFactor() )
    else:
        protForMedianLig = protMedian

    return ligMedian, protForMedianLig, ligAvg, protAvg


def RotorConfFluctuationStats( mol):
    '''For each heavy-atom rotor, calculate the statistics of a unique torsion
    associated with the rotor. The list of per-conformer torsions, the mean,
    and standard deviation are attached to the rotor bond as Generic Data.
    Input: a multiconformer OEMol
    Returns: the same molecule with the following Generic Data attached to
        each rotor, listed here by Generic Data Tag:
        TorIdxs: a list of four atoms tags (atom Generic Data tag"cIdx")
            identifying the atoms of the measured torsion for this rotor.
        TorAngles: a list of the per-conformer torsion angles in radians
        TorWrapped: a periodicity-corrected (wrapped) version of TorAngles
        TorMean: the mean value in radians for this torsion
        TorStdev: the standard deviation in radians for this torsion'''
    # Generate the stats for fluctuations around heavy-atom rotors
    tors = FindUniqueRotorTorsions( mol)
    tlist = torsListByConf( mol, tors)
    twrapped, tmeans, tdevs = torsWrapAroundMean( tlist)
    # attach persisting atom Idxs to atoms for reference by Generic Data TorIdxs
    for atom in mol.GetAtoms():
        atom.SetData("cIdx", atom.GetIdx())
    # attach to each rotor the stats for fluctuations
    for tor, torAngles, torWrapped, torMean in zip(
      tors, tlist, twrapped, tmeans):
        torIdxs = [tor[0].GetIdx(), tor[1].GetIdx(), tor[2].GetIdx(), tor[3].GetIdx()]
        torStdev = np.std(torWrapped)
        #print( torIdxs[0], torIdxs[1], torIdxs[2], torIdxs[3], torStdev)
        for bond in mol.GetBonds():
            atmBgnIdx = bond.GetBgn().GetIdx()
            atmEndIdx = bond.GetEnd().GetIdx()
            if ((atmBgnIdx==torIdxs[1] or atmBgnIdx==torIdxs[2])
                   and (atmEndIdx==torIdxs[1] or atmEndIdx==torIdxs[2])):
                #print( 'Matches bonded atoms', atmBgnIdx, atmEndIdx)
                bond.SetData("TorIdxs", torIdxs)
                bond.SetData("TorAngles", torAngles)
                bond.SetData("TorWrapped", torWrapped)
                bond.SetData("TorMean", torMean)
                bond.SetData("TorStdev", torStdev)
    return mol

def BFactorsFromTrajMolAndAverageMol( trajMol, avgMol, findMedian=True):
    '''Returns an OEMol with isotropic B-factors on each atom under
    Generic Data tag "Traj_BFactor". With findMedian=True (default)
    the OEMol returned is the robust estimate of the closest conformer
    to avgMol; if False avgMol is returned.
    Input:
        trajMol: OEMol with one conformer for each frame of the trajectory
        avgMol: OEMol with the trajectory average xyz coords for trajMol
        findMedian: bool if True (default) return the estimate of the closest
            conformer to avgMol
    Returns: an OEMol with calculated isotropic B-factors on each atom'''
    medianConfIdx = 0
    xyzDevsAbsSumMin = 1.0E+10
    # place avgMol xyz into numpy array
    xyzOEarray = oechem.OEFloatArray(3*avgMol.GetMaxAtomIdx())
    avgMol.GetCoords(xyzOEarray)
    xyzMean = np.array( xyzOEarray)
    # accumulate squared xyz deviations from the mean for each conformer
    xyzDevsSumSq = np.zeros( 3*trajMol.GetMaxAtomIdx() )
    for conf in trajMol.GetConfs():
        conf.GetCoords(xyzOEarray)
        xyz = np.array( xyzOEarray)
        xyzDevs = xyz-xyzMean
        xyzDevsAbsSum = np.sum( np.absolute(xyzDevs))
        if findMedian and xyzDevsAbsSum<xyzDevsAbsSumMin:
            medianConf = conf
            medianConfIdx = conf.GetIdx()
            xyzDevsAbsSumMin = xyzDevsAbsSum
        xyzDevsSq = np.multiply( xyzDevs, xyzDevs)
        xyzDevsSumSq = np.add( xyzDevsSumSq, xyzDevsSq)
    xyzAvgDevsSq = xyzDevsSumSq/float(trajMol.NumConfs())
    atomAvgDevsSq = xyzAvgDevsSq.reshape( trajMol.GetMaxAtomIdx(),3)
    atomAvgFlucsSq = atomAvgDevsSq.sum(axis=1)
    AtomBFactorsIso = (8*np.pi*np.pi/3.0)*atomAvgFlucsSq
    if findMedian:
        medianMol = oechem.OEMol( medianConf)
        medianMol.SetData( 'Traj_MedianConfIdx', medianConf.GetIdx() )
        if not oechem.OEHasResidues(medianMol):
            oechem.OEPerceiveResidues(medianMol, OEPreserveResInfo_All)
        for atom, bfac in zip( medianMol.GetAtoms(), AtomBFactorsIso):
            val = min( bfac, 150)
            res = oechem.OEAtomGetResidue(atom)
            res.SetBFactor( val)
        return medianMol
    else:
        if not oechem.OEHasResidues(avgMol):
            oechem.OEPerceiveResidues(avgMol, OEPreserveResInfo_All)
        for atom, bfac in zip( avgMol.GetAtoms(), AtomBFactorsIso):
            val = min( bfac, 150)
            res = oechem.OEAtomGetResidue(atom)
            res.SetBFactor( val)
    return avgMol

def ConfAverage( mol):
    '''Return a new OEMol having xyz coords being the as-is average of those of
    the conformations in the input multi-conformer OEMol
    input: a multi-conformer OEMol
    output: a single-conformer OEMol containing averaged xyz coords'''
    # iterate over all conformational xyz coords calculating the average in a numpy array
    xyzSum = np.zeros( 3*mol.GetMaxAtomIdx() )
    xyzOEarray = oechem.OEFloatArray(3*mol.GetMaxAtomIdx())
    for conf in mol.GetConfs():
        conf.GetCoords(xyzOEarray)
        xyz = np.array( xyzOEarray)
        xyzSum = np.add( xyzSum, xyz)
    xyzMean = xyzSum/float(mol.NumConfs())
    # copy the means to the OEFloatArray
    for i in list( range(3*mol.GetMaxAtomIdx())):
        xyzOEarray[i] = xyzMean[i]
    # make a new OEMol containing just those confs
    molAvg = oechem.OEMol( mol)
    molAvg.DeleteConfs()
    molAvg.NewConf( xyzOEarray)
    return molAvg

def ConfCOMs( mol):
    '''Returns the numpy array of the per-conformer molecular Center Of Mass (COM).
    input: a multi-conformer OEMol
    returns: numpy array of one float per conformer'''
    com = oechem.OEFloatArray(3)
    trjCOMlist = []
    for conf in mol.GetConfs():
        oechem.OEGetCenterOfMass( conf,com)
        trjCOMlist.append( (com[0],com[1],com[2]) )
    return np.array( trjCOMlist)

def XyzDistFromRefXyz( xyz, RefXyz):
    '''Calculates deviations in xyz coords from reference xyz coords.
    Input:
        xyz: nX3 numpy array of xyz coords
        RefXyz: list or array of length 3 of the reference xyz coords
    Returns: nX3 numpy array of deviations from the reference xyz coords'''
    devs = xyz-RefXyz
    devsSq = np.multiply( devs, devs)
    return np.sqrt( devsSq.sum(axis=1))

def torsListByConf( mol, torsAtoms):
    '''For each torsion torsAtoms (quadruple of OEAtoms) we iterate through every
    conformer building a list of the torsion angle measurements.
    input:
        mol: the multiconformer molecule whose torsions are being measured.
        torsAtoms: the list of OEAtom quadruples representing the torsions to be measured.
    output:
        torsList: the list of lists with each list of floats containing the torsion angle
        measurements over all conformers for one of the torsions in torsAtoms.'''
    torsList = []
    for tor in torsAtoms:
        torList = [ oechem.OEGetTorsion(conf,tor[0], tor[1], tor[2], tor[3]) for conf in mol.GetConfs() ]
        torsList.append( torList)
    return torsList

def minAnglePeriodicDistance( angle, refAngle):
    '''To get the minimum angular distance, this function compares angle to refAngle,
    wrapping the angle if necessary so that its value is within pi radians of the refAngle
    input:
        angle (in radians) to be wrapped if necessary
        refAngle (in radians) for comparison
    ouput: angle (in radians), wrapped if it was necessary'''
    ang = angle
    absDiff = abs( ang-refAngle)
    if absDiff>np.pi:
        if angle<0:
            ang += 2*np.pi
        else:
            ang -= 2*np.pi
    return ang

def torsWrapAroundMean( torsList):
    '''For each list of torsion angles (in radians) this function uses circular statistics
    to find the circular mean, then wraps all the angles around the mean correcting for the
    periodicity so that no angle is more than pi radians away from the mean
    input:
        torsList: a list of lists of angles (in radians)
    output:
        torsMean: a list of means, one for each list of angles
        torsWrapped: a list of lists of angles (in radians), with angles wrapped around the mean'''
    torsMeans = []
    torsWrapped = []
    torsDevs = []
    for i, torList in enumerate(torsList):
        tor = np.array(torList)
        sinTor = np.sin( tor)
        cosTor = np.cos( tor)
        meanSinTor = np.mean( np.sin( tor))
        meanCosTor = np.mean( np.cos( tor))
        # good approximation of the circular mean
        meanTor = np.arctan2(meanSinTor,meanCosTor)
        # wrap angles based on this mean
        wrapped = np.array( [ minAnglePeriodicDistance( startAngle, meanTor) for startAngle in tor] )
        torsWrapped.append( wrapped)
        # recalculate a more accurate mean
        meanTor = np.mean(wrapped)
        torsMeans.append( meanTor )
        # calculate deviations of angles from the mean
        torsDevs.append( wrapped-meanTor )
    return torsWrapped, torsMeans, torsDevs


def HeavyConnections(atom):
    '''Determines the summed atomic numbers of the atom and all its connected atoms
    input: an OEatom
    returns: the summed atomic numbers (an int)'''
    nbrSum = atom.GetAtomicNum()
    for nbr in atom.GetAtoms():
        nbrSum += nbr.GetAtomicNum()
        #print('atom neighbours:', nbr.GetIdx(), nbr.GetAtomicNum(), nbrSum)
    return nbrSum

def FindHeaviestTorsion( torList):
    '''From a torsion list of OEAtom quadruples, this function returns the quadruple
    having the heaviest extended terminal atoms as determined by summing the
    atomic numbers of the terminal atoms and all their connections.
    input: a list of OEAtom quadruples, each quadruple defining a torsion.
    output: a single OEAtom quadruple, defining the heaviest torsion.
    '''
    if len( torList)<1:
        return None
    best = torList[0]
    bestTorScore = HeavyConnections( best[0]) + HeavyConnections( best[3])
    for tor in torList:
        torScore = HeavyConnections( tor[0]) + HeavyConnections( tor[3])
        if torScore>bestTorScore:
            best = tor
            bestTorScore = torScore
    return best

def FindUniqueRotorTorsions( mol):
    '''For each rotor in the molecule, this function finds the heaviest torsion
    returning a list of one torsion for each rotor. Heaviest is defined as having
    the highest summed atomic numbers for the terminal atoms of the torsion,
    including all atoms connected to the terminal atoms.
    input: an OEMol
    output: a list of OEAtom quadruples, one for each rotor in the molecule.
    '''
    uniqueRotors = []
    firstpass = True
    ctr1 = None
    ctr2 = None
    grp = []
    for tor in oechem.OEGetTorsions(mol):
        if firstpass or tor.b!=ctr1 or tor.c!=ctr2:
            # this is a torsion with a new central bond
            firstpass = False
            if len( grp)>0:
                heavyTor = FindHeaviestTorsion( grp)
                uniqueRotors.append( heavyTor)
            grp = [[tor.a, tor.b, tor.c, tor.d]]
            #print('New central bond:', tor.a.GetIdx(), tor.b.GetIdx(), tor.c.GetIdx(), tor.d.GetIdx())
            ctr1 = tor.b
            ctr2 = tor.c
        else:
            # add to the list
            grp.append( [tor.a, tor.b, tor.c, tor.d])
            #print('  Adding new tor to group:', tor.a.GetIdx(), tor.b.GetIdx(), tor.c.GetIdx(), tor.d.GetIdx())
    # process the last group of torsions upon exiting the loop
    if len( grp)>0:
        heavyTor = FindHeaviestTorsion( grp)
        uniqueRotors.append( heavyTor)
    return uniqueRotors


