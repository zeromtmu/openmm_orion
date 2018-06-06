#############################################################################
# Copyright (C) 2018 OpenEye Scientific Software, Inc.
#############################################################################
import numpy as np
import openeye.oechem as oechem
import hdbscan


def ClusterTrajRmsd2D(rmsd2D, alpha):
    '''
    Perform Sci-Kit Learn clustering algorithm HDBSCAN on the precomputed
    2D array of distance metrics, returning a list of clusters and a
    list of outliers.
    Inputs arguments:
    rmsd2D: the 2D all-by-all RMSD array for all conformers in the trajectory.
    alpha: the tunable parameter for HDBSCAN, set based on the standard
        deviation of the rmsd2D distribution.
    Returns:
    orderedClusLists: a list of clusters, each cluster a list of conformer
        indices, in order of decreasing size of cluster.
    outliers: a list of conformer indices which are the outliers from HDBSCAN.
    '''
    # Cluster using HDBSCAN
    clusterer = hdbscan.HDBSCAN(metric='precomputed', alpha=alpha)
    clusterer.fit(rmsd2D)
    #
    # Analyze the cluster array
    clusArray = clusterer.labels_
    clusNames = set(clusArray)
    nClus = len(clusNames)-(1 if -1 in clusNames else 0)
    outliers = [ i for i,n in enumerate(clusArray) if n==-1 ]
    #
    # return if no clusters (ie all outliers)
    if nClus<1:
        return clusArray, [], outliers
    # Make list of clusters as lists of conformer indices
    clusLists = []
    for clus in clusNames:
        if clus>-1:
            clusLists.append([])
    for i,clus in enumerate(clusArray):
        if clus>-1:
            clusLists[clus].append(i)
    # sort list from biggest to smallest
    orderedClusLists = clusLists
    if len(orderedClusLists)>1:
        def clusSize(n):
            return len(n)
        orderedClusLists = sorted( clusLists, key=clusSize, reverse=True)

    return orderedClusLists, outliers


def ClusterLigTrajReport(rmsd2D, rmsdArray, clusLists, outliers, reportPrefix):
    '''
    Write a pylab/matplotlib-based report as part of function ClusterLigTraj for
    clustering ligand trajectory conformations.
    Inputs:
    rmsd2D: the 2D all-by-all RMSD array for all conformers in the trajectory.
    rmsdArray: a 1D flattened version of rmsd2D
    clusLists: a list of clusters, each cluster a list of conformer indices.
    outliers: a list of conformer indices corresponding to outliers in the clustering.
    reportPrefix: the filename prefix used to make the report reportPrefix+'clusReport.png'
    '''
    cluster_method = 'HDBSCAN'
    alpha2Std = min( 2.0*rmsdArray.std(), 1.5)
    #
    # Generate text lines for report
    textLines = []
    textLines.append('Trajectory name: %s' % reportPrefix)
    textLines.append( 'HDBSCAN clustering %d items:' % len(rmsd2D))
    textLines.append( '  - using alpha %.1f (min(2*stdev(RMSD),1.5))' % (alpha2Std))
    textLines.append( '  - produced %d clusters with %d outliers' % (len(clusLists), len(outliers)))
    if len(clusLists)>0:
        textLines.append( '  - largest cluster is of size %d' % len(clusLists[0]) )
    if len(clusLists)>1 and len(clusLists[1])>100:
        textLines.append( '  - the second cluster of size %d is also large' % len(clusLists[1]))
    #
    # generate graphic
    # figure layout as 6 rows by one column on 8.5 by 11
    fig = pylab.figure(figsize=(8.5,11.))
    pylab.rcParams['font.size'] = 14
    gspec = pylab.GridSpec(nrows=6, ncols=1, hspace=0.5)
    #
    # RMSD flat distribution
    pltRMS = fig.add_subplot(gspec[0,:])
    pltRMS.hist( rmsdArray)
    pylab.yticks([])
    #
    # RMSD heatmap
    pltHeat = fig.add_subplot(gspec[1:4,:])
    pltHeat.pcolor( rmsd2D)
    #pltHeat.colorbar()
    #
    # Print cluster info to figure
    perLine = 0.022
    ypos = 0.21
    for line in textLines[::-1]:
        ypos = ypos+perLine
        pylab.figtext( 0.15, ypos, line)
    #
    # Strip plot of clusters
    pltClus = fig.add_subplot(gspec[5:,:])
    msize = 10
    for i, clus in enumerate(clusLists):
        pylab.scatter( clus, [i]*len(clus), color='SteelBlue',s=msize)
    outlierPlotVals = [-0.5]*len(outliers)
    pltClus.scatter(outliers,outlierPlotVals, color='k',marker='*',s=msize)
    pylab.xlim( 0,len(rmsd2D))
    pylab.ylim( -.8,)
    yticks = list(range(0,len(clusLists)))
    if len(clusLists)>6:
        yticks = list(range(0,len(clusLists)+1,2))
    pylab.yticks(yticks)
    #
    #print('saving figure under name', reportPrefix+'clusReport.png')
    pylab.savefig(reportPrefix+'clusReport.png')

    return True


def ClusterLigTraj(ligTraj, reportPrefix=None):
    '''
    Clusters conformers based on as-is RMSD and returns a list of
    clusters and the all-by-all 2D array of RMSDs.
    Clustering is done using Sci-Kit Learn's HDBSCAN. A cluster centroid for
    each cluster is identified as the conformer index under the molecule
    Generic Data tag "ClusterCentroid".
    Input arguments:
    ligTraj: a multiconformer OEMol of the trajectory
    reportPrefix: if a .png report is desired, this filename prefix is used.
    Returns:
    rmsd2D: the 2D all-by-all RMSD array for all conformers in the trajectory.
    clusLists: a list of clusters, each cluster a list of conformer indices.
    '''
    automorph = True
    heavyOnly = True
    overlay = False
    rmsd1D = oechem.OEDoubleArray(ligTraj.GetMaxConfIdx())
    #
    # all by all rmsds; this is the time-consuming step
    rmsd2D = []
    for conf in ligTraj.GetConfs():
        oechem.OERMSD(conf, ligTraj, rmsd1D, automorph, heavyOnly, overlay)
        rmsd2D.append( list(rmsd1D))
    #
    # Flatten the 2D array to get global statistics
    flatRmsd2D = []
    for rmsdVec in rmsd2D:
        flatRmsd2D += rmsdVec
    flatRmsd2Darray = np.array(flatRmsd2D)
    #
    # cluster based on the 2D rmsd array
    # Note use of 2.0*rmsdStdev as a parameter in the clustering
    cluster_method = 'HDBSCAN'
    alpha2Std = min( 2.0*flatRmsd2Darray.std(), 1.5)
    clusLists, outliers = ClusterTrajRmsd2D(rmsd2D, alpha2Std)
    #
    # Generate report if reportPrefix provided
    if reportPrefix is not None:
        # print('printing report with reportPrefix', reportPrefix)
        ClusterLigTrajReport(rmsd2D, flatRmsd2Darray, clusLists, outliers, reportPrefix)

    return rmsd2D, clusLists


def FindClusterCentroidFromRMSDArray( rmsd2D, cluster):
    '''
    Find the cluster centroid based on the 2D RMSD array. The method used is
    to examine the all-by-all RMSDs of the conformers in the cluster,
    summing the RMSDs of each conformer with all the others. The conformer
    with the minimum summed RMSDs is chosen to be the centroid.
    Input arguments:
    rmsd2D: the 2D all-by-all RMSD array for all conformers in the trajectory.
    cluster: a list of indices in rmsd2D corresponding to the cluster.
    Returns:
    minRmsdIdx: the cluster centroid's conformer index in the cluster'''
    # Find Idx of cluster centroid from rmsd2D array
    #     method: min of summed rmsds for each conf
    if len(cluster)<1:
        return None
    elif len(cluster)<2:
        return 0
    minRmsdSum = 1.0e10
    minRmsdIdx = -1
    for idx1 in cluster:
        rmsdSum = 0
        for idx2 in cluster:
            rmsdSum += rmsd2D[idx1][idx2]
        if rmsdSum<minRmsdSum:
            minRmsdIdx = idx1
            minRmsdSum = rmsdSum
    return minRmsdIdx


def OEMolFromCluster( trajMol, cluster, centroidIdx=None):
    '''Generate a multiconformer OEMol from a cluster. If centroidIdx is
    specified and that idx is in the cluster list, the OEMol will have
    the new cluster centroid conformer idx in Generic Data Tag "ClusterCentroid".
    Input arguments:
    cluster: a list of conformer indices in multiconformer trajMol.
    trajMol: a multiconformer OEMol of the trajectory; contains the cluster.
    Returns:
    clusMol:  a multiconformer OEMol for the cluster.
    '''
    clusMol = oechem.OEMol(trajMol)
    clusMol.DeleteConfs()
    for i, conf in enumerate(trajMol.GetConfs()):
        if i in cluster:
            if conf.GetIdx()==centroidIdx:
                conf.SetData('ClusterCentroid', 0)
            clusMol.NewConf( conf)
    # Now set the ClusterCentroid Idx (if there is one) on the molecule
    for conf in clusMol.GetConfs():
        if conf.HasData('ClusterCentroid'):
            clusMol.SetData('ClusterCentroid', conf.GetIdx())
    return clusMol


def ProtLigClustersFromTraj( trajLig, trajProt, filePrefix=None):
    '''
    Clusters ligand conformers based on as-is RMSD and returns two lists
    of OEMols: one list of ligand OEMols, one per large cluster, and one list
    of protein OEMols, clustered corresponding to the ligand clusters
    (so each ordered ligand conformer has a matching ordered protein conformer).
    Clustering is done using Sci-Kit Learn's HDBSCAN. A cluster centroid for
    each cluster is identified as the conformer index under the molecule
    Generic Data tag "ClusterCentroid".
    Input arguments:
    trajLig: a multiconformer OEMol of the ligand trajectory
    trajProt: the multiconformer OEMol of the protein trajectory
    filePrefix: if a .png report is desired, this filename prefix is used.
    Returns a pair of lists of OEMols:
    clusterLigs: a list of multiconformer ligand OEMols, one per large cluster.
    clusterProts: a list of multiconformer protein OEMols, matching the ligand clusters.
    '''
    # cluster (and write report .png if given a filePrefix)
    if filePrefix is not None:
        rmsd2d, clusters = ClusterLigTraj(trajLig, filePrefix)
    else:
        rmsd2d, clusters = ClusterLigTraj(trajLig)
    #
    # Create OEMols from large clusters (>10% of trajectory)
    # Tag the cluster centroid
    clusterLigs = []
    clusterProts = []
    for clus in clusters:
        if len(clus)/len(rmsd2d)<0.1:
            #print('cluster of size', len(clus), 'is too small to bother with')
            continue
        centroidIdx = FindClusterCentroidFromRMSDArray( rmsd2d, clus)
        if centroidIdx in clus:
            ligClus = OEMolFromCluster( trajLig, clus, centroidIdx)
            protClus = OEMolFromCluster( trajProt, clus, centroidIdx)
        else:
            ligClus = OEMolFromCluster( trajLig, clus)
            protClus = OEMolFromCluster( trajProt, clus)
        clusterLigs.append( ligClus)
        clusterProts.append( protClus)
    #
    return clusterLigs, clusterProts


