#############
Release Notes
#############


v0.13.6
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to ``OpenEye-orionplatform==0.1.14``

New Cubes
--------------------------------------------------------------------------------
* New :ref:`cube_ConformerSplitterCube` and :ref:`cube_ConformerMergerCube` cubes which preserve fields added during
  calculations on the split conformers.
    * previous cubes for Splitting and Merging conformers which did not preserve fields have been removed.
* New :ref:`cube_MolTitleToField` cube to create a field with the title of the molecule in a chosen MolField.
* New :ref:`cube_GameplanCube` added that performs Gameplan calculations.
* New :ref:`cube_DGConformerCube` added that generates molecule conformers using distance geometry.
* New :ref:`cube_SerialMacrocycleConvCube` added that along with the :ref:`cube_DGConformerCube` can be used to create
an efficient Macrocycle Omega floe that can generate macrocycle conformers in parallel mode.
* New :ref:`cube_OmegaFastRocsCube` added that generates conformers for FastROCS.
* New :ref:`cube_ConfAlignCube` added that can be used to align all conformers of a molecule based on RMSD.
* New :ref:`cube_OmegaFastRocsCube` added that generates conformers for FastROCS.
* New :ref:`cube_CopyRecordCube` added that can be used to generate multiple copies of a record.
* New :ref:`cube_TorsionLabelerCube` that identifies torsions on each input molecule and generates copies of the record.
with a different torsion labeled by an atom tag on each copy.
* New :ref:`cube_ActiveSiteLigandEnergyCube` added that can be used to calculate intermolecular energies between a
ligand and a protein.
* New :ref:`cube_BroodCube` added that performs Brood calculations.

Cube Updates
--------------------------------------------------------------------------------
* Added new parameters `mode` and `rotor_offset` in the :ref:`cube_OmegaCube`.
* Fixed a parameter mismatch issue with the :ref:`cube_AssignChargesCube`.
* The :ref:`cube_TorsionScanCube` is updated to use a single tag to identify a torsion, instead of four tags. It also
now has options to generate conformations using Omega prior to torsion driving or use a provided conformer ensemble. By
default the cube reports both the MMFF/Sheffield energy as well as the MMFF/PB energy for each output conformer.
* Fixed :ref:`cube_ZapLigandPrepCube` to fail gracefully when it cannot prep a molecule.
* The :ref:`cube_TorsionScanCube` is updated to use a signle tag to identify a torsion, instead of four tags.
* New parameters `in_strain_mol_field` and `out_strain_mol_field` added to the :ref:`cube_FreeformStrainCube` and the
:ref:`cube_FreeformEneStrainCube` cubes to enable using separate input conformers for freeform free energy and strain calculations.
* New parameter `use_inp_ens` have been added, to enable using user provided conformer ensemble in calculations, to the
following freeform cubes:

  * :ref:`cube_FreeformConfEntropyCube`
  * :ref:`cube_FreeformConfFreeEnergiesCube`
  * :ref:`cube_FreeformStrainCube`
  * :ref:`cube_FreeformEneStrainCube`
  * :ref:`cube_FreeFormConfPrepEnsembleCube`
  * :ref:`cube_FreeFormConfPreOptimizeCube`
  * :ref:`cube_FreeFormConfRemoveDupsCube`
  * :ref:`cube_FreeFormConfOptimizeCube`
  * :ref:`cube_FreeFormConfEnergiesCube`


v0.13.4 November 2018
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to use ``OpenEye-orionplatform==0.1.9`` and the `2018.Oct`_ release of the `OpenEye Toolkits`_

Cube Updates
--------------------------------------------------------------------------------
* Minor fix to :ref:`cube_OmegaCube` to preserve molecule titles.
* Speed up fix for many of the FreeFormConf cubes.
* Updated the Hermite shape cubes for compatibility with the `2018.Oct`_ release of the `OpenEye Toolkits`_


v0.13.3 September 2018
======================

General Notice
--------------------------------------------------------------------------------
* Upgrades to use ``OpenEye-orionplatform==0.1.7``

New Cubes
--------------------------------------------------------------------------------

* The following serial cubes have been added to merge fields from two records:

  * :ref:`cube_JoinCube`
  * :ref:`cube_MergeCube`

* The following serial logic cubes and their corresponding parallel version
  have been added to compare on string fields:

  * :ref:`cube_CompareStringFieldsCube`
  * :ref:`cube_CompareStringFieldToConstantCube`
  * :ref:`cube_CompareStringFieldToRegexCube`

Cube Updates
--------------------------------------------------------------------------------

* :ref:`cube_CalculateXLogPCube` now calculates XLogP on the neutral form of
  the molecule. A new parameter has been added to the cube, called
  :ref:`Keep Neutral Form<cube_param_CalculateXLogPCube_keep_neutral>`,
  that controls whether to keep the neutral form of the molecule after the calculation
  or restore the original input form.

v0.13.0 August 2018
===================

General Notice
--------------------------------------------------------------------------------
* Upgrades to use ``OpenEye-orionplatform==0.1.0``


.. _2018.Oct: https://docs.eyesopen.com/toolkits/python/releasenotes/releasenotes/index.html
.. _OpenEye Toolkits: https://docs.eyesopen.com/toolkits/python/index.html