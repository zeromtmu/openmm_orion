
.. raw:: html

   <style type="text/css">
      span.red { color: #cd4f39; }
   </style>

.. role:: red

.. raw:: html

   <style type="text/css">
      span.redbold { color: #cd4f39  ; font-weight: bold;}
   </style>

.. role:: redbold

.. raw:: html

   <style type="text/css">
      span.green { color: darkgreen; }
   </style>

.. role:: green

.. raw:: html

   <style type="text/css">
      span.greenbold { color: darkgreen; font-weight: bold;}
   </style>

.. role:: greenbold

.. raw:: html

   <style type="text/css">
      span.blue { color: darkblue; }
   </style>

.. role:: blue

.. raw:: html

   <style type="text/css">
      span.bluebold { color: darkblue; font-weight: bold;}
   </style>

.. role:: bluebold

.. raw:: html

   <style type="text/css">
      span.blacktitle { font-size: 175%; font-weight: 700;
        font-family: "Roboto Slab","ff-tisa-web-pro","Georgia",Arial,sans-serif; }
   </style>

.. role:: blacktitle


Short Trajectory MD
-------------------


NOTE: this is an Alpha Test version.
We are actively working on improving the MD sampling.

The Short Trajectory MD (STMD) protocol performs MD simulations given a set of
prepared ligands and a prepared protein as input.
The ligands need to have coordinates, all atoms, and correct chemistry. Each
ligand can have multiple conformers but each conformer will be run separately
as a different ligand.
The protein needs to be prepared to an MD standard: protein chains must be capped,
all atoms in protein residues (including hydrogens) must be present, and missing
protein loops resolved. Crystallographic internal waters should be retained where
possible. The parametrization of some common nonstandard residues is partially supported.
The STMD floe requires as inputs the protein and a set of ligands correctly posed
in the protein binding site. A complex is formed with each ligand and conformer
separately, and the complex is solvated and parametrized according to the
selected force fields. A minimization stage is peformed on the system followed
by a warm up (NVT ensemble) and three equilibration stages (NPT ensemble). In the
minimization, warm up and equilibration stages positional harmonic restraints are
applied on the ligand and protein. At the end of the equilibration stages a short
(default 2ns) production run is performed on the unrestrained system.

Required Input Parameters:
--------------------------
ligands (file): dataset of prepared ligands posed in the protein active site.
protein (file): dataset of the prepared protein structure.

Outputs:
--------
out:  OERecords (one per ligand) of MD and Analysis results.
floe report: html page of the Analysis of each ligand.


:bluebold:`Promoted Parameters`

   * | **ligands**   (data_source) :  Ligand Input File - Ligand file name 

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **protein**   (data_source) :  Protein Input File - Protein file name 

   * | **temperature**   (decimal) :  Temperature (Kelvin) 
     | *Default:* :blue:`300.0`  

   * | **pressure**   (decimal) :  Pressure (atm) 
     | *Default:* :blue:`1.0`  

   * | **hmr**   (boolean) :  On enables Hydrogen Mass Repartitioning. Not currently implemented in Gromacs 
     | *Default:* :blue:`False`  

   * | **md_engine**   (string) :  Select the MD available engine 
     | *Default:* :blue:`OpenMM`  
     | *Choices:* :green:`OpenMM`, :green:`Gromacs`

   * | **protein_prefix**   (string) :  Protein prefix used to identify the protein 
     | *Default:* :blue:`PRT`  

   * | **density**   (decimal) :  Solution density in g/ml 
     | *Default:* :blue:`1.03`  

   * | **salt_concentration**   (decimal) :  Salt concentration (Na+, Cl-) in millimolar 
     | *Default:* :blue:`50.0`  

   * | **protein_ff**   (string) :  Force field parameters to be applied to the protein 
     | *Default:* :blue:`amber99sbildn.xml`  
     | *Choices:* :green:`amber99sbildn.xml`, :green:`amberfb15.xml`

   * | **ligand_ff**   (string) :  Force field to be applied to the ligand 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **other_ff**   (string) :  Force field used to parametrize other molecules not recognized by the protein force field like excipients 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **prod_ns**   (decimal) :  Length of MD run in nanoseconds 
     | *Default:* :blue:`2.0`  

   * | **prod_trajectory_interval**   (decimal) :  Trajectory saving interval in ns 
     | *Default:* :blue:`0.002`  

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **charge_ligands**   (boolean) :  Charge the ligand or not 
     | *Default:* :blue:`True`  


