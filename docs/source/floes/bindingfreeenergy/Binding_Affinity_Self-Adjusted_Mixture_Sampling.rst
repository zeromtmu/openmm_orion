
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


Binding Affinity Self-Adjusted Mixture Sampling
-----------------------------------------------


NOTE: this is an Alpha Test version. 
We are actively working on improving Yank in Orion

The Absolute Binding Affinity Free Energy protocol (ABFE) performs Binding Affinity calculations
on a set of provided ligands and related receptor by using YANK Self-Adjusted Mixture Sampling
( http://getyank.org/latest/ ). The ligands need to have coordinates and correct chemistry.
Each ligand can have multiple conformers, but each conformer will be treated as a different ligand
and prepared to run ABFE. The protein needs to be prepared at MD preparation standard. This includes
capping the protein, resolve missing atoms in protein residues and resolve missing protein loops.
The parametrization of some "known unknown" non standard protein residues is partially supported.
Ligands need to be already posed in the protein binding site. A complex (Bonded State) is formed,
solvated and parametrized accordingly to the selected force fields. In a similar fashion the Unbounded
state is also prepared. Minimization, Warm up (NVT) and Equilibration (NPT) stages are performed
an the Bonded and Unbounded states. In order to minimize Molecular Dynamics (MD) failures along
these stages, positional harmonic restraints are applied on the ligand and protein with different
force constants. At the end of the equilibration stages the ABFE calculations are run by YANK with
the selected parameters. Calculated Binding Affinities for each ligand are output with the related
floe reports.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands
protein: Dataset of the prepared protein

Outputs:
--------
out : Dataset of the solvated systems with the calculated binding free energies and
floe reports


:bluebold:`Promoted Parameters`

   * | **ligands**   (data_source) :  Ligand Input File - Ligand file name 

   * | **iterations**   (integer) :  Total number of Yank iterations 
     | *Default:* :blue:`1000`  

   * | **temperature**   (decimal) :  Temperature (Kelvin) 
     | *Default:* :blue:`300.0`  

   * | **pressure**   (decimal) :  Pressure (atm) 
     | *Default:* :blue:`1.0`  

   * | **hmr**   (boolean) :  On enables Hydrogen Mass Repartitioning. NOTE:Not currently implemented in Gromacs 
     | *Default:* :blue:`False`  

   * | **restraints**   (string) :  Select the restraint types to apply to the ligand during the alchemical decoupling. Choices: harmonic, boresch 
     | *Default:* :blue:`boresch`  
     | *Choices:* :green:`harmonic`, :green:`boresch`

   * | **verbose**   (boolean) :  Yank verbose mode on/off 
     | *Default:* :blue:`False`  

   * | **protocol_sams**   (string) :  Select the sams protocol type 
     | *Default:* :blue:`windows_sams`  
     | *Choices:* :green:`auto_protocol`, :green:`windows_sams`

   * | **charge_ligands**   (boolean) :  Charge the ligand or not 
     | *Default:* :blue:`True`  

   * | **ligand_forcefield**   (string) :  Force field to be applied to the ligand 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **other_forcefield**   (string) :  Force field used to parametrize other molecules not recognized by the protein force field like excipients 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **protein_ff**   (string) :  Force field parameters to be applied to the protein 
     | *Default:* :blue:`amber99sbildn.xml`  
     | *Choices:* :green:`amber99sbildn.xml`, :green:`amberfb15.xml`

   * | **density**   (decimal) :  Solution density in g/ml 
     | *Default:* :blue:`1.03`  

   * | **salt_concentration**   (decimal) :  Salt concentration (Na+, Cl-) in millimolar 
     | *Default:* :blue:`50.0`  

   * | **protein**   (data_source) :  Protein Input File - Protein file name 

   * | **protein_prefix**   (string) :  Protein prefix used to identify the protein 
     | *Default:* :blue:`PRT`  


