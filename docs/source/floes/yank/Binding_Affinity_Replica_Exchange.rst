
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


Binding Affinity Replica Exchange
---------------------------------


NOTE: this is an Alpha Test version. 
We are actively working on improving Yank in Orion

The Absolute Binding Affinity Free Energy protocol (ABFE) performs Binding Affinity calculations
on a set of provided ligands posed in a receptor by using YANK Replica Exchange ( http://getyank.org/latest/ ).
The ligands need to have coordinates and correct chemistry. Each ligand can have multiple conformers,
but each conformer will be prepared and treated as a different ligand.
The protein needs to be prepared to MD standard: This includes capping the protein,
resolving missing atoms in protein residues and resolving missing protein loops.
The parametrization of some common non-standard protein residues is partially supported.
Though input separately from the protein, Ligands need to be already posed in the
protein binding site.
A bound complex is formed, solvated and parametrized according to the selected force fields.
The unbound state is similarly prepared. For both bound and unbound states the ABFE
calculation is preceded by minimization, warm up, and equilibration in the presence of
positional harmonic restraints.
The ABFE calculation is then run by YANK with the selected parameters.
The output floe report for each ligand contains the calculated binding affinity and health checks.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands
protein: Dataset of the prepared protein

Outputs:
--------
* out : Dataset of the solvated systems with the calculated binding free energies
* floe report : An analysis of the results for each ligand


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

   * | **protocol_repex**   (string) :  Select the Repex window schedule protocol 
     | *Default:* :blue:`windows_29`  
     | *Choices:* :green:`auto_protocol`, :green:`windows_29`

   * | **charge_ligands**   (boolean) :  Calculate ligand partial charges 
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


