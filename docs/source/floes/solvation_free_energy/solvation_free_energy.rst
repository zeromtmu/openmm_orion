
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


Solvation Free Energy
---------------------


The Solvation Free Energy protocol performs Solvation Free Energy Calculations (SFEC) on
a set of input ligands using YANK ( http://getyank.org/latest/ ). The ligands need
to have coordinates, correct chemistry and must be neutral. Each ligand can have multiple
conformers, but each conformer will be prepared and treated as a different ligand.
The ligands are solvated in water (or other solvent or solvent mixture) and parametrized
by the selected force field.
Preceding the SFEC is minimization, warm up, and equilibration in the presence of
positional harmonic restraints. The SFEC is then run by YANK with the selected parameters.
The output floe report contains the Solvation Free Energy values and health checks.

Required Input Parameters:
-----------
ligands: Dataset of the prepared ligands

Outputs:
--------
* out : Dataset of the solvated systems with the calculated solvation free energies
* floe report : An analysis of the results for each ligand


:bluebold:`Promoted Parameters`

   * | **hmr**   (boolean) :  Hydrogen Mass Repartitioning 
     | *Default:* :blue:`False`  

   * | **iterations**   (integer) :  Total Number of Yank iterations for the entire floe. A Yank iteration is 500 MD steps 
     | *Default:* :blue:`1000`  

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **temperature**   (decimal) :  Temperature (Kelvin) 
     | *Default:* :blue:`300.0`  

   * | **pressure**   (decimal) :  Pressure (atm) 
     | *Default:* :blue:`1.0`  

   * | **Ligand ForceField**   (string) :  Force field to be applied to the ligand 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **verbose**   (boolean) :  Print verbose YANK logging output 
     | *Default:* :blue:`False`  

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **density**   (decimal) :  Solution density in g/ml - Solution Density in g/ml 
     | *Default:* :blue:`1.0`  

   * | **solvents**   (string) :  Solvent components - Comma separated smiles strings of solvent components 
     | *Default:* :blue:`[H]O[H]`  

   * | **molar_fractions**   (string) :  Molar fractions - Comma separated strings of solvent molar fractions 
     | *Default:* :blue:`1.0`  

   * | **ligands**   (data_source) :  Ligand Input File - Ligand file name 

   * | **charge_ligands**   (boolean) :  Calculate ligand partial charges 
     | *Default:* :blue:`True`  


