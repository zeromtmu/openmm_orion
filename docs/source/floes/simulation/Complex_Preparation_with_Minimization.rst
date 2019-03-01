
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


Complex Preparation with Minimization
-------------------------------------


Complex Preparation Workflow

Ex. python floes/openmm_complex_prep.py --protein protein.oeb
--ligands ligands.oeb  --ofs-data_out complex.oeb

Parameters:
-----------
protein (file): OEB file of the prepared protein
ligands (file): OEB file of the prepared ligands


Outputs:
--------
ofs: Output file


:bluebold:`Promoted Parameters`

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **ligands**   (data_source) :  Ligand Input File - Ligand file name 

   * | **density**   (decimal) :  Solution density in g/ml 
     | *Default:* :blue:`1.03`  

   * | **salt_concentration**   (decimal) :  Salt concentration (Na+, Cl-) in millimolar 
     | *Default:* :blue:`50.0`  

   * | **steps**   (integer) :  Number of minimization steps.
                  If 0 the minimization will continue
                  until convergence 
     | *Default:* :blue:`0`  

   * | **md_engine**   (string) :  Select the MD Engine 
     | *Default:* :blue:`OpenMM`  
     | *Choices:* :green:`OpenMM`, :green:`Gromacs`

   * | **protein**   (data_source) :  Protein Input File - Protein file name 

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **charge_ligands**   (boolean) :  Charge the ligand or not 
     | *Default:* :blue:`True`  


