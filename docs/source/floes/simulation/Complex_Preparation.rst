
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


Complex Preparation
-------------------


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

   * | **protein_ff**   (string) :  Force field parameters to be applied to the protein 
     | *Default:* :blue:`amber99sbildn.xml`  
     | *Choices:* :green:`amber99sbildn.xml`, :green:`amberfb15.xml`

   * | **ligand_ff**   (string) :  Force field to be applied to the ligand 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **other_ff**   (string) :  Force field used to parametrize other molecules not recognized by the protein force field like excipients 
     | *Default:* :blue:`GAFF2`  
     | *Choices:* :green:`GAFF`, :green:`GAFF2`, :green:`SMIRNOFF`

   * | **protein**   (data_source) :  Protein Input File - Protein file name 

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **charge_ligands**   (boolean) :  Charge the ligand or not 
     | *Default:* :blue:`True`  


