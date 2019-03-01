
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


Minimize
--------


Minimize an OpenMM-ready solvated complex

Ex: python floes/openmm_prepMDminimize.py --system complex.oeb --ofs-data_out min.oeb --steps 1000`

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
steps (int): Number of MD steps to minimize the system. If 0 until convergence will be reached

Outputs:
--------
ofs: Outputs the minimized system


:bluebold:`Promoted Parameters`

   * | **steps**   (integer) :  Number of minimization steps.
                  If 0 the minimization will continue
                  until convergence 
     | *Default:* :blue:`0`  

   * | **md_engine**   (string) :  Select the MD Engine 
     | *Default:* :blue:`OpenMM`  
     | *Choices:* :green:`OpenMM`, :green:`Gromacs`

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **system**   (data_source) :  System Input File - System input file 


