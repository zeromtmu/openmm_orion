
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


NPT Simulation
--------------


NPT simulation of an OpenMM-ready System

Ex: python floes/openmm_MDnpt.py --system complex.oeb --nanoseconds 0.01

Parameters:
-----------
complex (file): OEB file of the prepared system

Optional:
--------
picosec (float): Number of picoseconds to warm up the complex
temperature (decimal): target final temperature in K
pressure (decimal): target final pressure in atm

Outputs:
--------
ofs: Outputs the constant temperature and pressure system


:bluebold:`Promoted Parameters`

   * | **nanoseconds**   (decimal) :  Length of MD run in nanoseconds 
     | *Default:* :blue:`0.01`  

   * | **temperature**   (decimal) :  Selected temperature in K 
     | *Default:* :blue:`300.0`  

   * | **pressure**   (decimal) :  Selected pressure in atm 
     | *Default:* :blue:`1.0`  

   * | **md_engine**   (string) :  Select the MD Engine 
     | *Default:* :blue:`OpenMM`  
     | *Choices:* :green:`OpenMM`, :green:`Gromacs`

   * | **trajectory_interval**   (decimal) :  Trajectory saving interval in ns 
     | *Default:* :blue:`0.0005`  

   * | **reporter_interval**   (decimal) :  Reporter saving interval in ns 
     | *Default:* :blue:`0.001`  

   * | **out**   (dataset_out) :  Output dataset to write to 

   * | **fail**   (dataset_out) :  Output dataset to write to 

   * | **system**   (data_source) :  System Input File - System input file 


