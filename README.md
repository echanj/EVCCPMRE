# EVCCPMRE

Author: Eric J. Chan
Repo create date: 1/3/2023 
Latest version upload: 23/11/2023

Short Description:

EVCCPMRE is a python based wrapper program to perform replica exchange MC for Crystal Structure Prediction using 
a modified version of the fortran program UPACK. EVCCP stands for "Extended Variables Coupled to Crystal Polymorphs".   

version info:
 - EVCCPMRE_UPACK_rigid_run.py : Original stable version of the wrapper that we used for inital trials.
   (it will only work in py27)

 - EVCCP_module_v8_rigid.py : latest py37 compatable version which has been setup to run on headless server.
   Inputs and log outputs are modified so that each bath now requires input displacement magnitudes and force constant. The algoritham has been modified so that the effect of on-the-fly updates of these inputs can be tested.         
 
