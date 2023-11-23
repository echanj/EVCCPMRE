#!/usr/bin/env python

'''
Author Eric Chan 2021


framework for working with co-crystals and rigid molecules 
it will bypass the need to create the cor.001 on the fly from a Z-matrix description 

 local dev version  

   - NOV 2023 - changing input so that each bath reads in a specified shift and spring con  

   - in version 8 com_shift eul_shift and kbias  are now arrays with an entry for each bath
   - they can be updated on-the-fly


  v8.0 
  
   - this version built to run via notebook style and then update the repo
     adapted directly from run_UPACK_PT_coumarin_module_v8_rigid.py which was the local version 
     running using python 2.7 

   - plan was to test finite difference bit to update shifts 
     but for now best to just fix shifts and temperatures as input variables 

   - for FEP shifts will be fixed 
     modifying code accordingly 
     disabled dumping of packfiles for FEP 

   - bugfix we needed to np.round() euler shifts prior to astype(int) 
     this affects smaller values of shifts  i.e. +/- 1 degree  

   - adding in scaling feature based on polynomial fit to the data 
     the polynomial is degree=2 and interaction terms only 
      see /Users/echanj/Work/upack/data/data_coumarin_test/server_jobs 
      despite it being fitted to coumarin data 
      it should translate ok but one should check 
      the subroutine is get_shifts() 
      this part of the routine is still isotropic assumption. this should be ok, MC occurs on gas phase     


   - turning off the shift factor scaling to test effect on temp, sprcon and aceptance ratios
      line 464 in run_upack_search() 
      line 869 

       -  subtley in the code 
          line 1103 :     eul_shift=().astype(float)  was originally  set to be .astype(int)
           reasoning is that integer values for the euler shift are used as input to UPACK
           should change this in a later c++ verison if it ever happens. 


   - when doing random search needed to change output comvar formating since too many decimal places
     added in condintional statement for this  

  v7.2 - initializing bath at a proper minimum causes non-ergodic behavior for strong coupling 
         adding option to accept inital shift otherwise we would be relying too much on acc2 
         line 870 and 1012


  v7.1 - upack is modifed to ouput the bias energy in the file ebias.17
         this is read to get proper energies since it is possible that  
         upack with output symetry equivelant trasformed eulers and com in the ecv.18
         this is ok for restarting but not for stepwise MC when now we are using very strong couplings
         so the trasformed coordinates appear very far from input.    
         still need to read in the coordinates from ecv.18
         but no longer need to read the structure energy from pack.10. 
         still need the skiprow number for  manage_PT_data_multiprocess_cocryst.py

  v7.0 - modified since coupling term with Kspr now affect the SD  and CG routines 
         of the structure generation. The energy in the pack.10 includes the coupling term.
         the couling term can still be evaluated from the output ecv.18 and thus can be seperated.     



  v6.1  - august 30 2021 - modified so that input for flexible search must include  torsions 
        for all molecules.  
        for a rigid search it should still only require one set of torsions in the input   
        becasue there is only one z-matrix.

        also the hydfix string is now only print out as seperate argument if not a random search  
        this is becase reading in alot of 1000.0 casues erorrs 

       - modified to run prep on the molecule to optimize the config prior to packing
         -  this si becasue shifts from the z-matrix  casues intra atomic clasehes during rigid searches,
             this shoudl prevent this by making sure we always start from a nice rigid molecule 
            this means the torsions wil corespond to the configuration prior to optimization.    



 v6.0 - 1st July 2021 - added in a hardwired option so that at each cycle 
       an added overrelaxtion update is perfroed whereby in each bath the CV's 
        for the MinE config will over write the current state. 
        right now we only print this info on the fly. we want to use this each cycle 
        because the system is intrinsically non-ergodic, we think this will help 
        the CV states diffuse better becasue certain "un-ideal" states are prone to an 
        over-rejection error. An ideal state has samll varance even at low T.        
        this is done by adding in the ecv.18 variables which are optionaly printed to the PT return pipes 
        as eulvar2, comvar2. 
        so the over relaxation update can occur just at the end of the cycle 
        these are tied to the current minE_old config for each bath and wont actually partake in any exchange directly.
        The idea is that each bath is monitoring minE and eulvar2 and comvar2.  
        and the importance is to maintian adiabaticity bettwen polymorph and ref config. within the ith steps of each cycle.
        for over relaxation we can imageing an inner MCloop that instantly forces the refernce confg. to adopt the minE config. 
        as though we have drastcally increased the force constant for the steps of our imaginary MC loop.
        this is anecdotally just like asking neuclar positions to quickly update to some given electron density.    


 v5.102 patch made to allow use of placing a STOP file in a search folder that is taking a long time  


 v5.101 server  version
  May 28 - fixed a critical bug in how the echange move occurs. it appears the coordinates were not properly 
         swaping if the normal metropolis MC move was rejected. see lines 871 onwards
         modified so that the current configuration is properly updated
         this version has alot more print statements
          it is essentailly the same but better for diagnosing errors 

 v5.1
May 13 2021 - added hardwired option so that only one molecue can change per step
              the aim is that at the C.V. will always reduce the bias energy 


 v5.0
April 22 2021 - adjust algoritham so that the choice of minE is based on U(X,H,S) and not just U(X,H)
                as this ensures coupling is properly taken into account. The idea is that the minE config will 
                be as close as possible to the corresponding refercne system config and not just the minE 
                from each batch which. we must then examine the effect of the spring const.
              - adjust output logs to include the actual final system value of minE, aveE and bias (i.e. not just the test value)
              - adjust options menu to enable/disable metaV and also multi bath swaps.       
 


 v4.0
Mar 24 2021 - updated control of range for com coordinates - now specify cellmin and cell max
              note: it is possible to get a good idea of what ranges are suitable by perfroming a preliminary 
                    pack12 run with a set DFIX and PRINT 4. 
                    you can then grep the output for the listing of trial unitcell values.
                    i.e. look for the lines begining with " Trial values for dexp, vexp, ax, by, cz ...."  
                     
                     

Mar 16 2021 - introduce a history depedent biasing potential (metaV) which is similar to 
                Wang-landau scheme but using gaussian kernals.
                this get updated at the end of each step due to parellization 
            - note : there were alot of other bug fixes - the bias calcualtion part was not reading in the correct corrdinate
                    this actually affects the bias calcualtion but not aveE or minE
                     it was crtical for correctly evaluating metaV 
            - now allowing for nbath/2 exchange moves per step i.e. this better allows possibility that low energy strcture will
                 be able to fall into lower T bath, since number of steps in this algoritham is limited.  

 v3.0
Mar 16 2021 - now just testing  MinE and bias between referece CV state and MinE state.
               this is beacsue the average energy and biasing over all samples in the mini-batch 
               has high variance    



Mar 11 2021 -  hardwire mod so that option to drive/monitor eulvar and comvar instead of eulvar2 and comvar2 
            - add back in the aveage energy offset 

Feb 24 2021 -  this code abstracted from the xxix project and designed to run 
               on any molecule without 
               the Z-matrix capability for torison. i.e. for now just rigid body search then option for flexible optimization
               modified code so that it will use the minE+bias as the basis for
               the rejection scheme (the bias here is an energy penalty from the minE config.) 

 v2.0 
Feb 1 2021: modified the code to update extended variables from minimum energy of the distributino 
            this happens inthe read_pack_files() function
            note: this modification has onluy been made on the eulers and the COM 
                  it is not currently implemented for the torsions 
            - should also monitor minimum energy

jan 20 2021: added fix to scale shifts by kB*T
Jan 19 2021: bugfix made to mean energy evaulation - big problems when sampling for small lots 
             where there is no relevant information. the average energy is very high 
             for most FF the case where ene<0 has some information. so this is preseved.  
             - fix up Kb factor for both PT and metropMC moves - double check working correctly    

modified for even high z' but also  increased stde and conj grad for pack.10 generation
it will take longer to generate the pack 10 

modified version for z'=4 on xxix 1/13/21
# hydfix %s # removed for xxix 
# no need for Z-matrix handling since rigid body 

Simualted annealing version of the previous MC module 
- allows for multiple baths
- ramping cycles are linear. each ramp will go up and down 
- we can use multiprocessing to paralleize multiple baths 

note: the dihedral output during PACK12 may not be the same value as that specified for the Z-matrix
      but the Z-matrix is controling the same torsion.
       eg. for ethanol the Z-matrix controls torsion-H9 O3 C2 H8 which is defaulted to 176.95 degrees 
           but UPACK controls the same torsion with H9 O3 C2 C1  which is -62.965 degrees

     when using torsion as a CV it is best to run in flexible search mode with hydfix -1000
       but you shoudl always use hydfix -1000 in rigid molecule search mode if you have specified any torsions.   

 update 1/5/21 - ypack12 has been modified so that if option "PRINT > 1" then the file ecv.18 will be output to
the folder and contain information on com,eulers and dihedrals. In theory we must use this to re-wieght 
 the energy of the distributions , thus the read_pack_files has been modified now to process the biasing energy.    

 update 1/6/21 - removed unessesary unit cell data extraction from read_pack_files() which was useful in earlier MC versions
                  but only worked for monoclinic space groups
                  
 update 1/7/21 - search allow for CV lines in the inputs to aslo be counted as random variables, this is for eg. in the case
                 Z'=2 and we only want to tether one set of euler and com for biasing.
                 
'''

import re
import sys
import os

import numpy as np
from numpy import linalg as la
import subprocess 

from shutil import copyfile, copy2
from shutil import move,rmtree

import multiprocessing
import time

# import openbabel

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# import scipy
# import matplotlib


def get_shifts(eu,co,Tf,Ks):

# this is the model from fitting 
 w0=3.938529191636282 
 W=np.array([-9.12320281e-01, 5.42184012e+01, -3.74806002e-01, -2.50966109e-04,
              3.84889514e+00, 6.11086728e-02, -1.51648423e-05, -5.92904311e+00, 
              1.29726400e-02,  9.99481459e-06])
 print ("input  eul,com shifts  for T=%.1f : %.6f %.6f"%(Tf,eu,co))

 XIN= np.array([eu, co, Tf, Ks, eu*co, eu*Tf, eu*Ks, co*Tf, co*Ks, Tf*Ks])
 Ediff_pred = w0 + np.dot(W,XIN)   
 
 Pacc=np.exp(-Ediff_pred/Tf)
 print ("initial predicted ediff Pacc: %.6f %.6f" %(Ediff_pred,Pacc))
 
 scw0=1.0
 if Pacc<0.5 : 
 # print ("uniform scaling of input shifts for the given temperature") 
  while np.all([Pacc<0.5,scw0>0,eu>0.51,co>0.00001])  : 
    scw0-=0.001
    eu=eu*scw0 
    co=co*scw0 
    XIN= np.array([eu, co, Tf, Ks, eu*co, eu*Tf, eu*Ks, co*Tf, co*Ks, Tf*Ks])
    Ediff_pred = w0 + np.dot(W,XIN)   
    Pacc=np.exp(-Ediff_pred/Tf)
    # print (scw0,eu,co,Ediff_pred, Pacc)
    
 # print ("scaling completed") 
  print ("scaled eul,com shifts, predicted Ediff and Pacc  for T=%.1f : %.6f %.6f %.6f %.6f "%(Tf, eu, co, Ediff_pred, Pacc))

 elif Pacc>0.5 : 
 # print ("uniform scaling of input shifts for the given temperature") 
  while np.all([Pacc>0.5])  : 
    scw0+=0.001
    eu=eu*scw0 
    co=co*scw0 
    XIN= np.array([eu, co, Tf, Ks, eu*co, eu*Tf, eu*Ks, co*Tf, co*Ks, Tf*Ks])
    Ediff_pred = w0 + np.dot(W,XIN)   
    Pacc=np.exp(-Ediff_pred/Tf)
    # print (scw0,eu,co,Ediff_pred, Pacc)
    
 # print ("scaling completed") 
  print ("scaled eul,com shifts, predicted Ediff and Pacc  for T=%.1f : %.6f %.6f %.6f %.6f "%(Tf, eu, co, Ediff_pred, Pacc))
 else: # Pacc==0.5
  print ("scaled eul,com shifts, predicted Ediff and Pacc  for T=%.1f : %.6f %.6f %.6f %.6f "%(Tf, eu, co, Ediff_pred, Pacc))

 return eu,co





def gaussNd(S,S0,w,SIG):
 # S is the position on the bias-potential we are evaluating  
 # S0 is a vector of coordinates for the center of an N-d gaussian with corresponding widths (SIG)  
 X=((S-S0)**2.0)/(2.0*SIG**2.0)  # linearized
 return w*np.exp(-np.sum(X))

def metabias(eulvar_hills,comvar_hills,eulfix,comfix,zpr,nmpg):

 w=1.0 # gaussian hieght in kJ/mol 
 eulsig=np.ones(3*zpr*nmpg)*10.0 # gaussian widths for eulers in degrees
 comsig=np.ones(3*zpr*nmpg)*0.2 # gaussian widths for coms in nanometers
 all_sig=np.r_[eulsig,comsig]
 all_var=np.r_[eulfix,comfix]
 nst,nbs,nvar = np.shape(eulvar_hills)  # should be current number of steps, baths, vars 

 if not nvar==nmpg*zpr*3: 
  print ("Error number of euler variables %i incorrect for this Z'. should be %i  " %(nvar,nmpg*3*zpr))
  exit()

 all_hills=np.c_[np.array(eulvar_hills).reshape((-1,3*zpr*nmpg)),np.array(comvar_hills).reshape((-1,3*zpr*nmpg))]

 #  each row of all_hills is an N-dimenisonal gaussian we need to evaluate 
 Vmb=0.0
 for hill in all_hills:
  Vmb=Vmb+gaussNd(all_var,hill,w,all_sig)

 return Vmb

def read_pack_files(wdir_path,sr,zpr,nmpg,eulvar,comvar,dihed,MC_other,Kbias,eulvar_hills,comvar_hills):

 isFlex = MC_other[3]
 # if isFlex==0:
 # ntor=len(dihed) # number of torsions per molecue
 # dihed=(np.zeros((zpr,ntor))+dihed).reshape(1,-1)[0]  
 #else: 
 #ntor= len(dihed)/zpr 

#  dihed=dihed+di_offset # convert to upack setting
#  if doing rigid body the torsion bias from each step does not differ between structures
#  we can keep track of the angles each step for metaV
 
# print "in read pack files"
# check where variables are random and have no effect on bias
# eulerIsRand=np.where(MC_invar[0][:]>999)
# comIsRand=np.where(MC_invar[2][:]>999)
 eulerIsRand=np.where(eulvar>999)
 comIsRand=np.where(comvar>999)
# torIsRand=np.where(dihed>999)
# torIsFix=np.where(dihed<-999)
 # above also fixs up for the case of random or fixed Torsions
 #  this is so the bias will not be affected by magnitude of deviation from fixed input values. 
 # i.e.  assume fixed (eg. dimeric H-bond ) interactions are less likely to affect chocie of minimum    

 eulfix=np.array(eulvar) # (MC_invar[0][:]).astype(float)
 comfix=np.array(comvar) # (MC_invar[2][:]).astype(float)
 
# evaualte history dependant bias  using gaussian kernel
#  this conditional will aloow us to completly disable the meta potential withnin the algoritham
#  becasue this get slower and slower as the algoritham evolves  
 if np.shape(eulvar_hills)[0] > 0: 
  metaV=metabias(eulvar_hills,comvar_hills,eulfix,comfix,zpr,nmpg)
 else: 
  metaV=0.0

# out12name=wdir_path+'out12'
#with open(out12name) as file:
#       for line in file:
#         if re.search('Dihedrals from COR file', line):
#          line = line.replace("\n", " ")   # get rid of newline 
#          uptorfix=np.array(line.split()[-ntor:]).astype(float) # the value of dihed that corresponds to COR.001 file
#          break
 
 nR=MC_other[0]
 eng=[]     
 Ecv_com2=[]
 Ecv_eul2=[]
 Centers=[]                     # each column is molecule
 Eulers=[]
 # Dihedrals=[]

 packname=wdir_path+'ebias.17'
 fhpack = open(packname,'r')
 dat = fhpack.readlines() 
 fhpack.close()

 ecvname=wdir_path+'ecv.18'
 fhecv = open(ecvname,'r')
 cvdat = fhecv.readlines() 
 fhecv.close()
 
 cnl=zpr+(nmpg-1)  # differnt structure entry every cnl lines 
 moldat=[]
 # tordat=[]
 #for z in range(zpr): 
 # tordat.append(cvdat[(nmpg*zpr)+z::cnl]) 
 for z in range(zpr*nmpg): 
  moldat.append(cvdat[z::cnl]) 
 

 for n in range(nR): 
  eng.append(float(dat[n].split()[5])) 
  Ecv_com2.append(float(dat[n].split()[1])) 
  Ecv_eul2.append(float(dat[n].split()[3]))
  Centers.append([np.array(moldat[z][n].split()[3:6]).astype(float) for z in range(zpr*nmpg)]) 
  Eulers.append([np.array(moldat[z][n].split()[7:10]).astype(float) for z in range(zpr*nmpg)]) 
  # Dihedrals.append([tordat[z][n].split()[-ntor:] for z in range(zpr)]) 

 eng=np.array(eng).astype(float)      
 Ecv_com2=np.array(Ecv_com2).astype(float)      
 Ecv_eul2=np.array(Ecv_eul2).astype(float)      
 Centers=np.reshape(Centers,(nR,nmpg*zpr*3))
 Eulers=np.reshape(Eulers,(nR,nmpg*zpr*3))
 # Dihedrals=np.reshape(Dihedrals,(nR,ntor*zpr)).astype(float)

# print eng
# print Ecv_com2
# print Ecv_eul2
# print Centers
# print Eulers
# print Dihedrals
 

# v5.0 modify this so that the choice of minE now takes into acount the biasing  

# Bias1 -  reletive to unbiased minE configuration as refernece - currently disabled for v5.0 

# COMfix=np.zeros((nR,zpr*3))+Centers[idx_minE]  # this is the V2 modification, we no longer use  MC_invar
# EULfix=np.zeros((nR,zpr*3))+Eulers[idx_minE]

# TORfix=np.zeros((nR,zpr*ntor))+np.vstack([uptorfix for i in range(zpr)]).flatten()
# Ecv_com=Kbias[0]*0.5*(Centers-COMfix)**2.0     
# Ecv_eul=Kbias[1]*0.5*(1.0-np.cos(np.radians(Eulers-EULfix)))   
# Ecv_tor=Kbias[2]*0.5*(1.0-np.cos(np.radians(Dihedrals-TORfix)))  
# for idx in comIsRand: Ecv_com[:,idx]=0
# for idx in eulerIsRand: Ecv_eul[:,idx]=0     # if variable >999 i.e. random sampling then set bias component to zero  
# Ecv_com=np.sum(Ecv_com,axis=1)
# Ecv_eul=np.sum(Ecv_eul,axis=1)
# Ecv_tor=np.sum(Ecv_tor,axis=1)

# Bias2 - reletive to Extended Variable configuration

# COMfix2=np.zeros((nR,nmpg*zpr*3))+comfix   # this is from original version. takes into account the test variable
# EULfix2=np.zeros((nR,nmpg*zpr*3))+eulfix
# TORfix2=np.zeros((nR,ntor*zpr))+dihed
# Ecv_com2=Kbias[0]*0.5*(Centers-COMfix2)**2.0     
# Ecv_eul2=Kbias[1]*0.5*(1.0-np.cos(np.radians(Eulers-EULfix2)))   
# Ecv_tor2=Kbias[2]*0.5*(1.0-np.cos(np.radians(Dihedrals-TORfix2)))  
 if np.any(comIsRand==0): Ecv_com2==np.zeros(nR)
 if np.any(eulerIsRand==0): Ecv_eul2==np.zeros(nR)     # if variable >999 i.e. random sampling then set bias component to zero  
# for idx in torIsRand: Ecv_tor2[:,idx]=0
# for idx in torIsFix: Ecv_tor2[:,idx]=0
# Ecv_com2=np.sum(Ecv_com2,axis=1)
# Ecv_eul2=np.sum(Ecv_eul2,axis=1)
# Ecv_tor2=np.sum(Ecv_tor2,axis=1)

 if isFlex==0:
  ubias=Ecv_com2+Ecv_eul2  # dont include tor for flex=0
 else:
  ubias=Ecv_com2+Ecv_eul2+Ecv_tor2  

 engR=eng+ubias  #  refernce min energies in v7.1 add the bias  to the energy i.e. eng is unbiased energies  
# v7.1 modify this so that the choice of minE takes into acount the biasing.  
 idx_minE=np.where(engR==np.min(engR))[0][0] # figure out which structure had minE in the distrib of biased energies.
 minE=eng[idx_minE]# this is the unbiased minE that corresponds with min(eng)
 
# Bias3 - CV's now reletive to just MinE config only - single points 

# Ecv_com3=Kbias[0]*0.5*np.array(Centers[idx_minE]-comfix)**2.0     
# Ecv_eul3=Kbias[1]*0.5*(1.0-np.cos(np.radians(np.array(Eulers[idx_minE]-eulfix))))   
 Ecv_com3=Ecv_com2[idx_minE]
 Ecv_eul3=Ecv_eul2[idx_minE]
# Ecv_com3=np.sum(Ecv_com3)
# Ecv_eul3=np.sum(Ecv_eul3)
# Ecv_tor3=Ecv_tor2[idx_minE]   

# To obtain biasing energy we will use harmonic spring for centers 
# and harmonic cosine term for eulers and dihedrals  

 eng_fix=eng[eng<0.0]  # make sure energies below a cutoff  (default zero)

# ave_bias = np.mean(Ecv_com+Ecv_eul+Ecv_tor)  # with torsion 
# ubias=Ecv_com+Ecv_eul+Ecv_com2+Ecv_eul2  #  # ths is where we control which parts contribute towards the bias
# ubias=Ecv_com2+Ecv_eul2  #  Use the original bias which tethers to the CV
# ubias_fix=ubias[eng<0.0] 
 if isFlex==0:
  ubias_fix=Ecv_com3+Ecv_eul3 # the values here are only based on the minE structure which should have minE < 0 
  cfg_dihed=dihed
 else:
  ubias_fix=Ecv_com3+Ecv_eul3+Ecv_tor3
  cfg_dihed=Dihedrals[idx_minE]

# print Ecv_com3,Ecv_eul3,Ecv_tor3,ubias_fix

 if len(eng_fix)==0:
  ave_eng,var_eng,ave_bias=0.0,0.0,0.0 # this should reject the move if we didnt sample anything significnt
 else:
  ave_bias = np.mean(ubias_fix)         
  ave_eng = np.mean(eng_fix)
  var_eng = np.var(eng_fix)

 return ave_eng,var_eng,ave_bias,Eulers[idx_minE],Centers[idx_minE],cfg_dihed,minE,metaV 


def zmat2corr(rootname,targetdir,dfix,shift,nmpg):

 torsions=np.zeros(np.shape(dfix))

 tnat=0 # total number of atoms
 for nm in range(nmpg):
  zmfhin = open(rootname+'_'+str(nm+1)+'.zmat.up','r')  # use fixed label for zmatrix file  
  zdat = zmfhin.readlines()
  zdat= [re.sub(r'("(?:[^"]+|(?<=\\)")*")|#[^\n]*','',line) for line in zdat] # get rid of any commentry after '#'
  zmfhin.close()
  
  qxyz= np.array(zdat[-2].split()).astype(float)
  cxyz= np.array(zdat[-1].split()).astype(float)
 
  fhzmat=open(targetdir+'new%i.zmat'%(nm+1),'w')
  fhzmat.write(zdat[0])       # write the title  
  fhzmat.write(zdat[1])       # write the number of atoms 
  tnat=tnat+int(zdat[1])        # update total for final corr.001  
  for line in zdat[2:-2]:     # writeing the zmatrix 
   ztor= float(line.split()[-1]) # check if we have fvar on torsions  
   if ztor > 1000.0 :  
    torID=int(ztor-1000)-1 
    itor=dfix[torID]+((2*np.random.rand()-1.0)*shift[torID])  
    itor= itor+360 if itor < -180 else itor                                     
    itor= itor-360 if itor >= 180 else itor                                   
    torsions[torID]=itor 
    zstr=line.split()
    fhzmat.write("%s         %s   %s    %s   %s    %s  %.5f\n"%(zstr[0],zstr[1],zstr[2],zstr[3],zstr[4],zstr[5],torsions[torID]))
   else: 
    fhzmat.write(line)
  fhzmat.close()
 
  z2xlog = open(targetdir+'coords%i.log'%(nm+1), 'w')
  z2xrun = subprocess.Popen(['zmat2xyz','--q=%.6f,%.6f,%.6f,%.6f'%(qxyz[0],qxyz[1],qxyz[2],qxyz[3]),
                                      '--com=%.6f,%.6f,%.6f'%(cxyz[0],cxyz[1],cxyz[2]),'new%i.zmat'%(nm+1)],cwd=targetdir, stdout=z2xlog)
  z2xrun.wait()
  z2xlog.close()
 
 
 fhcor=open(targetdir+'cor.001','w')
 fhcor.write(' %s\n'%rootname)
 fhcor.write('%5i\n'%tnat)

 anum=0
 for nm in range(nmpg):
  coords = open(targetdir+'coords%i.log'%(nm+1),'r')
  coords_dat = coords.readlines() 
  coords.close()
  nat=int(coords_dat[0])
  for a in range(nat):
   anum+=1
   alabel=coords_dat[a+2].split()[0]
   carts=np.array(coords_dat[a+2].split()[1:4]).astype(float)
   fhcor.write("%5i%8s %6i%8.4f%8.4f%8.4f\n" %(nm+1,alabel,anum,carts[0]*0.1,carts[1]*0.1,carts[2]*0.1))    # the FFtype is just a dummy      
 fhcor.close()
 # os.remove(targetdir+'coords.log')
 # os.remove(targetdir+'new.zmat')

 return torsions   

def run_upack_search(rootname,label1,sgr,zpr,nmpg,label2,MC_invar,MC_other,DenExp,Tfactor,Kbias):

 EnablePackFiles=False # v8.0 hardwired option that enables writing of packfiles   
 isONEMOL=False    #  v5.1 hardwired option which only shifts one molecule each step 
 runPPOPT=False     #  v6.1 hardwired option to run preopt of the cor.001 output from zmat2cor()
 Tfactor=1.0        #  v7 dev onwards hardwired option so that shifts do not scale with input Tfactor.  

 UPACK_path='/home/ubuntu/upack/scripts/'
# UPACK_path='/Users/echanj/Work/upack/scripts/'
 mccyc = MC_other[1]
 mcstep = MC_other[2]
 isFlex = MC_other[3]

# define path names 
 dir_path=rootname+'_'+label1+'_'+sgr+'/' 
 wdir_path=dir_path+label2+'/' 
 packfile_path=dir_path+'packfiles/' 
 
# check to see if any params should be treated as random variable i.e. > 999  
 eulerIsRand=np.where(MC_invar[0][:]>999)
 comIsRand=np.where(MC_invar[2][:]>999)
 
# print "in bath %s"%label2 
# print MC_invar[1]
 # construct pack12 input file 
 # sort out shifts on mc params 
 np.random.seed()   # seeding generator on each process
 if isONEMOL : 
  molidx=np.random.randint(zpr*nmpg) # select the molecule 
  eulshifts=np.zeros(nmpg*zpr*3)
  eulshifts[(molidx-1)*3:molidx*3]=(2.0*np.random.rand(3)-1.0)*Tfactor*MC_invar[1][(molidx-1)*3:molidx*3]
  eulvar= (MC_invar[0][:]+np.round(eulshifts)).astype(int)
 else:
  eulvar= (MC_invar[0][:]+np.round((2.0*np.random.rand(nmpg*zpr*3)-1.0)*Tfactor*MC_invar[1][:])).astype(int)
 # make sure eulvar contains numbers  between 0-360
 eulvar= np.array([ev+360 if ev<0 else ev for ev in eulvar])                                    
 eulvar= np.array([ev-360 if ev>=360 else ev for ev in eulvar])                                    
# print eulvar
# also need to keep within certain box extents 
# in UPACK if generated randomly the COM is e.g. ax*rand(),by*rand(),cz*rand()
 cellmax=12.0   # maximum com hardwired variable  
 cellmin=-1.0   # minimum unitcell hardwired variable  
 if isONEMOL : 
  comshifts=np.zeros(nmpg*zpr*3)
  comshifts[(molidx-1)*3:molidx*3]=(2.0*np.random.rand(3)-1.0)*Tfactor*MC_invar[3][(molidx-1)*3:molidx*3]
  comvar= (MC_invar[2][:]+comshifts).astype(float)
 else:
  comvar= (MC_invar[2][:]+((2.0*np.random.rand(nmpg*zpr*3)-1.0)*Tfactor*MC_invar[3][:])).astype(float)

 comvar= np.array([com-(cellmax-cellmin) if com > cellmax  else com for com in comvar])                                    
 comvar= np.array([com+(cellmax-cellmin) if com <= cellmin  else com for com in comvar])                                    
 nrand = MC_other[0]
 
 for idx in eulerIsRand: eulvar[idx]=1000
 for idx in comIsRand: comvar[idx]=1000

 rseed = np.random.randint(100000,300000)+ord(label2)   # adds on a bit to doubly make sure seed is different between baths 

#####################
# deal with torsions
###################

# ONETOR=False
 # these values are hardwired and traslate the zmatrix setting into upack setting
 # di_offset=np.array([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])  # XXXII
 # xxxb
#di_offset=np.array([  0.00000000e+00,  -1.00000000e-01,   0.00000000e+00,
#                    -1.00000000e-01,   1.72800000e+02,   0.00000000e+00,
#                     0.00000000e+00,   0.00000000e+00,   0.00000000e+00,
#                     1.71600000e+02])
 Tor_fix = MC_invar[4][:]
 Tor_shift = MC_invar[5]*Tfactor

 # option isFlex off
 # we are controling the rigid config with a cor.001
 # only one set of CV is tracked 
 if isFlex==0 : 
 # generate a cor.001 file on the fly 
  torIsRand=np.where(Tor_fix>999)
  torIsFix=np.where(Tor_fix<-999)
  # manually bipassing zmatrix calc  #  dihed=zmat2corr(rootname,wdir_path,Tor_fix,Tor_shift,nmpg)
  dihed=Tor_fix.copy()
  copy2('./cor.001',wdir_path)
  # finished manual bypass
  for idx in torIsRand: dihed[idx]=1000
  for idx in torIsFix: dihed[idx]=-1000
  zp_dihed=dihed # +di_offset
  up_dihed=(np.zeros((zpr,len(Tor_fix)))+zp_dihed).reshape(1,-1)[0]  # this is the input value read by UPACK

# in this modified 6.1 verison - choice optimizing the flexibe molecule 
#  to generate a geometry optimized version of the cor.001 to prevent clashes in the molecule 
#  from rnadom shifts 

  if runPPOPT:

   setpp='''cutoff 2.4 3.4 3.4 5.0 5.0 0.0 0.0
flex 1
'''
   
   fhoutpp=open(wdir_path+"%s.pp" %(rootname),'w')
   fhoutpp.write(setpp)
   fhoutpp.close()
   
   pplog = open(wdir_path+'outpp', 'w')
   pprun = subprocess.Popen([UPACK_path+'runpp','s',rootname,'cor.001',rootname  ],cwd=wdir_path, stdout=pplog)
   pprun.wait()
   pplog.close()
   
   copyfile(wdir_path+'0000.mdi',wdir_path+'cor.001')

 
 else :  # isFlex on
  # the torsion will be monitored in upack formating 
  # copy over the coords for the molecule - use instead of zmatrix if rigid
  copy2('./cor.001',wdir_path)
  # generate the cor.001 file on the fly 
  Tor_fix = MC_invar[4]
  Tor_shift = MC_invar[5]*Tfactor
  # print Tor_fix
  # print Tor_shift

  # check if the input array has size multiplied by zpr 
 # if len(Tor_fix)/zpr < len(Tor_shift): 
  # dihed=Tor_fix # +di_offset
  # zp_dihed=(np.zeros((zpr,len(Tor_fix)))+dihed).reshape(1,-1)[0]  # this is the input value read by UPACK
  #else:
  zp_dihed=Tor_fix

  torIsRand=np.where(zp_dihed>999)
  torIsFix=np.where(zp_dihed<-999)

  # sort out the shift arrays
  # if ONETOR:
  dshifts=(2*np.random.rand(len(Tor_shift))-1.0)*Tor_shift  
  #  dshifts=(np.zeros((zpr,len(Tor_shift)))+dshifts).reshape(1,-1)[0]  
  #else: 
  #  zp_shifts=(np.zeros((zpr,len(Tor_shift)))+Tor_shift).reshape(1,-1)[0]  
  #  dshifts=(2*np.random.rand(zpr*len(Tor_shift))-1.0)*zp_shifts  

  zp_dihed=zp_dihed+dshifts

  zp_dihed= np.array([di+360 if di<180 else di for di in zp_dihed])                                    
  zp_dihed= np.array([di-360 if di>=180 else di for di in zp_dihed])                                    
  for idx in torIsRand: zp_dihed[idx]=1000
  for idx in torIsFix: zp_dihed[idx]=-1000
  up_dihed=zp_dihed

 str_hydfix = ["%.1f" % h for h in up_dihed] 
 str_hydfix = ' '.join(str_hydfix) 

 setpack12='''cutoff 1.2 3 3 2 2 0.001 0.001
spgr %s
ngr %i
rand %i
flex %i
dexp %.4f 1.2
angles -0.5 0.5 -0.5 0.5 -0.5 0.5
kspr %i %i
eulfix %s
grvfix %s
emax 40 1000
ebig 80000
seed %i
print 2
ntrymax 50000
'''

#
# axes 0.80 0.90 0.87 1.07 2.55 2.75  
# angles -0.5 0.5 -0.5 0.5 -0.5 0.5
# hydfix %s
# intra 0
# cutoff 1.5 2.4 2.4 3.0 3.0 0.0 0.0
# cutoff 3.0 5 5 5 5 0.001 0.001
# pres 1.0
# hydfix %s # removed for xxix 
# additional options below
# cutoff 1.2 3 3 2 2 0.001 0.001 (current default for PT mode)
# cutoff 1.2 3 3 3 3 0.0 0.0
# cutoff 0.8 3 3 2 2 0.02 0.02
# eulfix 217 84 10 322 113 252 
# axes 0.4 0.6 0.5 0.8 0.7 0.9   # better to not use this and let cells be as random as possible, use for testing only
# grvfix 0.014  0.163  0.118  -0.106  0.152  0.334
# hydfix  -63.4  179.1
# dexp 1.15 1.2    - when running PT it is always much faster to pre-calulate the expected density - we can do this in the init stage 
# dexp -20 1.2     - it appears the more structures you volume check it helps remove any issue with high energy artifacts in the dist   
# ebig 100000
# emax 40 10000
# defailt settings for SD and CG used in  pack3 for pack 12 they are usually smaller 
# but for higher z'  we want to ensure configs are true minima
# sd 100    0.01       0.0001   5   0.02
# cg 5000  0.0005  0.005  0.3  0 .02

 str_eulfix =' '.join(eulvar.astype(str))
# str_grvfix = ["%.6f" % c for c in comvar]
 str_grvfix = ["%.3f" % c for c in comvar] #  truncated decimal presion for inital cond  
# str_hydfix = ' '.join((np.ones(zpr*len(MC_invar[4]))-1001).astype(int).astype(str)) # the output here can be controled by isFlex

# change the string output formating if random variables
 if np.any(comIsRand): str_grvfix = ["%.1f" % c for c in comvar]
 str_grvfix = ' '.join(str_grvfix)

 # these are set to be the same for now
 Kspr_com=Kbias 
 Kspr_eul=Kbias

 fhout=open(wdir_path+"%s.pa12" %(rootname),'w')
 fhout.write(setpack12 %(sgr,zpr,nrand,isFlex,DenExp,Kspr_com,Kspr_eul,str_eulfix,str_grvfix,rseed))
 if np.shape(torIsRand)[1] != len(zp_dihed) and isFlex==1: fhout.write("hydfix %s\n"%(str_hydfix))
 fhout.close()

 # exit()

 print ("SA ramp cycle %i step %i running PACK12 in %s with seed %i" %(mccyc,mcstep,wdir_path,rseed))
 pack12log = open(wdir_path+'out12', 'w')
 pack12run = subprocess.Popen([UPACK_path+'runpa12','s',rootname,rootname  ],cwd=wdir_path, stdout=pack12log)
 pack12run.wait()
 pack12log.close()
# print "done"


# make a backup of the pack.10 and pack.19 for that step 
# dont do this for Diagnostic calcs eg FEP
 if EnablePackFiles:
        copyfile(wdir_path+'pack.10',packfile_path+'bath_%s_cyc_%i_step_%i.pack.10'%(label2,mccyc,mcstep))
        copyfile(wdir_path+'pack.19',packfile_path+'bath_%s_cyc_%i_step_%i.pack.19'%(label2,mccyc,mcstep))

# convert zp_dihed back into original coordinates only for isFlex=0
 if isFlex == 0: 
  zp_dihed=zp_dihed# -di_offset      # these are settings based on .zmat format 

 return wdir_path,eulvar,comvar,zp_dihed

def manage_job_folders(rootname,sgr,label1,label2):

 dir_path=rootname+'_'+label1+'_'+sgr+'/' 
 if not os.path.exists(dir_path): os.mkdir(dir_path)
# for this monte carlo will will make a first branch and then have a differnt folder for each parrelel task 
# test prototype of single cpu job is just folder labeled 'a'
 wdir_path=dir_path+label2+'/' 
 if not os.path.exists(wdir_path): 
  os.mkdir(wdir_path)
 else: 
  rmtree(wdir_path)
  os.mkdir(wdir_path)

# place to store all pack.10 files from the current MC run 
 packfile_path=dir_path+'packfiles/' 
 if not os.path.exists(packfile_path): 
  os.mkdir(packfile_path)
 else: 
  rmtree(packfile_path)
  os.mkdir(packfile_path)

 png_path=dir_path+'png/' 
 if not os.path.exists(png_path): 
  os.mkdir(png_path)
 else: 
  rmtree(png_path)
  os.mkdir(png_path)

 log_path=dir_path+'log/' 
 if not os.path.exists(log_path): 
  os.mkdir(log_path)
 else: 
  rmtree(log_path)
  os.mkdir(log_path)
 
 return log_path,png_path

def do_ntrytest(wdir_path):
 ntrytest=False
 ntry=0
 with open(wdir_path+'out12') as file:
         for line in file:
              if re.search('Stopped - change fix params', line):
               ntry= int(line.split()[-1])
               ntrytest=True
 if os.path.exists(wdir_path+'STOP'):
  ntrytest=True
  ntry=-1

 # check pack19 for errors 
 pack19name=wdir_path+'pack.19'
 with open(pack19name) as file:
   for line in file:
          if re.search('\*', line):
             ntrytest=True
             ntry=-2
             break

 return ntry,ntrytest               


def ptrr_pack12(rootname,label1,sgr,zpr,nmpg,Skiprows,label2,MC_invar,MC_other,Kbias,DenExp,Tfactor,eulvar_hills,comvar_hills,send_end):

  lowEthresh=-108.0
  minE=-1000.0
  lowErep=0


  # check for runing upack when  debuging
  # wdir_path,eulvar,comvar,dihed = run_upack_search(rootname,label1,sgr,zpr,nmpg,label2,MC_invar,MC_other,DenExp,Tfactor,Kbias)    


  while minE < lowEthresh :  # adhoc placed a condition on the minE so that dont get abnormal low since we know the global min    
   if lowErep != 0 :  print  ("detected an abnormally low energy of %.3f in bath %s !! restarting this step !!"%(minE,label2))
   # print ("exiting for debug purpose ! while minE < lowEthresh ")
   # exit()
   try:
    ntrytest=True
    while ntrytest:
     wdir_path,eulvar,comvar,dihed = run_upack_search(rootname,label1,sgr,zpr,nmpg,label2,MC_invar,MC_other,DenExp,Tfactor,Kbias)    
     ntry,ntrytest= do_ntrytest(wdir_path)  # return the number of tries and if it satisfied the test 
     if ntry > 0:
      print ("Number of attempts = %i for these minibatch settings in bath %s. restarting !!"%(ntry,label2))
      print (eulvar,comvar,dihed)
     elif ntry < 0:
      print ("Manual STOP or problem with pack.19 in bath %s. restarting !!"%(label2))
     else:
     #  print ("exiting for debug purpose !")
     #  exit()
      continue
   
#       eulvar2, comvar2 and cfg_dihed is not the test coordinate and is read from ecv.18 - we must be careful that they can be transformed it is only used to update each cyle   
    ave_eng,var_eng,ave_bias,eulvar2,comvar2,cfg_dihed,minE,metaV = read_pack_files(wdir_path,Skiprows,zpr,nmpg,eulvar,comvar,dihed,MC_other,Kbias,eulvar_hills,comvar_hills)  # ?? bug fix  
    # print "minE = %.3f " %minE

# hardwired option to print minE config details to screen as part of each step
    do_print_minE_config=False
    if do_print_minE_config:
       print ("MinE eulers bath %s "%label2 +' '.join(eulvar2.astype('str')))   #
      # fhout_mclogs[bath].write(' '.join(eulvar.astype('str')))   # old way
   
       str_comvar2 = ["%.3f" % c for c in comvar2]   # fix made here - the old version is the one we keep
       str_comvar2 = ' '.join(str_comvar2)
       print("MinE centers bath %s %s " %(label2,str_comvar2))
   
       str_cfg_dihed = ["%.1f" % d for d in cfg_dihed]   # fix made here - the old version is the one we keep
       str_cfg_dihed = ' '.join(str_cfg_dihed)
       print("MinE dihedrals bath %s %s " %(label2,str_cfg_dihed))
   
    eulerIsRand=np.where(MC_invar[0][:]>999)
    comIsRand=np.where(MC_invar[2][:]>999)
    for idx in eulerIsRand: eulvar2[idx]=1000
    for idx in comIsRand: comvar2[idx]=1000

   except:
     print ("\n !!! ERROR - there is a problem with PACK12 for %s and bath %s - check manually !!!\n " %(label1,label2))
     eulvar,comvar,dihed,ave_eng,var_eng=[0,0,0,0,0]
     pass

   lowErep+=1

  #send_end.send([eulvar2,comvar2,dihed,ave_eng,var_eng,ave_bias,minE])  # return minE var
  send_end.send([eulvar,comvar,dihed,ave_eng,var_eng,ave_bias,minE,metaV,eulvar2,comvar2,cfg_dihed])  # return ext var  

def run_PT_baths(rootname,label1,label2s,sgr,Skiprows,nrand,natoms,zpr,nmpg,mc_iter_ramp,nrampcyc,Kb,Tinit,Tfini, 
                          eul_init,eul_shift,com_init,com_shift,Tor_fix,Tor_shift,logdir,pngdir,nbaths,isFlex,reFreq,
                          Kbias,DenExp,metaV_ON,excM_ON,upd_sk):                

  MC_other=np.array([nrand,0,0,isFlex]) 

  # in version 8 com_shift eul_shift and kbias  are now arrays with an entry for each bath

  #Tcyc sort out inital variables and Temperatures for the baths per cycle
  Tcyc=[]
  MC_invar=[]
  fhout_mclogs=[]
  for bath in range(nbaths): 
   Tcyc.append(np.r_[np.linspace(Tinit[bath],Tfini[bath],mc_iter_ramp),np.linspace(Tfini[bath],Tinit[bath],mc_iter_ramp)])
   MC_invar.append(list([eul_init,np.zeros(zpr*nmpg*3),com_init,np.zeros(zpr*nmpg*3),Tor_fix,np.zeros(np.shape(Tor_fix))]))
   fhout_mclog=open(logdir+'bath_%s_Tin%.3f_Tfi%.3f.log'%(label2s[bath],Tinit[bath],Tfini[bath]),'w')
   fhout_mclog.write('# bath:%s ,  Tinit(%.3fK):Tfini(%.3fK) , distribution size %i # \n'%(label2s[bath],Tinit[bath],Tfini[bath],nrand))
   fhout_mclog.write('# cycle, step, T, ave_eng, var_eng, ave_bias, e_diff, acc1, acc2, rej, eulers[zpr*nmpg*3], coms[zpr*nmpg*3], torsions[ntor], minE_test, metaV_test, minE, bias, metaV, com_shift, eul_shift, kspr \n')
   fhout_mclogs.append(fhout_mclog)
  Tcyc=np.vstack(Tcyc)

#  for bath in range(nbaths): 
#   fhout_mclogs[bath].write("this is a test %i\n"%bath)

  # keeping track of seperate bath varables  
  eulvar_old = []
  comvar_old = []
  dihed_old = []
  eng_old = []
  bias_old = []
  minE_old = []
  metaV_old = []
  eulvar2_old = []
  comvar2_old = []
  dihed2_old = []

  # store list for history dependant biasing 
  eulvar_hills = []
  comvar_hills = []
  tor_hills = []


  jobs = []
  pipe_list = []
  # get inital params for unshifted cv's
  print ("geting inital params for unshifted cv's")
  
  for bath in range(nbaths): 
   recv_end, send_end = multiprocessing.Pipe(False)
   Tfactor=Kb*Tinit[bath]
   p = multiprocessing.Process(target=ptrr_pack12, args=(rootname,label1,sgr,zpr,nmpg,Skiprows,label2s[bath],MC_invar[bath],MC_other,Kbias[bath],DenExp,Tfactor,eulvar_hills,comvar_hills,send_end,)) 
   jobs.append(p)
   pipe_list.append(recv_end)
   p.start()

  for job in jobs: # use this to parrelelize jobs   
     job.join()

  result_list = [x.recv() for x in pipe_list]

  # treatment for scaling of shifts - original parameters for each bath 
  escal0=eul_shift
  cscal0=com_shift
  tscal0=Tor_shift

  # declare arrays for shift variables using the scalar variable name
  # eul_shift=np.zeros((nbaths,len(eul_init)))
  # com_shift=np.zeros((nbaths,len(com_init)))

  nel=len(eul_init)
  nco=len(com_init)

  # declare a shift parameter for every molecule in each bath  
  eul_shift=np.array(nel*eul_shift).reshape((nel,-1)).transpose()
  com_shift=np.array(nco*com_shift).reshape((nco,-1)).transpose()
  Tor_shift=np.zeros((nbaths,len(Tor_fix)))  # set to 0 for rigid body 


  for bath in range(nbaths): 
    eulvar=result_list[bath][0]  
    comvar=result_list[bath][1]  
    dihed=result_list[bath][2]  
    ave_eng=result_list[bath][3]  
    var_eng=result_list[bath][4]  
    ave_bias=result_list[bath][5]  
    minE=result_list[bath][6]  
    metaV=result_list[bath][7]  
    eulvar2=result_list[bath][8]  
    comvar2=result_list[bath][9]  
    dihed2=result_list[bath][10]  

    eulvar_old.append(eulvar)
    comvar_old.append(comvar)
    eng_old.append(ave_eng)   
    bias_old.append(ave_bias)   
    dihed_old.append(dihed)
    minE_old.append(minE)   
    metaV_old.append(metaV)   
    eulvar2_old.append(eulvar2)
    comvar2_old.append(comvar2)
    dihed2_old.append(dihed2)

 # have to update shift scaling at this point becasue it is bath dependent. its isotropic atm  llikly benefit from aniso treatment  
 # need to be carful here becasue can overlaod the eul_shift and com_shift variable structure - need to clean this as well
 #   escal,cscal = get_shifts(escal0,cscal0,Kb*Tinit[bath],Kbias[0]) 
    escal,cscal = escal0,cscal0 # bypassing scale function , plan to re-implement at later stage
 
    # just duplicate of earlier code - ineffiecnt but works can change later !!! 
#    eul_shift[bath,:]=(np.zeros((len(eul_init)))+np.round(escal)).astype(float)  # this should be a float -see run_upack_search()
#    com_shift[bath,:]=(np.zeros((len(com_init)))+cscal).astype(float)  
#    Tor_shift[bath,:]=(np.zeros((len(Tor_fix)))+tscal0).astype(float)    # not implemented yet for rigid body  

#    eul_shift[bath,:]=eul_shift[bath,:]+np.round(escal)  # this should be a float -see run_upack_search()
#    com_shift[bath,:]=com_shift[bath,:]+cscal
#    Tor_shift[bath,:]=Tor_shift[bath,:]+tscal0   # not implemented yet for rigid body  

#    print  "initial Euler,COM and Tor shifts on bath %s at %.3f with Tfactor %.3f" %(label2s[bath],Tinit[bath],Kb*Tinit[bath])
#    print  eul_shift*(Kb*Tinit[bath]),com_shift*(Kb*Tinit[bath]),Tor_shift*(Kb*Tinit[bath])
    print(  "initial Euler,COM and Tor shifts on bath %s at %.3f with Tfactor %.3f --- final T scaling eulers rounded to nint " %(label2s[bath],Tinit[bath],Kb*Tinit[bath]))
    print(  eul_shift[bath],com_shift[bath],Tor_shift[bath])

#   in version 7 updating shift array on fly  be carful to duplicate array when assigning to parrelel baths
    MC_invar[bath]=[eulvar,eul_shift[bath],comvar,com_shift[bath],dihed,Tor_shift[bath]]

  # updte hills list with bath variables as lists   
  eulvar_hills.append(np.vstack(eulvar_old).copy()) # trick is need to use the .copy() function to make sure not just storing pointer to current object
  comvar_hills.append(np.vstack(comvar_old).copy())

  print( "End initializing baths")
  
  accept_first_move=True

  print ("Begin Parallel Tempering")
  # contiune SA
  acc1=np.zeros(nbaths)
  acc2=np.zeros(nbaths)
  rej=np.zeros(nbaths)
  np.random.seed()   # seeding RNG 
  tostep=0     
  for cyc in range(nrampcyc):
   MC_other[1]=cyc   
   for step in range(len(Tcyc[0])):
    MC_other[2]=step   
    tostep+=1

    # we should be able to perfrom the rep exchane here by flipping MC_invar
    # Propose a replica exchange REFreq part of the time
    # algoritham assumes replica IDs increase with bath temperature
    maxswaps=1
    if excM_ON: maxswaps=nbaths/2
    if (np.random.rand() < reFreq):
     nswaps=0 # number of swaping trails
     while nswaps < maxswaps:
      # pick a random replica
        rxch0=np.random.random_integers(0,nbaths-1)   #
        if rxch0==0 : 
          rxch1=1               # if it is the lowest temp then only one nieghbor
        elif rxch0==nbaths-1 : 
          rxch1=nbaths-2     # if it is the highest temp then only neighbor
        else : 
         rxch1= rxch0+(np.random.random_integers(0,1)*2-1)   # it can be one of its niegbors 

        # Calculate beta differences and potential energy differences
        deltaBeta = 1.0/(Kb*Tcyc[rxch1,step]) - 1.0/(Kb*Tcyc[rxch0,step])

        if metaV_ON:
         deltaU = (minE_old[rxch1]+bias_old[rxch1]+metaV_old[rxch1])-(minE_old[rxch0]+bias_old[rxch0]+metaV_old[rxch0])   # use_minE + metaV
        else:
         deltaU = (minE_old[rxch1]+bias_old[rxch1])-(minE_old[rxch0]+bias_old[rxch0])   # use_minE + bias
        

        # Follow the Metropolis algorithm
        if (deltaBeta*deltaU > 0.0):
            # Accept - thereforce, exchange the systems
            # Exchange x positions and energies  
            MC_invar[rxch0],MC_invar[rxch1]=MC_invar[rxch1],MC_invar[rxch0]    # use comma assignment to swap items on same list inline 
            eng_old[rxch0],eng_old[rxch1]=eng_old[rxch1],eng_old[rxch0]
            bias_old[rxch0],bias_old[rxch1]=bias_old[rxch1],bias_old[rxch0]    # we must also keep track of the bias component
            minE_old[rxch0],minE_old[rxch1]=minE_old[rxch1],minE_old[rxch0]    # and the minE
            metaV_old[rxch0],metaV_old[rxch1]=metaV_old[rxch1],metaV_old[rxch0]    
        
            print ("acc1 exchange at cyc %i step %i between baths %i(T=%.3f) : %i(T=%.3f) with dE %.3f " %(cyc,step,rxch0,Tcyc[rxch0,step],rxch1,Tcyc[rxch1,step],deltaU) )
            # print "debug metaV for baths T=%.3f and T=%.3f was %.3f and %.3f " %(Tcyc[rxch0,step],Tcyc[rxch1,step],metaV_old[rxch0],metaV_old[rxch1]) 
            
        # if deltaBeta*deltaU is positive, accept with the corresponding probability    
        elif (np.random.rand() < np.exp(deltaBeta*deltaU)):
        # elif (np.random.rand() <= 1.0):                 # test  
            # Accept - thereforce, exchange the systems
            # Exchange x positions and energies  
            MC_invar[rxch0],MC_invar[rxch1]=MC_invar[rxch1],MC_invar[rxch0]    
            eng_old[rxch0],eng_old[rxch1]=eng_old[rxch1],eng_old[rxch0]
            bias_old[rxch0],bias_old[rxch1]=bias_old[rxch1],bias_old[rxch0]    # we must also keep track of the bias component
            minE_old[rxch0],minE_old[rxch1]=minE_old[rxch1],minE_old[rxch0]    # and the minE
            metaV_old[rxch0],metaV_old[rxch1]=metaV_old[rxch1],metaV_old[rxch0]    
        
            print ("acc2 exchange at cyc %i step %i between baths %i(T=%.3f) : %i(T=%.3f) with dE %.3f, dBet %.3f and expProb %.3f  " %(cyc,step,rxch0,Tcyc[rxch0,step],rxch1,Tcyc[rxch1,step],deltaU,deltaBeta,np.exp(deltaBeta*deltaU)) )
           # print "debug metaV for baths T=%.3f and T=%.3f was %.3f and %.3f " %(Tcyc[rxch0,step],Tcyc[rxch1,step],metaV_old[rxch0],metaV_old[rxch1]) 
        # else :
            # print "no exchange at cyc %i step %i between baths %i(T=%.3f) : %i(T=%.3f) with dE %.3f, dBet %.3f and expProb %.3f  " %(cyc,step,rxch0,Tcyc[rxch0,step],rxch1,Tcyc[rxch1,step],deltaU,deltaBeta,np.exp(deltaBeta*deltaU)) 
            # print "debug metaV for baths T=%.3f and T=%.3f was %.3f and %.3f " %(Tcyc[rxch0,step],Tcyc[rxch1,step],metaV_old[rxch0],metaV_old[rxch1]) 

        # Otherwise - each replica remains with its own position and potential energy            
        nswaps+=1
########## finished exchange moves #################
    # Perform  regular Monte Carlo SA moves for each bath

    jobs = []
    pipe_list = []
    for bath in range(nbaths): 
      recv_end, send_end = multiprocessing.Pipe(False)
      Tfactor=Kb*Tcyc[bath,step]
      p = multiprocessing.Process(target=ptrr_pack12, args=(rootname,label1,sgr,zpr,nmpg,Skiprows,label2s[bath],MC_invar[bath],MC_other,Kbias[bath],DenExp,Tfactor,eulvar_hills,comvar_hills,send_end,)) 
      jobs.append(p)
      pipe_list.append(recv_end)
      p.start()
     
    for job in jobs: # use this to parrelelize jobs   
        job.join()
     
    result_list = [x.recv() for x in pipe_list]
 
    for bath in range(nbaths): 
      eulvar=result_list[bath][0] 
      comvar=result_list[bath][1]
      dihed=result_list[bath][2]  
      ave_eng=result_list[bath][3]  
      var_eng=result_list[bath][4]  
      ave_bias=result_list[bath][5]  
      minE=result_list[bath][6]  
      metaV=result_list[bath][7]  
      eulvar2=result_list[bath][8] 
      comvar2=result_list[bath][9]
      dihed2=result_list[bath][10]  

      fhout_mclogs[bath].write('%i ' %(cyc))
      fhout_mclogs[bath].write('%i ' %(step))
      fhout_mclogs[bath].write('%.3f ' %(Tcyc[bath,step]))
      fhout_mclogs[bath].write('%.6f %.6f %.6f ' %(ave_eng, var_eng, ave_bias))

      if metaV_ON:
       eng_diff=(minE+ave_bias+metaV)-(minE_old[bath]+bias_old[bath]+metaV_old[bath])   # uses metaV
      else: 
       eng_diff=(minE+ave_bias)-(minE_old[bath]+bias_old[bath])  # turn off metaV
      # eng_diff=(ave_eng+ave_bias)-(eng_old[bath]+bias_old[bath])     # use aveE
      # eng_diff=(ave_bias)-(bias_old[bath])     # use bias
      if eng_diff < 0.0 :          
       acc1[bath]+=1
       MC_invar[bath]=[eulvar,eul_shift[bath],comvar,com_shift[bath],dihed,Tor_shift[bath]]
       eulvar_old[bath]=eulvar
       comvar_old[bath]=comvar
       dihed_old[bath]=dihed
       eulvar2_old[bath]=eulvar2
       comvar2_old[bath]=comvar2
       dihed2_old[bath]=dihed2
       eng_old[bath]=ave_eng   
       bias_old[bath]=ave_bias   
       minE_old[bath]=minE   
       metaV_old[bath]=metaV   
      elif np.random.rand() <= np.exp(-eng_diff/(Kb*Tcyc[bath,step])): 
       acc2[bath]+=1
       MC_invar[bath]=[eulvar,eul_shift[bath],comvar,com_shift[bath],dihed,Tor_shift[bath]]
       eulvar_old[bath]=eulvar
       comvar_old[bath]=comvar
       dihed_old[bath]=dihed
       eulvar2_old[bath]=eulvar2
       comvar2_old[bath]=comvar2
       dihed2_old[bath]=dihed2
       eng_old[bath]=ave_eng   
       bias_old[bath]=ave_bias   
       minE_old[bath]=minE   
       metaV_old[bath]=metaV   
      elif  accept_first_move and tostep==1 : 
       print ("accepted first shift")
       acc2[bath]+=1
       MC_invar[bath]=[eulvar,eul_shift[bath],comvar,com_shift[bath],dihed,Tor_shift[bath]]
       eulvar_old[bath]=eulvar
       comvar_old[bath]=comvar
       dihed_old[bath]=dihed
       eulvar2_old[bath]=eulvar2
       comvar2_old[bath]=comvar2
       dihed2_old[bath]=dihed2
       eng_old[bath]=ave_eng   
       bias_old[bath]=ave_bias   
       minE_old[bath]=minE   
       metaV_old[bath]=metaV   
      else: 
       rej[bath]+=1
       # v5.101 fix here needed to make sure whatever in the bath is actually transfered properly and then printed to the logs 
       # MC_invar[bath]=[eulvar_old[bath],eul_shift,comvar_old[bath],com_shift,dihed_old[bath],Tor_shift] # if we rejected a move we actaully dont need to update  
       eulvar_old[bath]=MC_invar[bath][0]
       comvar_old[bath]=MC_invar[bath][2]
       dihed_old[bath]=MC_invar[bath][4]
       metaV_old[bath]=metabias(eulvar_hills,comvar_hills,eulvar_old[bath],comvar_old[bath],zpr,nmpg) # still need to update metaV on-the-fly even for old positions 
       
      fhout_mclogs[bath].write('%.3f %.2f %.2f %.2f ' %(eng_diff,float(acc1[bath])/tostep,float(acc2[bath])/tostep,float(rej[bath])/tostep))

      str_eulvar_old = ["%i" % e for e in eulvar_old[bath]]   # fix made here - the old version is the one we keep
      fhout_mclogs[bath].write(' '.join(str_eulvar_old))   # fix made here  now outputs correct variable
      str_comvar = ["%.3f" % c for c in comvar_old[bath]]   # fix made here - the old version is the one we keep
      str_comvar = ' '.join(str_comvar)
      str_dihed = ["%.1f" % d for d in dihed_old[bath]]
      str_dihed = ' '.join(str_dihed)
      fhout_mclogs[bath].write(' %s ' %(str_comvar))
      fhout_mclogs[bath].write(' %s ' %(str_dihed))      # fix this for multiple dihedrals  
      fhout_mclogs[bath].write('%.6f ' %(minE))          # this is for the test position 
      fhout_mclogs[bath].write('%.6f ' %(metaV))         # this is for the test position
      fhout_mclogs[bath].write('%.6f ' %(minE_old[bath]))          # this is for the current state v5.0
      fhout_mclogs[bath].write('%.6f ' %(bias_old[bath]))         # this is for the current state
      fhout_mclogs[bath].write('%.6f ' %(metaV_old[bath]))         # this is for the current state

      fhout_mclogs[bath].write('%.3f ' %(com_shift[bath,0]))         # extra info for updating meta-params
      fhout_mclogs[bath].write('%3i ' %(eul_shift[bath,0]))         # extra info for updating meta-params
      fhout_mclogs[bath].write('%.3f ' %(Kbias[bath]))         # extra info for updating meta-params

      fhout_mclogs[bath].write('\n')



      # we need this for next section on updates 

 
      # the isofunc decides how much the lr is step-wise attenuated  by the Tfactor
      # 
      Tfactor=Kb*Tinit[bath]
      gamma=2.0   #  set higher to converge to 1.0 faster i.e low temp will change at a slower rate per step 
      isofunc=(gamma*Tfactor)/(1.0+gamma*Tfactor)   # this si the same function as langmuir isotherm
      

      # update shifts in each bath based on the aceptance 
      update_shifts = True if int(upd_sk[0])==1 else False  # hardwire switch
      slr=float(upd_sk[2])
      if ((step%2 == 0) and (step != 0) and (update_shifts)) : 
            if float(rej[bath])/tostep > 0.5:
                if eul_shift[bath,0] > np.ceil(escal[bath]*slr*Tfactor) : eul_shift[bath,:] = np.floor(eul_shift[bath,:] * (1.0 - slr*isofunc))
                if com_shift[bath,0] > cscal[bath]*slr*Tfactor : com_shift[bath,:] *= (1.0 - slr*isofunc)
            else:
                eul_shift[bath,:] = np.ceil(eul_shift[bath,:] * (1.0 + slr*isofunc))
                com_shift[bath,:] *=  (1.0 + slr*isofunc)        # only applying lr on com for testing

      #   this way of incremting shift values works but is slow to converge  
      #    the above algo is modified in simialr manner to treatment of kspr see if it will converge faster 

      #     if float(rej[bath])/tostep > 0.5:
      #         if eul_shift[bath,0] > 1: eul_shift[bath,:]=eul_shift[bath,:] - 1.0
      #         if com_shift[bath,0] > cscal[bath]*lr: com_shift[bath,:]=com_shift[bath,:]-cscal[bath]*lr  
      #     else:
      #         eul_shift[bath,:]=eul_shift[bath,:]+1.0
      #         com_shift[bath,:]=com_shift[bath,:]+cscal[bath]*lr        # only applying lr on com for testing

            print("bath:", bath , "rejection ratio: ", float(rej[bath])/tostep, " adjusted shifts :",  eul_shift[bath], com_shift[bath])

      # rescale kspr in each bath based on the aceptance 
      update_kspr = True if int(upd_sk[1])==1 else False  # hardwire switch
      klr=float(upd_sk[3])  # increment by some percent of Tfactor
      if ((step%2 == 1) and (step > 1) and update_kspr) : 
            if float(rej[bath])/tostep > 0.5:
                 Kbias[bath] *= (1 - klr*isofunc) 
            else:
                 Kbias[bath] *= (1 + klr*isofunc) 

            print("bath:", bath , "rejection ratio: ", float(rej[bath])/tostep, " adjusted kspr :",  Kbias[bath])
     


  # updte hills list at the end of every step with bath variables as lists    
    eulvar_hills.append(np.vstack(eulvar_old).copy())    # trick is need to use the .copy() function to make sure not just storing pointer to current object
    comvar_hills.append(np.vstack(comvar_old).copy())

   #perfrom overrelaxation update at the end of every cycle  
   # print "end cycle %i" %cyc
   for bath in range(nbaths):           # reason need extra loop here is becasue this works with list and not numpy array
    MC_invar[bath][0]=eulvar2_old[bath]
    MC_invar[bath][2]=comvar2_old[bath]
    if isFlex: MC_invar[bath][4]=dihed2_old[bath] # only need to over relax the dihedrals for isFlex 

  fhout_mclog.close()






def EVCCPMC(fin):


 fhin = open(fin,'r')
 in_vars = fhin.readlines()
 in_vars= [re.sub(r'("(?:[^"]+|(?<=\\)")*")|#[^\n]*','',line) for line in in_vars] # get rid of any commentry after '#'
 fhin.close()

 print ("reading input varibles from %s" %fin )

 rootname=in_vars[0].split()[0]
 label1=in_vars[0].split()[1]
 sgr=in_vars[0].split()[2]

 print ("rootname:%s label1:%s spacegroup:%s" %(rootname,label1,sgr))  
 
 zpr,nrand,natoms,Skiprows,nmpg=np.array(in_vars[1].split()).astype(int)
 mc_iter_ramp,nrampcyc=np.array(in_vars[2].split()).astype(int)
 Tinit=np.array(in_vars[3].split()).astype(float)
 Tfini=np.array(in_vars[4].split()).astype(float)
 eul_init=np.array(in_vars[5].split()).astype(int)
 com_init=np.array(in_vars[6].split()).astype(float)
 Tor_fix=np.array(in_vars[7].split()).astype(float)
 init_shifts = np.array(in_vars[8].split()).astype(float)  # all shifts to be read in  

 # modified for development version 8 - parameters as lists
 e_sv = init_shifts[:len(Tinit)] 
 c_sv = init_shifts[len(Tinit):-1] 
 t_sv = init_shifts[-1]

 eul_shift=list(e_sv) 
 com_shift=list(c_sv) 
 Tor_shift=float(t_sv) 

# print (type(eul_shift))
# print (type(com_shift))
# print (type(Tor_shift))

 # eul_shift=(np.zeros((len(eul_init)))+e_sv).astype(float)  
 # com_shift=(np.zeros((len(com_init)))+c_sv).astype(float)  
 # Tor_shift=(np.zeros((len(Tor_fix)))+t_sv).astype(float)  

 isFlex =int(in_vars[11].split()[0])  # do flexible search option - will automaticall render dihdral as dummy var 
 reFreq =float(in_vars[12].split()[0])  # replica exchange frequency 
 Kbias =np.array(in_vars[12].split()[1:len(Tinit)+1]).astype(float)  # biasing force consts [one for each bath] 
 DenExp =float(in_vars[12].split()[-1])  # expected density
 metaV_ON =int(in_vars[13].split()[0])  # switch for metaV
 excM_ON  =int(in_vars[13].split()[1])  # switch for multiple exchanges per step
 upd_sk  = in_vars[14].split()  # switchs and lr for updating meta params
 
 print (Kbias)
 print (upd_sk)

 # hardwired  variables 
 Kb=0.0083144621  # bolzman factor in kJ/(mol.K) 
 label2s=['a','b','c','d','e','f','g','h','i','j','k','l','m','n','o','p','q'] # label of parellel chains perfromed in seperate folders 
 nbaths=np.shape(Tinit)[0]

 print ("there should be no greater than %i configs generated from this run" %(nrand*nbaths*mc_iter_ramp*2*nrampcyc))
 print ("there are %i baths" %nbaths)

# refresh working folders for each bath
 for bath in range(nbaths): 
  logdir,pngdir = manage_job_folders(rootname,sgr,label1,label2s[bath])

 extime = time.time()
 print("start  time for workflow execution  --- %s seconds ----" % ((time.time() - extime)))


 run_PT_baths(rootname,label1,label2s,sgr,Skiprows,nrand,natoms,zpr,nmpg,mc_iter_ramp,nrampcyc,Kb,Tinit,Tfini,   
              eul_init,eul_shift,com_init,com_shift,Tor_fix,Tor_shift,logdir,pngdir,nbaths,isFlex,reFreq,
              Kbias,DenExp,metaV_ON,excM_ON,upd_sk)                     

 # below was for degugging  
 #return (rootname,label1,label2s,sgr,Skiprows,nrand,natoms,zpr,nmpg,mc_iter_ramp,nrampcyc,Kb,Tinit,Tfini,   
 #             eul_init,eul_shift,com_init,com_shift,Tor_fix,Tor_shift,logdir,pngdir,nbaths,isFlex,reFreq,
 #             Kbias,DenExp,metaV_ON,excM_ON)                     

 print("end  time for workflow execution  --- %s seconds ----" % ((time.time() - extime)))


def main():
    
 from sys import argv
 script, fin = argv
    
 EVCCPMC(fin)

if __name__ == "__main__":
    main()

