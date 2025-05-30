#!/usr/bin/env python2
# -*- coding: utf-8 -*-
from __future__ import division, print_function, absolute_import
import math
import csv
import pvlib
import os
#import sys
#import pytz
import numpy as np
import pandas as pd
from tqdm import tqdm

from bifacialvf.vf import getBackSurfaceIrradiances, getFrontSurfaceIrradiances, getGroundShadeFactors
from bifacialvf.vf import getSkyConfigurationFactors, trackingBFvaluescalculator, rowSpacing
from bifacialvf.sun import  perezComp,  sunIncident, sunrisecorrectedsunposition #, hrSolarPos, solarPos,

#from bifacialvf.readepw import readepw

# Electrical Mismatch Calculation 
from bifacialvf.analysis import analyseVFResultsBilInterpol, analyseVFResultsPVMismatch
#import bifacialvf.analysis as analysis



#def simulate(myTMY3, meta, writefiletitle=None, tilt=0, sazm=180, 
#             clearance_height=None, hub_height = None, 
#             pitch=None, rowType='interior', transFactor=0.01, sensorsy=6, 
#             PVfrontSurface='glass', PVbackSurface='glass', albedo=None,  
#             tracking=False, backtrack=True, limit_angle=45,
#             calculatePVMismatch=False, cellsnum= 72, 
#             portraitorlandscape='landscape', bififactor = 1.0,
#             calculateBilInterpol=False, BilInterpolParams=None,
#             deltastyle='TMY3', agriPV=False):
# IO Files

# DEFAULTS:
writefiletitle=None
tilt=0
sazm=180, 
clearance_height=None
hub_height = None
pitch=None
rowType='interior'
transFactor=0.01
sensorsy=6, 
PVfrontSurface='glass'
PVbackSurface='glass'
albedo=None
tracking=False
backtrack=True
limit_angle=45
calculatePVMismatch=False
cellsnum= 72, 
portraitorlandscape='landscape'
bififactor = 1.0,
calculateBilInterpol=False
BilInterpolParams=None,
deltastyle='TMY3'
agriPV=False
 
TMYtoread=r"C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\data\724010TYA.csv"   # VA Richmond
writefiletitle="data/Output/Test_RICHMOND_1.0.csv"

# Variables
tilt = 10                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
albedo = 0.62               # ground albedo
clearance_height=0.4
pitch = 1.5                   # row to row spacing in normalized panel lengths. 
rowType = "interior"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
sensorsy = 6                # sensorsy(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)

 # Calculate PV Output Through Various Methods    
calculateBilInterpol = True   # Only works with landscape at the moment.
calculatePVMismatch = True
portraitorlandscape='landscape'   # portrait or landscape
cellsnum = 72
bififactor = 1.0

# Tracking instructions
tracking=False
backtrack=True
limit_angle = 60

import bifacialvf
# read input
myTMY3, meta = bifacialvf.readInputTMY(TMYtoread)
deltastyle = 'TMY3'

# In[2]:
# Function



if (clearance_height == None) & (hub_height != None):
    clearance_height = hub_height
    if tracking == False:
        print('Warning: hub_height passed and is being used as ',
              'clearance_height for the fixed_tilt routine.')
elif (clearance_height == None) & (hub_height == None):
    raise Exception('No row distance specified in either D or pitch') 
elif (clearance_height != None) & (hub_height == None): 
    if tracking == True:
        print('Warning: clearance_height passed and is being used as ',
              'hub_height for the tracking routine')
else:
    print('Warning: clearance_height and hub_height passed in. Using ' 
          + ('hub_height' if tracking else 'clearance_height') )
    if tracking == True:
        clearance_height = hub_height

C=clearance_height
heightlabel = 'Clearance_Height'

# In[3]:
if tracking == True:
     axis_tilt = 0  # algorithm only allows for zero north-south tilt with SAT
     #limit_angle = 45  # maximum tracker rotation 
     axis_azimuth=sazm    # axis_azimuth is degrees east of North
     tilt = 0            # start with tracker tilt = 0
     hub_height = C      # Ground clearance at tilt = 0.  C >= 0.5
     stowingangle = 90
     if hub_height < 0.5:
         print('Warning: tracker hub height C < 0.5 may result in ground clearance errors')
     heightlabel = 'Hub_Height'

# In[4]:

D = pitch - math.cos(tilt / 180.0 * math.pi)

if writefiletitle == None:
     writefiletitle = "data/Output/TEST.csv"
 
# In[5]:

noRows, noCols = myTMY3.shape
 
lat = meta['latitude']; lng = meta['longitude']; tz = meta['TZ']
 
try:
     name = meta['Name'] #TMY3
except KeyError:  
     name = meta['city'] #EPW
 
## infer the data frequency in minutes
dataInterval = (myTMY3.index[1]-myTMY3.index[0]).total_seconds()/60
 
 # In[6]:
if not (('azimuth' in myTMY3) and ('zenith' in myTMY3) and ('elevation' in myTMY3)):
     solpos, sunup = sunrisecorrectedsunposition(myTMY3, meta, deltastyle = deltastyle)
     myTMY3['zenith'] = np.radians(solpos['zenith'])
     myTMY3['azimuth'] = np.radians(solpos['azimuth'])
     myTMY3['elevation']=np.radians(solpos['elevation'])
 
 # In[7]:

if tracking == True:        
                 
     if not (('trackingdata_surface_tilt' in myTMY3) and ('trackingdata_surface_azimuth' in myTMY3)):
         gcr=1/pitch  
         trackingdata = pvlib.tracking.singleaxis(np.degrees(myTMY3['zenith']), 
                                                  np.degrees(myTMY3['azimuth']),
                                                  axis_tilt, axis_azimuth, 
                                                  limit_angle, backtrack, gcr)
         
         trackingdata.surface_tilt.fillna(stowingangle, inplace=True)
         myTMY3['trackingdata_surface_tilt'] = trackingdata['surface_tilt']         
         myTMY3['trackingdata_surface_azimuth'] = trackingdata['surface_azimuth']      
     
     [myTMY3['C'], myTMY3['D']] = trackingBFvaluescalculator(myTMY3['trackingdata_surface_tilt'], hub_height, pitch)
else: ## HACK
     myTMY3['C'] = C
     myTMY3['D'] = D
     myTMY3['sazm']  = sazm
     myTMY3['tilt'] = tilt

# In[8]:
 # Check what Albedo to se:
if albedo == None:
     if 'Alb' in myTMY3:
         print("Using albedo from TMY3 file.")
         print("Note that at the moment, no validation check is done",
               "in the albedo data, so we assume it's correct and valid.\n")
         useTMYalbedo = True
     else:
         print("No albedo value set or included in TMY3 file", 
               "(TMY Column name 'Alb (unitless)' expected)",
               "Setting albedo default to 0.2\n ")
         albedo = 0.2
         useTMYalbedo=False
else:
     if 'Alb' in myTMY3:
         print("Albedo value passed, but also present in TMY3 file. ",
               "Using albedo value passed. To use the ones in TMY3 file",
               "re-run simulation with albedo=None\n")
     useTMYalbedo=False


 
 # In[9]:
# REAR TEST
     
if tracking==False:        
     ## Sky configuration factors are the same for all times, only based on geometry and row type
     rearSkyConfigFactors, frontSkyConfigFactors, sky1, sky2, sky3 = getSkyConfigurationFactors(rowType, tilt, C, D)       ## Sky configuration factors are the same for all times, only based on geometry and row type
 
WANTED = rearSkyConfigFactors
WANTED1 = sky1
WANTED2 = sky2
WANTED3 = sky3

# In[10]:
 ## Create WriteFile and write labels at this time
 
 #check that the save directory exists, unless it's in root
savedirectory = os.path.dirname(writefiletitle)
if ( (not os.path.exists(savedirectory)) and (savedirectory != '')):
     os.makedirs(savedirectory)
 
# In[11]:

myTMY3['dni'] = myTMY3.DNI            
myTimestamp = myTMY3.index
hour = myTimestamp.hour
minute = myTimestamp.minute
myTMY3['dhi'] = myTMY3.DHI
Tamb=myTMY3.DryBulb
VWind = myTMY3.Wspd
         
if useTMYalbedo:
    albedo = myTMY3.Alb
                                         
zen = myTMY3['zenith']
azm = myTMY3['azimuth']
elv = myTMY3['elevation']

# In[12]:
# CHANGE THIS                

#[WANTEDRear, frontSkyConfigFactors] = getSkyConfigurationFactors(rowType, tilt, C, D)       ## Sky configuration factors are the same for all times, only based on geometry and row type
import matplotlib.pyplot as plt
plt.plot(WANTED, label='rearconfigfactors')
plt.plot(sky1, label='sky1'); plt.plot(sky2, label='sky2'); plt.plot(sky3, label='sky3')
plt.legend()

# In[13]:

daylighthours = myTMY3[myTMY3['zenith']<0.5*math.pi]

#tilt = myTMY3['trackingdata_surface_tilt']
#sazm = myTMY3['trackingdata_surface_azimuth']
C = myTMY3['C']                        
D = myTMY3['D']
     
# In[14]:        
    
# Reduce TMY3
df = myTMY3[['Date (MM/DD/YYYY)', 'Time (HH:MM)','GHI','DNI','DHI','DryBulb','RHum','Pressure','Wdir','Wspd',
            'zenith','azimuth','elevation','dni','dhi','C','D','sazm','tilt']]

df['albedo'] = 0.3
df['pitch'] = 1.5


# In[15]:
# TEST 2

rearConfigFactors, frontSkyConfigFactors = bifacialvf.vf.getSkyConfigurationFactors2(rowType, df, pitch)  











# In[15]:
#getSkyConfigurationFactors2 REDO
import pandas as pd

groundresolution = 100

rearSkyConfigFactors = []
frontSkyConfigFactors = []

beta = df['tilt']
# Tilt from horizontal of the PV modules/panels, in radians
DTOR = math.pi / 180.0  # HACK
beta = beta * DTOR
# Vertical height of sloped PV panel (in PV panel slope lengths)                    
h = beta.apply(math.sin)
# Horizontal distance from front of panel to rear of panel (in PV panel
# slope lengths)
x1 = beta.apply(math.cos)

# Forced fix for case of C = 0
# FIXME: for some reason the Config Factors go from 1 to 2 and not 0 to 1.
# TODO: investigate why this is happening in the code. #Sil: Did we fix this?        
df.loc[df['C']==0,'C']= 0.0000000001
C = df['C']
D = df['D']

if not C[C<0].empty:    
    LOGGER.error("ERROR: Clearance height is below ground level."
                 " Function GetSkyConfigurationFactors "
                 " will continue but results might be unreliable")
    df.loc[df['C']<0,'C']= 0.0000000001

# Divide the row-to-row spacing into 100 intervals and calculate
# configuration factors
delta = pitch / groundresolution # This was 101... why?
x=np.linspace(delta, pitch-delta, num=groundresolution, endpoint=True)
df_x1 = pd.DataFrame(zip(*[x1 for i in range(groundresolution)]))

# In[16]:
# ATAN investigation
 
clear = h+C
df_clear = pd.DataFrame(zip(*[clear for i in range(groundresolution)]))
df = df_x1.subtract(x)  #x1-x
df = df+pitch*2.0
df_yx = df_clear/df


# In[16]:
# If interior:
    
clear = h+C
df_clear = pd.DataFrame(zip(*[clear for i in range(groundresolution)]))
df = df_x1.subtract(x)  #x1-x
df = df+pitch*2.0
#df_yx = df_clear/df
#df_yx.applymap(math.atan)
angA = np.arctan2(df_clear, df)

#angB = math.atan2(C, (2.0 * pitch - x))
df_y = pd.DataFrame(zip(*[C for i in range(groundresolution)]))
df = 2.0*pitch - x
#df_yx = df_y/df
#angB = df_yx.applymap(math.atan)
angB = np.arctan2(df_y, df)

beta1 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

## check 1 rows away
#angA = math.atan2(h + C, (pitch + x1 - x))
df = df_x1.subtract(x)  #x1-x
df = df+pitch
#df_yx = df_clear/df
#angA = df_yx.applymap(math.atan)
angA = np.arctan2(df_clear, df)


#angB = math.atan2(C, (pitch - x))
df = pitch - x
#df_yx = df_y/df
#angB = df_yx.applymap(math.atan)
angB = np.arctan2(df_y, df)

# beta2 = min(angA, angB)
beta2 = angA.where(angA < angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

#beta3 = max(angA, angB)
beta3 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

# check 0 rows away
#beta4 = math.atan2(h + C, (x1 - x))
#df_yx = df_clear / df_x1.subtract(x)  #x1-x
#beta4 = df_yx.applymap(math.atan)
df_x1_mx = df_x1.subtract(x)  #x1-x
beta4 = np.arctan2(df_clear, df_x1_mx)

#beta5 = math.atan2(C, (-x))
#df_yx = df_y / -x
#beta5 = df_yx.applymap(math.atan) 
beta5 = np.arctan2(df_y, -x)

#beta6 = math.atan2(h + C, (-D - x))
df_D = pd.DataFrame(zip(*[D for i in range(groundresolution)]))
df_D = -df_D-x
df_yx = df_clear / df_D
#beta6 = df_yx.applymap(math.atan2)
beta6 = np.arctan2(df_clear, df_D)


# In[17]:
#SIL CONTINUE FROM HERE
sky1 =0; sky2 =0; sky3 =0
rearSkyConfigFactors=[]
frontSkyConfigFactors=[]

#if (beta2 > beta1):
#   sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2))
mask = (beta2-beta1) > 0
sky1 = 0.5*(np.cos(beta1[mask]) - np.cos(beta2[mask]))

#if (beta4 > beta3):
#    sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4))
mask = (beta4-beta3) > 0
sky2 = 0.5*(np.cos(beta3[mask]) - np.cos(beta4[mask]))

#if (beta6 > beta5):
#    sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6))
mask = (beta6-beta5) > 0
sky3 = 0.5*(np.cos(beta5[mask]) - np.cos(beta6[mask]))

sky1.fillna(0, inplace=True)
sky2.fillna(0, inplace=True)
sky3.fillna(0, inplace=True)

skyAll = sky1 + sky2 + sky3

# Save as arrays of values, same for both to the rear and front
rearSkyConfigFactors.append(skyAll)
frontSkyConfigFactors.append(skyAll)        
# End of if "interior"
# In[17]:
plt.plot(skyAll.iloc[0]), plt.plot(WANTED)
# In[17]:
plt.plot(sky1.iloc[0]), plt.plot(sky2.iloc[0]), plt.plot(sky3.iloc[0])
plt.plot(WANTED1, label='sky1'); plt.plot(WANTED2, label='sky2'); plt.plot(WANTED3, label='sky3')
plt.legend()

# In[17]:
plt.plot(sky1.iloc[0])
plt.plot(WANTED1, label='sky1')
plt.legend()


# In[17]:
plt.plot(sky2.iloc[0])
plt.plot(WANTED2, label='sky1')
plt.legend()

# In[17]:
plt.plot(sky3.iloc[0])
plt.plot(WANTED3, label='sky1')
plt.legend()


# In[PART B]:
rearGroundGHI=[]
frontGroundGHI=[]
pvFrontSH1, pvBackSH1, maxShadow1, rearGroundSH1, frontGroundSH1 = getGroundShadeFactors (rowType, tilt, C[0], D[0], elv[0], azm[0], sazm)



# In[Decomposing]:


rearGroundSH = []
frontGroundSH = []

sazm = myTMY3['sazm'].astype(float)
rtr = D + x1;                # Row-to-row distance (in PV panel slope lengths)

# Divide the row-to-row spacing into 100 intervals for calculating ground shade factors
delta = rtr / (groundresolution); # Shouldn't this be + 1
x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals

Lh = (h / np.tan(elv)) * np.cos(sazm - azm); # Horizontal length of shadow perpindicular to row from top of module to bottom of module
Lhc = ((h + C) / np.tan(elv)) * np.cos(sazm - azm); # Horizontal length of shadow perpindicular to row from top of module to ground level
Lc = (C / math.tan(elv)) * math.cos(sazm - azm); # Horizontal length of shadow perpindicular to row from bottom of module to ground level

ss1 = 0.0; se1 = 0.0; ss2 = 0.0; se2 = 0.0;  # Initialize shading start (s) and end (e) to zeros for two potential shading segments
pvFrontSH = 0.0;
pvBackSH = 0.0;

if (rowType == "interior"):

    if (Lh > D): # Front side of PV module partially shaded, back completely shaded, ground completely shaded
    
        pvFrontSH = (Lh - D) / (Lh + x1);
        pvBackSH = 1.0;
        ss1 = 0.0;      # Ground shaded from 0.0 to rtr
        se1 = rtr;
    
    elif (Lh < -(rtr + x1)):  # Back side of PV module partially shaded, front completely shaded, ground completely shaded
    
        pvFrontSH = 1.0;
        pvBackSH = (Lh + rtr + x1) / (Lh + x1);
        ss1 = 0.0;      # Ground shaded from 0.0 to rtr
        se1 = rtr;
    
    else:        # Ground is partially shaded (I assume)            
    
        if (Lhc >= 0.0):     # Shadow to rear of row, module front unshaded, back shaded
        
            pvFrontSH = 0.0;
            pvBackSH = 1.0;
            Ss = Lc;         # Shadow starts at Lc
            Se = Lhc + x1;   # Shadow ends here
            while (Ss > rtr):
            
                Ss -= rtr;          # Put shadow in correct rtr space if needed
                Se -= rtr;
            
            ss1 = Ss;
            se1 = Se;
            if (se1 > rtr):          # then need to use two shade areas
            
                se1 = rtr;
                ss2 = 0.0;
                se2 = Se - rtr;
                if (se2 > ss1):
                    # This would mean ground completely shaded, does this occur?
                    ss1 = 0.0;      # Ground shaded from 0.0 to rtr
                    se1 = rtr;
                
            
        
        else:                # Shadow to front of row, either front or back might be shaded, depending on tilt and other factors
        
            Ss = 0.0;         # Shadow starts at Lc, initialize
            Se = 0.0;         # Shadow ends here, initialize
            if (Lc < Lhc + x1):
            
                pvFrontSH = 0.0;
                pvBackSH = 1.0;
                Ss = Lc;         # Shadow starts at Lc
                Se = Lhc + x1;   # Shadow ends here
            
            else:
            
                pvFrontSH = 1.0;
                pvBackSH = 0.0;
                Ss = Lhc + x1;  # Shadow starts at Lhc + x1
                Se = Lc;        # Shadow ends here
            
            while (Ss < 0.0):
            
                Ss += rtr;          # Put shadow in correct rtr space if needed
                Se += rtr;
            
            ss1 = Ss;
            se1 = Se;
            if (se1 > rtr):         # then need to use two shade areas
            
                se1 = rtr;
                ss2 = 0.0;
                se2 = Se - rtr;
                if (se2 > ss1):
                    # This would mean ground completely shaded, does this occur?
                    ss1 = 0.0;      # Ground shaded from 0.0 to rtr
                    se1 = rtr;
                
            
        
        # End of if (Lh > D) else branching

    delta = rtr / 100.0;
    x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals            
    #for (i = 0; i <= 99; i++)
    for i in range(0,100):
        
        x += delta;
        #if ((x >= ss1 && x < se1) || (x >= ss2 && x < se2)):
        if ((x >= ss1 and x < se1) or (x >= ss2 and x < se2)):
        
            rearGroundSH.append(1);        # x within a shaded interval, set groundSH to 1 to indicate shaded
            frontGroundSH.append(1);       # same for both front and rear
        
        else:
        
            rearGroundSH.append(0);        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
            frontGroundSH.append(0);       # same for both front and rear
        
        #Console.WriteLine("x = 0,6:0.0000 groundSH = 1", x, groundSH[i]);
    
    # End of if row type == "interior"





























     
# In[17]:
# SECOND STEP: GET GROUND SHADE FACTORS


             # Sum the irradiance components for each of the ground segments, to the front and rear of the front of the PV row
             #double iso_dif = 0.0, circ_dif = 0.0, horiz_dif = 0.0, grd_dif = 0.0, beam = 0.0   # For calling PerezComp to break diffuse into components for zero tilt (horizontal)                           
             ghi, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, zen, 0.0, zen)
             
             
             for k in range (0, 100):
             
                 rearGroundGHI.append(iso_dif * rearSkyConfigFactors[k])       # Add diffuse sky component viewed by ground
                 if (rearGroundSH[k] == 0):
                     rearGroundGHI[k] += beam + circ_dif                    # Add beam and circumsolar component if not shaded
                 else:
                     rearGroundGHI[k] += (beam + circ_dif) * transFactor    # Add beam and circumsolar component transmitted thru module spacing if shaded
     
                 frontGroundGHI.append(iso_dif * frontSkyConfigFactors[k])     # Add diffuse sky component viewed by ground
                 if (frontGroundSH[k] == 0):
                     frontGroundGHI[k] += beam + circ_dif                   # Add beam and circumsolar component if not shaded 
                 else:
                     frontGroundGHI[k] += (beam + circ_dif) * transFactor   # Add beam and circumsolar component transmitted thru module spacing if shaded
             
     
             # b. CALCULATE THE AOI CORRECTED IRRADIANCE ON THE FRONT OF THE PV MODULE, AND IRRADIANCE REFLECTED FROM FRONT OF PV MODULE ***************************
             #double[] frontGTI = new double[sensorsy], frontReflected = new double[sensorsy]
             #double aveGroundGHI = 0.0          # Average GHI on ground under PV array
             aveGroundGHI, frontGTI, frontReflected = getFrontSurfaceIrradiances(rowType, maxShadow, PVfrontSurface, tilt, sazm, dni, dhi, C, D, albedo, zen, azm, sensorsy, pvFrontSH, frontGroundGHI)
 
             #double inc, tiltr, sazmr
             inc, tiltr, sazmr = sunIncident(0, tilt, sazm, 45.0, zen, azm)	    # For calling PerezComp to break diffuse into components for 
             save_inc=inc
             gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen)   # Call to get components for the tilt
             save_gtiAllpc=gtiAllpc
             #sw.Write(strLine)
             #sw.Write(",{0,6:0.00}", hour - 0.5 * dataInterval / 60.0 + minute / 60.0)
             #sw.Write(",{0,6:0.0},{1,5:0.0},{2,5:0.0},{3,5:0.0},{4,4:0.00},{5,6:0.0},{6,6:0.0}",
                 #dni * Math.Cos(zen) + dhi, inc * 180.0 / Math.PI, zen * 180.0 / Math.PI, azm * 180.0 / Math.PI, pvFrontSH, aveGroundGHI, gtiAllpc)
     
             # CALCULATE THE AOI CORRECTED IRRADIANCE ON THE BACK OF THE PV MODULE
             #double[] backGTI = new double[sensorsy]
             backGTI, aveGroundGHI = getBackSurfaceIrradiances(rowType, maxShadow, PVbackSurface, tilt, sazm, dni, dhi, C, D, albedo, zen, azm, sensorsy, pvBackSH, rearGroundGHI, frontGroundGHI, frontReflected, offset=0)
        
             inc, tiltr, sazmr = sunIncident(0, 180.0-tilt, sazm-180.0, 45.0, zen, azm)       # For calling PerezComp to break diffuse into components for 
             gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen)   # Call to get components for the tilt
             
             
             ## Write output
             decHRs = hour - 0.5 * dataInterval / 60.0 + minute / 60.0
             ghi_calc = dni * math.cos(zen) + dhi 
             incd = save_inc * 180.0 / math.pi
             zend = zen * 180.0 / math.pi
             azmd = azm * 180.0 / math.pi
             outputvalues=[myTimestamp, dni, dhi, albedo, decHRs, 
                           ghi_calc, incd, zend, azmd, pvFrontSH, aveGroundGHI, 
                           save_gtiAllpc, pvBackSH, aveGroundGHI, 
                           gtiAllpc, maxShadow, Tamb, VWind]
             frontGTIrow=[]
             backGTIrow=[]
             
             # INVERTING Sensor measurements for tracking when tracker
             # facing the west side.
             # TODO: Modify so it works with axis_azm different of 0 
             #        (sazm = 90 or 270 only)
             if tracking == True:                                   
                 if sazm == 270.0:
                     rangestart = sensorsy-1
                     rangeend = -1
                     steprange = -1
                     rearGroundGHI.reverse()
                 else:
                     rangestart = 0
                     rangeend = sensorsy
                     steprange = 1
             else:
                     rangestart = 0
                     rangeend = sensorsy
                     steprange = 1
                     
             for k in range(rangestart, rangeend, steprange):
                 frontGTIrow.append(frontGTI[k])
                 backGTIrow.append(backGTI[k])      
             outputvalues+=frontGTIrow
             outputvalues+=backGTIrow
             
             
             if tracking==True:
                 outputvalues.append(tilt)
                 outputvalues.append(sazm)
                 outputvalues.append(C)
                 outputvalues.append(D)

             if agriPV:
                 outputvalues.append(str(rearGroundGHI).replace(',', ''))
                 
             sw.writerow(outputvalues)
 
 	# End of daylight if loop 
 
# End of myTMY3 rows of data

 if calculateBilInterpol==True:
     analyseVFResultsBilInterpol(filename=writefiletitle, portraitorlandscape=portraitorlandscape, bififactor=bififactor, writefilename=writefiletitle)

 if calculatePVMismatch==True:
     analyseVFResultsPVMismatch(filename=writefiletitle, portraitorlandscape=portraitorlandscape, bififactor=bififactor, numcells=cellsnum, writefilename=writefiletitle)

  

