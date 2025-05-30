#!/usr/bin/env python
# coding: utf-8

# Cell 1 (MARKDOWN)
# # Bilinear Interpolation EU PV SEC
# This code uses bifacialvf_Mismatch version, which has bilinear interpolation routines
# 
# Bilinear Interpolatoin which must be included in a subfolder in the bifacial_vf folder.
# bifacialvf > bifacialvf > BF_BifacialIrradiances
# 

# In[1]:


import sys

sys.path.insert(0, '../bifacialvf')

try:
    from bifacialvf import *
    print " Successful import of bifacialvf_mismatch version . "
except ImportError:
    raise RuntimeError('bifacialvf_mismatch is required. download distribution')
    # Simple example system using Radiance.
    


# In[ ]:


# FIXED S"YSTEMS.
# Using random values for the setup based on Richmond VA.
# Looping for SHANGHAI and CAIRO, at different clearance heights. 

beta = 10                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
C = 0.75                      # GroundClearance(panel slope lengths). For tracking this is tilt = 0 hub height 
D = 0.51519                 # DistanceBetweenRows(panel slope lengths)
rowType = "interior"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
cellRows = 6                # CellRows(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)
albedo = 0.62               # ground albedo
calculateDetailedMismatch = False

 #BILINEAR INTERPOLATION VALUES    
calculateBilInterpol = True
calculatePVMismatch = True
# PORTRAIT BILINEAR INTERPOLATION NOT WORKING ATM!!!!! 
portraitorlandscape='landscape'   # portrait or landscape
mat_contents = sio.loadmat(r'BF_BifacialIrradiances\BilinearInterpParams\IVArrayYingli.mat')
IVArray=mat_contents['IVArray']        
mat_contents = sio.loadmat(r'BF_BifacialIrradiances\BilinearInterpParams\newBilinearParamsYingLi.mat')
beta_voc_all=mat_contents['beta_voc_all']        
m_all=mat_contents['m_all']
bee_all=mat_contents['bee_all']
interpolA = 0.005  # More accurate interpolation. Do 0.01 as an option.

#TRACKING VALUES
tracking=False
backtrack=True
rtr=1.5
max_angle=60

'''
#BILINEAR INTERPOLATION VALUES    
calculateBilInterpol = True
calculateDetailedMismatch = True
mat_contents = sio.loadmat(r'C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\BF_BifacialIrradiances\BilinearInterpParams\IVArrayYingli.mat')
IVArray=mat_contents['IVArray']        
mat_contents = sio.loadmat(r'C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\BF_BifacialIrradiances\BilinearInterpParams\newBilinearParamsYingLi.mat')
beta_voc_all=mat_contents['beta_voc_all']        
m_all=mat_contents['m_all']
bee_all=mat_contents['bee_all']
interpolA = 0.005  # More accurate interpolation. Do 0.01 as an option.
portraitorlandscape='landscape'
'''

for i in range (0, 5):
    if i == 0:
        C = 0.15
    if i == 1:
        C = 0.25
    if i == 2:
        C = 0.55
    if i == 3:
        C = 0.75
    if i == 4:
        C = 1

    for j in range (0,2):
        if j == 0:
            TMYtoread=r"C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\data\CHN_Shanghai.Shanghai.583670_IWEC.epw"   # China
            City='SHANGHAI'

        if j ==1:
            TMYtoread=r"C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\data\EGY_Cairo.623660_IWEC.epw"   # China
            City='CAIRO'

        if j ==2:
            TMYtoread=r"C:\Users\sayala\Documents\GitHub\bifacialvf\bifacialvf\data\724010TYA.csv"   # China
            City='Richmond'

        writefiletitle="data/Output/"+City+'_C_'+str(C)+".csv"


        simulate(TMYtoread, writefiletitle, beta, sazm, C, D, 
                 rowType, transFactor, cellRows, 
                 PVfrontSurface, PVbackSurface, albedo, 
                 tracking, backtrack, rtr, max_angle,
                 calculatePVMismatch, portraitorlandscape, 
                 calculateBilInterpol, interpolA, IVArray, beta_voc_all,
                 m_all, bee_all)

print (" ")
print ("FINISHED Finished")


# In[ ]: