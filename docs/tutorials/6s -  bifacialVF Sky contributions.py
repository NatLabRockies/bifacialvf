#!/usr/bin/env python
# coding: utf-8

# ## 6* - BIFACIALvf Sky contributions
# 
# Runs bifacialvf for a day, and compares to same day without DNI, without DHI, and without DNI no albedo, without DHI no albedo to capture the contributions of Ground Reflected vs other sources of DNI and DHI.
# 

# In[5]:


testfolder = r'TEMP'
Resultsfolder = r'TEMP'
exampleflag = False
debugflag = False


# In[6]:


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import style
import pvlib
import datetime
import pprint
import os


# ## 6.a: Run bifacialvf

# In[8]:


import bifacialvf

# Print bifacialvf Version:
bifacialvf.__version__


# In[9]:


TMYtoread=os.path.join(testfolder,'SRRL_WeatherFile_TMY3_15.csv')
writefiletitle=os.path.join(Resultsfolder,'bifacialVF_Results_SRRLAlbedo_60min.csv')

# Variables
tilt = 10                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
albedo = None               # Calculated in previous section from SRRL data. Value is 0.28 up to 11/18/19o
hub_height=1.5/2            #1.5m / 2m collector width
pitch = 2/0.35/2              # 1 / 0.35 where 0.35 is gcr --- row to row spacing in normalized panel lengths. 
rowType = "interior"        # RowType(first interior last single)
transFactor = 0             # TransmissionFactor(open area fraction)
sensorsy = 12                # sensorsy(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
PVfrontSurface = "glass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "glass"     # PVbackSurface(glass or ARglass)

 # Calculate PV Output Through Various Methods    
calculateBilInterpol = False   # Only works with landscape at the moment.
calculatePVMismatch = False
portraitorlandscape='portrait'   # portrait or landscape
cellsnum = 72
bififactor = 1.0

# Tracking instructions
tracking=True
backtrack=True
limit_angle = 65

# Using YEARLY Albedo:
myTMY3, meta = bifacialvf.bifacialvf.readInputTMY(TMYtoread)


# In[14]:


myTMY3.head()


# In[15]:


myTMY3, meta = bifacialvf.bifacialvf.fixintervalTMY(myTMY3,meta)


# In[16]:


myTMY3.head()


# In[19]:


filterdates = (myTMY3.index >= '2020-02-02 0:0:0 -7') & (myTMY3.index < '2020-02-03 0:0:0 -7')
myTMY3[filterdates]


# In[20]:


writefiletitle=os.path.join(Resultsfolder,'bifacialVF_02_02.csv')

bifacialvf.simulate(myTMY3[filterdates], meta, writefiletitle=writefiletitle, 
         tilt=tilt, sazm=sazm, pitch=pitch, hub_height=hub_height, 
         rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
         albedo=albedo, tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, calculatePVMismatch=calculatePVMismatch,
         cellsnum = cellsnum, bififactor=bififactor,
         calculateBilInterpol=calculateBilInterpol,
         portraitorlandscape=portraitorlandscape)


# In[23]:


noDNI = myTMY3[filterdates].copy()
noDNI['DNI'] = 0


# In[ ]:


noDNI.iloc[40:50]


# In[33]:


noDHI = myTMY3[filterdates].copy()
noDHI['DHI'] = 0


# In[35]:


noDNInoAlb = myTMY3[filterdates].copy()
noDNInoAlb['DNI'] = 0
noDNInoAlb['Alb'] = 0.00000000000001


# In[36]:


noDHInoAlb = myTMY3[filterdates].copy()
noDHInoAlb['DHI'] = 0
noDHInoAlb['Alb'] = 0.00000000000001


# In[ ]:





# In[41]:


cases = [noDNI, noDHI, noDNInoAlb, noDHInoAlb]
titles = ['noDNI', 'noDHI', 'noDNInoAlb', 'noDHInoAlb']

for i in range (0, len(cases)):
    title = 'bifacialVF_02_02_'+titles[i]+'.csv'
    case = cases[i]
    
    writefiletitle=os.path.join(Resultsfolder, title)

    bifacialvf.simulate(case, meta, writefiletitle=title, 
             tilt=tilt, sazm=sazm, pitch=pitch, hub_height=hub_height, 
             rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
             albedo=albedo, tracking=tracking, backtrack=backtrack, 
             limit_angle=limit_angle, calculatePVMismatch=calculatePVMismatch,
             cellsnum = cellsnum, bififactor=bififactor,
             calculateBilInterpol=calculateBilInterpol,
             portraitorlandscape=portraitorlandscape)


# In[65]:


# ANALYSE DATA

titles = ['','_noDNI', '_noDHI', '_noDNInoAlb', '_noDHInoAlb']

compiled_RearIrrad = []

for i in range (0, len(titles)):
    title = 'bifacialVF_02_02'+titles[i]+'.csv'
    writefiletitle=os.path.join(Resultsfolder, title)
    
    data, meta = bifacialvf.loadVFresults(writefiletitle)
    
    filterAzm = (data['sazm'] == 270)
    data['fixed_RowBackGTI_1'] = data['No_1_RowBackGTI'] 
    data['fixed_RowBackGTI_1'][filterAzm] = data['No_12_RowBackGTI'][filterAzm]
    compiled_RearIrrad.append(data['fixed_RowBackGTI_1'].sum())
    


# In[72]:


compiled_RearIrrad


# In[88]:


print("Full simulation rear irradiance", compiled_RearIrrad[0])
print("Dni + DHI modeled separately", compiled_RearIrrad[1]+compiled_RearIrrad[2])
print("% from original", (compiled_RearIrrad[1]+compiled_RearIrrad[2])*100/compiled_RearIrrad[0])


# In[89]:


#labels = 'noDNI', 'noDHI', '_noDNInoAlb', '_noDHInoAlb'
# labels aka DHI     DNI    DHI other sources     DNI other sources
DHIonly = compiled_RearIrrad[1] 
DNIonly = compiled_RearIrrad[2]
DHIotherSources = compiled_RearIrrad[3] 
DNIotherSources = compiled_RearIrrad[4] 
ground_reflected_DHI = DHIonly-DHIotherSources
ground_reflected_DNI = DNIonly - DNIotherSources


# In[ ]:





# In[90]:


labels = 'DHI (other sources)',  'DNI (other sources)', 'Ground Reflected DHI','Ground Reflected DNI'
sizes = [DHIotherSources, DNIotherSources, ground_reflected_DHI, ground_reflected_DNI]
    
fig1, ax1 = plt.subplots()
ax1.pie(sizes, labels=labels, autopct='%1.1f%%', startangle=90)
ax1.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.

plt.show()


# In[ ]:




