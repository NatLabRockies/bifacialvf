#!/usr/bin/env python
# coding: utf-8

# # Modeling_Sandia_Tracker_Final_2

# In[3]:


import pandas as pd
import scipy.io as sio
import matplotlib.pyplot as plt
import sys
import datetime
import os
sys.path.insert(0, '../bifacialvf')
from modelPaper import readpvdaqSNL_CORRECT, simulatepvDAQ

print "Loaded"


# In[3]:


import os


# In[4]:


TMYtoread=os.path.join(DATADIR, "USA_VA_Richmond.Intl.AP.724010_TMY.epw")   # VA Richmond EPW
writefiletitle=os.path.join(TESTDIR, "_RICHMOND_1axis.csv")

# Variables
tilt = 0                   # PV tilt (deg)
sazm = 180                  # PV Azimuth(deg) or tracker axis direction
albedo = 0.2               # ground albedo
clearance_height = 1
pitch = 2.86                # 0.35 GCR. row to row spacing in normalized panel lengths. 
rowType = "single"        # RowType(first interior last single)
transFactor = 0.013         # TransmissionFactor(open area fraction)
sensorsy = 6                # sensorsy(# hor rows in panel)   <--> THIS ASSUMES LANDSCAPE ORIENTATION 
PVfrontSurface = "ARglass"    # PVfrontSurface(glass or ARglass)
PVbackSurface = "ARglass"     # PVbackSurface(glass or ARglass)

# Tracking instructions
tracking=True
backtrack=True
limit_angle = 60

deltastyle = 'TMY3'
# Calculate PV Output Through Various Methods    


# read input
myTMY3, meta = bifacialvf.bifacialvf.readInputTMY(TMYtoread)
#deltastyle = 'TMY3'
myTMY3_2 = myTMY3.iloc[0:24].copy()
# Simulate just the first 24 hours of the Richmond data file. SINGLE
bifacialvf.simulate(myTMY3_2, meta, writefiletitle=writefiletitle, 
         tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
         rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
         PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
         albedo=albedo, tracking=tracking, backtrack=backtrack, 
         limit_angle=limit_angle, deltastyle=deltastyle)


# In[4]:


pvdaq_SNL_Filetitle = "..\\bifacialvf\data\pvdaq_data_1405_2017.csv"
pvDaq_SNL = readpvdaqSNL_CORRECT(pvdaq_SNL_Filetitle, decimate= False, decimatestyle = "Mean" , decimatevalue = 5 )    

timestamps_all = pvDaq_SNL.index
writefiletitle = "..\\bifacialvf\data\Output\SNL_Tracking_output_Jupyter.csv"

name = "[1429] RTC,SNL,Baseline"
lat = 35.0549
lng = -106.5433
elevation = 1658
transmittance = 0.75
dataInterval=1
offset=0
tzstr='US/Mountain'
Cv = 0.5
C = Cv
rtr = 2
tracking = True
backtrack = True
albedo = 0.6
tz = -7
rowType = 'single'
transFactor = 0.016
cellRows = 6
beta = 0
sazm = 180
PVfrontSurface = 'glass'
PVbackSurface = 'glass'
gcr = 2.0/7.0

simulatepvDAQ(pvDaq_SNL, writefiletitle,  beta= beta, sazm= sazm, C= Cv, D= None,
         rowType= rowType, transFactor= transFactor, cellRows= cellRows, 
         PVfrontSurface= PVfrontSurface, PVbackSurface= PVbackSurface, albedo= albedo,  
         tracking= tracking, backtrack = backtrack, rtr= rtr, Cv= Cv, offset= offset, 
         lat = lat, lng = lng, tz = tz, name = name, dataInterval= dataInterval)


mat_contents = sio.loadmat('data\SNL_TrackingBifi_SAT_Data.mat')

DateTime_E=mat_contents['DateTime_E']
DateTime_W=mat_contents['DateTime_W']
EETS_Eback=mat_contents['EETS_Eback']
EETS_Efront=mat_contents['EETS_Efront']
EETS_Wback=mat_contents['EETS_Wback']
EETS_Wfront=mat_contents['EETS_Wfront']
TrackerAngle_E=mat_contents['TrackerAngle_E']
TrackerAngle_W=mat_contents['TrackerAngle_W']

import datetime
timestamp_E = []
EETS_Eback_new=[]
EETS_Efront_new=[]
TrackerAngle_E_new=[]
timestamp_W = []
EETS_Wback_new=[]
EETS_Wfront_new=[]
TrackerAngle_W_new=[]


for i in range (0, len(DateTime_E)):
    A=str(DateTime_E[i])
    datehour, minute, second = A.split(':')
    yeartrash, month, dayhour = datehour.split('-')
    day, hour = dayhour.split(' ')
    year = 2017; month = int(month); day = int(day); hour = int(hour); minute = int(minute); second=0
    temptime=datetime.datetime(year, month, day, hour, minute, 0, 0)                
    timestamp_E.append(temptime)
    EETS_Eback_new.append(float(EETS_Eback[i]))
    EETS_Efront_new.append(float(EETS_Efront[i]))
    TrackerAngle_E_new.append(float(TrackerAngle_E[i]))


for i in range (0, len(DateTime_W)):
    A=str(DateTime_W[i])
    datehour, minute, second = A.split(':')
    yeartrash, month, dayhour = datehour.split('-')
    day, hour = dayhour.split(' ')
    year = 2017; month = int(month); day = int(day); hour = int(hour); minute = int(minute); second=0
    temptime=datetime.datetime(year, month, day, hour, minute, 0, 0)                
    timestamp_W.append(temptime)
    EETS_Wback_new.append(float(EETS_Wback[i]))
    EETS_Wfront_new.append(float(EETS_Wfront[i]))
    TrackerAngle_W_new.append(float(TrackerAngle_W[i]))

dataE=pd.DataFrame({'EETS_Eback': EETS_Eback_new, 'EETS_Efront': EETS_Efront_new, 'TrackerAngle_E': TrackerAngle_E_new})
dataW=pd.DataFrame({'EETS_Wback': EETS_Wback_new, 'EETS_Wfront': EETS_Wfront_new, 'TrackerAngle_W': TrackerAngle_W_new})
dataE = pd.DataFrame.from_records(dataE, index=timestamp_E)
dataW = pd.DataFrame.from_records(dataW, index=timestamp_W)

#Load the results from the resultfile
from loadVFresults import loadVFresults
#print data.keys()
(data, metadata) = loadVFresults(writefiletitle)
# calculate average front and back global tilted irradiance across the module chord
data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)
# Print the annual bifacial ratio.
frontIrrSum = data['GTIFrontavg'].sum()
backIrrSum = data['GTIBackavg'].sum()

# Choosing to compare datapoints where both E and W have data, by joining the tables and droping NAs
timestamps_model=[]
GTIFrontavg_new=[]
GTIBackavg_new=[]
for i in range (0, len(data)):
    timestamps_model.append(datetime.datetime(data['Year'][i], data['Month'][i], data['Day'][i], data['Hour'][i], data['Minute'][i], 0, 0))
    GTIFrontavg_new.append(float(data['GTIFrontavg'][i]))
    GTIBackavg_new.append(float(data['GTIBackavg'][i]))

dataModel = pd.DataFrame({'GTIFrontavg': GTIFrontavg_new, 'GTIBackavg': GTIBackavg_new})
dataModel = pd.DataFrame.from_records(dataModel, index=timestamps_model)

resultBoth = dataW.join(dataE, how="outer").join(dataModel, how="outer")
resultBoth=resultBoth.dropna()      

resultBoth['EETS_ERatio']=resultBoth['EETS_Eback']/resultBoth['EETS_Efront']
resultBoth['EETS_WRatio']=resultBoth['EETS_Wback']/resultBoth['EETS_Wfront']
resultBoth['GTIModelRatio']=resultBoth['GTIBackavg']/resultBoth['GTIFrontavg']
GainModel = resultBoth['GTIBackavg'].sum()/resultBoth['GTIFrontavg'].sum()*100
GainE = resultBoth['EETS_Eback'].sum()/resultBoth['EETS_Efront'].sum()*100
GainW = resultBoth['EETS_Wback'].sum()/resultBoth['EETS_Wfront'].sum()*100

plt.figure(figsize=(20,10))
plt.plot(resultBoth['GTIModelRatio'], label='Model, gain '+ str(round(GainModel,2))); plt.plot(resultBoth['EETS_ERatio'], label='EETS_E, gain '+ str(round(GainE,2))); plt.plot(resultBoth['EETS_WRatio'], label='EETS_W, gain '+ str(round(GainW,2))); plt.legend()
plt.show()
print('EAST TRACKER: The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(C,rtr,GainW))
print('WEST TRACKER: The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(C,rtr,GainE))
print('MODEL TRACKER: The bifacial ratio for ground clearance {} and rtr spacing {} is: {:.1f}%'.format(C,rtr,GainModel))


# In[ ]:


