!/usr/bin/env python
# coding: utf-8

# # BifacialvF- PVMISMATCH Testing Landscape and Portrait Results.
# 

# In[1]:


import bifacialvf
import pandas as pd


# In[4]:


frontGTIrow=[53.01720526, 55.32772694, 55.78247445, 55.69094834, 55.62417481, 55.59217422]
backGTIrow=[8.550246348, 4.626045218, 4.406054333, 6.141527447, 9.211464402, 12.87567548]

portraitorlandscape='portrait'
sensorsy=6
debug=True
                       
cellCenterPVM, stdpl, cellsx, cellsy = bifacialvf.analysis.setupforPVMismatch(portraitorlandscape=portraitorlandscape, sensorsy=sensorsy)

PowerAveraged, PowerDetailed, sunmatAveraged, sunmatDetailed = bifacialvf.analysis.calculateVFPVMismatch(cellCenterPVM=cellCenterPVM, stdpl=stdpl, cellsx=cellsx, cellsy=cellsy, sensorsy=sensorsy, frontGTIrow=frontGTIrow, backGTIrow=backGTIrow, debug=debug)

print("Results: ", round(PowerAveraged,2), round(PowerDetailed,2))

print(portraitorlandscape, cellsx, cellsy)
print(pd.DataFrame(stdpl))
print(pd.DataFrame(sunmatAveraged).round(3))
pd.DataFrame(sunmatDetailed).round(3)


# In[6]:


frontGTIrow=[53.01720526, 55.32772694, 55.78247445, 55.69094834, 55.62417481, 55.59217422]
backGTIrow=[8.550246348, 4.626045218, 4.406054333, 6.141527447, 9.211464402, 12.87567548]

portraitorlandscape='landscape'
sensorsy=6
                       
cellCenterPVM, stdpl, cellsx, cellsy = bifacialvf.analysis.setupforPVMismatch(portraitorlandscape=portraitorlandscape, sensorsy=sensorsy)

PowerAveraged, PowerDetailed, sunmatAveraged, sunmatDetailed = bifacialvf.analysis.calculateVFPVMismatch(cellCenterPVM, stdpl, cellsx, cellsy, sensorsy, frontGTIrow, backGTIrow, debug=True)

print("Results: ", round(PowerAveraged,2), round(PowerDetailed,2))

print(portraitorlandscape, cellsx, cellsy)
print(pd.DataFrame(stdpl))
print(pd.DataFrame(sunmatAveraged).round(3))
pd.DataFrame(sunmatDetailed).round(3)

