"""
Tests of the bifacialvf module. 
 
to run coverage test, py.test --cov-report term-missing --cov=bifacialvf
"""
import pytest
import numpy as np
import pandas as pd
import os
import bifacialvf
from bifacialvf.tests import (
    FIXED_ENDTOEND_GTIFRONT, FIXED_ENDTOEND_GTIBACK,
    TRACKED_ENDTOEND_GTIFRONT, TRACKED_ENDTOEND_GTIBACK)

FIXED_ENDTOEND_GTIFRONT = [8.4, 55.2, 71.1, 108.8, 124.7, 139.3,
111.9, 115.2, 58.6, 19.5]

TESTDIR = os.path.dirname(__file__)  # this folder
DATADIR = os.path.abspath(os.path.join(TESTDIR,'..','data'))
MARION_C_File = os.path.join(TESTDIR,'724010TYA_BillMarion_C_Results.csv')

print(TESTDIR)
print(DATADIR)
#os.path.abspath(os.path.join(os.path.dirname( __file__ ), '..', 'templates'))
def test_readInputTMY():
    ''' 
    read in 724010TYA.csv and USA_VA_Richmond.Intl.AP.724010_TMY.epw
    6/22 GHI is 909, sunup, sunrise at 7:00 and 18:00 on 12/31
    '''
    (myTMY3, meta) = bifacialvf.bifacialvf.readInputTMY(os.path.join(DATADIR,"724010TYA.CSV"))
    assert int(myTMY3.loc['1985-06-22 12:0:0'].GHI) == 909
    assert np.allclose(myTMY3.loc[pd.to_datetime(['1978-12-31 7:0:0-5:00','1978-12-31 18:0:0-5:00']),'GHI'].values, np.array([0,0]))
    
    (myEPW, metaEPW) = bifacialvf.bifacialvf.readInputTMY(os.path.join(DATADIR,"USA_VA_Richmond.Intl.AP.724010_TMY.epw"))
    assert int(myEPW.loc['1963-06-22 12:0:0'].GHI) == 858
    assert np.allclose(myTMY3.loc[pd.to_datetime(['1978-12-31 7:0:0-5:00','1978-12-31 18:0:0-5:00']),'GHI'].values, np.array([0,0]))
    
def test_endtoend():
    '''
    end to end test, first 24 hours of VA Richmond
    
    '''
    #TODO:  consolidate and improve this
    
    # IO Files
    TMYtoread=os.path.join(DATADIR,"724010TYA.CSV")   # VA Richmond 724010TYA.csv
    writefiletitle=os.path.join(DATADIR,"_RICHMOND_1.0.csv")
    
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
    #calculateBilInterpol = True   # Only works with landscape at the moment.
    #calculatePVMismatch = True
    #portraitorlandscape='landscape'   # portrait or landscape
    #cellsnum = 72
    #bififactor = 1.0
    
    # Tracking instructions
    tracking=False
    backtrack=False
    limit_angle = 60

    # read input
    myTMY3, meta = bifacialvf.bifacialvf.readInputTMY(TMYtoread)
    deltastyle = 'TMY3'
    myTMY3_2 = myTMY3.iloc[0:24].copy()
    # Simulate just the first 24 hours of the Richmond data file
    bifacialvf.simulate(myTMY3_2, meta, writefiletitle=writefiletitle, 
             tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
             rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
             albedo=albedo, tracking=tracking, backtrack=backtrack, 
             limit_angle=limit_angle, deltastyle=deltastyle)
                                        
    #Load the results from the resultfile
    from bifacialvf import loadVFresults
    (data, metadata) = loadVFresults(writefiletitle)
    data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
    data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)

    # Lad comparison data    
    (data_Marion, metadata2) = loadVFresults(MARION_C_File)
    # calculate average front and back global tilted irradiance across the module chord
    data_Marion['GTIFrontavg'] = data_Marion[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
    data_Marion['GTIBackavg'] = data_Marion[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)

    assert np.allclose(list(round(data['GTIFrontavg'][1:4],1)), list(round(data_Marion['GTIFrontavg'][1:4],1)))
    assert np.allclose(list(round(data['GTIFrontavg'][1:4],1)), FIXED_ENDTOEND_GTIFRONT[1:4])
    assert np.allclose(list(round(data['GTIBackavg'][1:4],1)), list(round(data_Marion['GTIBackavg'][1:4],1)))
    assert np.allclose(list(round(data['GTIbackBroadBand'][1:10],1)), list(round(data_Marion['GTIbackBroadBand'][1:10],1)))    


def test_1axis_endtoend():
    '''
    end to end test, first 24 hours of VA Richmond .EPW file
    this one uses ARGlass and a single row.
    
    '''
    #TODO:  consolidate and improve this
    
    # IO Files
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
                                        
    #Load the results from the resultfile
    from bifacialvf import loadVFresults
    (data, metadata) = loadVFresults(writefiletitle)
    #print data.keys()
    # calculate average front and back global tilted irradiance across the module chord
    data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
    data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)
    assert np.allclose(data['GTIFrontavg'].array, TRACKED_ENDTOEND_GTIFRONT)
    assert np.allclose(data['GTIBackavg'].array, TRACKED_ENDTOEND_GTIBACK)

def test_bilininterpol():
    import pandas as pd
    inputfile = os.path.join(TESTDIR, "Test_RICHMOND_mismatch.csv")
    outputfile = os.path.join(TESTDIR, "Test_RICHMOND_mismatch_BilInterpol.csv")
    portraitorlandscape='landscape'   # portrait or landscape
    bififactor = 1.0
    
    bifacialvf.bifacialvf.analyseVFResultsBilInterpol(filename=inputfile, portraitorlandscape=portraitorlandscape, bififactor=bififactor)
    results = pd.read_csv(outputfile, skiprows=2)
    assert results['BilInterpol FRONT + BACK (Averaged) PmaxIdeal [W]'].sum() == \
        pytest.approx(247.595, abs = 0.001)
    assert results['BilInterpol FRONT + BACK (Detailed) PmaxUnmatched [W]'].sum() == \
        pytest.approx(247.165, abs = 0.001)


def test_pvmismatch():
    import pandas as pd
    inputfile = os.path.join(TESTDIR, "Test_RICHMOND_mismatch.csv")
    outputfile = os.path.join(TESTDIR, "Test_RICHMOND_mismatch_PVMismatch.csv")
    portraitorlandscape='landscape'   # portrait or landscape
    bififactor = 1.0
    numcells = 72

    bifacialvf.bifacialvf.analyseVFResultsPVMismatch(filename=inputfile, portraitorlandscape=portraitorlandscape, bififactor=bififactor, numcells=numcells)
    results = pd.read_csv(outputfile, skiprows=2)
    assert results['PVMismatch FRONT + BACK (Averaged) PmaxIdeal [W]'].sum() == \
        pytest.approx(195.309, abs = 0.001)
    assert results['PVMismatch FRONT + BACK (Detailed) PmaxUnmatched [W]'].sum() == \
        pytest.approx(194.879, abs = 0.001)    
    

'''  FROM test_vf with nice test fixtures n stuff
@pytest.mark.parametrize('beta, C, D, expected',
    [(160, 0.5, 1, SKY_BETA160_C05_D1), (20, 0.5, 1, SKY_BETA20_C05_D1),
     (20, 0, 1, SKY_BETA20_C0_D1), (160, 0, 1, SKY_BETA160_C0_D1),
     (160, 1, 1, SKY_BETA160_C1_D1), (20, 1, 1, SKY_BETA20_C1_D1),
     (20, 1, 0, SKY_BETA20_C1_D0), (160, 1, 0, SKY_BETA160_C1_D0),
     (160, 0.5, 0, SKY_BETA160_C05_D0), (20, 0.5, 0, SKY_BETA20_C05_D0)])
def test_getSkyConfigurationFactors(beta, C, D, expected):
    """
    Benchmark against to the master branch on 2018-08-20 at 91e785d.
    """
    assert np.allclose(
        getSkyConfigurationFactors("interior", beta, C, D), expected)
'''