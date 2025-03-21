#!/usr/bin/env python2
# -*- coding: utf-8 -*-
        #          This program calculates irradiances on the front and back surfaces of bifacial PV modules.
        #          Key dimensions and nomenclature:
        #          tilt = PV module tilt angle from horizontal, in degrees
        #          sazm = PV module surface azimuth from north, in degrees
        #          1.0 = normalized PV module/panel slant height
        #          C = ground clearance of PV module, in PV module/panel slant heights
        #          D = distance between rows, from rear of module to front of module in next row, in PV module/panel slant heights
        #          h = sin(tilt), vertical PV module dimension, in PV module/panel slant heights
        #          x1 = cos(tilt), horizontal PV module dimension, in PV module/panel slant heights
        #          pitch = x1 + D, row-to-row distance, from front of module to front of module in next row, in PV module/panel slant heights
        #          sensorsy = number of horzontal results, usually corresponding to the rows of cells in a PV module/panel along the slope of the sampled axis.
        #          PVfrontSurface = PV module front surface material type, either "glass" or "ARglass"
        #          PVbackSurface = PV module back surfac ematerial type, either "glass" or "ARglass"
        #        
        #         Program flow consists of:
        #          a. Calculate irradiance distribution on ground
        #          b. Calculate AOI corrected irradiance on front of PV module, and irradiance reflected from front of PV module
        #          c. Calculate irradiance on back of PV module

# ensure python3 compatible division and printing
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
from bifacialvf.loadVFresults import loadVFresults

#from bifacialvf.readepw import readepw

# Electrical Mismatch Calculation 
from bifacialvf.analysis import analyseVFResultsBilInterpol, analyseVFResultsPVMismatch
#import bifacialvf.analysis as analysis

def readInputTMY(TMYtoread):
    '''
    ## Read TMY3 data and start loop ~  
    
    Parameters
    ----------
    TMYtoread: TMY3 .csv weather file, which can be downloaded at http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2005/tmy3/by_state_and_city.html
                   Also .epw weather files, which can be downloaded here: https://energyplus.net/weather and here: http://re.jrc.ec.europa.eu/pvg_tools/en/tools.html#TMY
    
    Returns
    dataframe, meta
        
    '''
    import pandas as pd
    def _tmy_reader(TMYtoread):
        try:
            (myTMY3,meta)=pvlib.iotools.read_tmy3(TMYtoread, map_variables=True)
        except TypeError:
            (myTMY3,meta)=pvlib.iotools.read_tmy3(TMYtoread)
        return(myTMY3,meta)
    
    if TMYtoread is None: # if no file passed in, the readtmy3 graphical file picker will open.
        (myTMY3,meta)=_tmy_reader(TMYtoread)  # , coerce_year=2001     
    elif TMYtoread.lower().endswith('.csv') :  
        (myTMY3,meta)=_tmy_reader(TMYtoread)  # , coerce_year=2001      
    elif TMYtoread.lower().endswith('.epw') : 
        (myTMY3,meta) = pvlib.iotools.read_epw(TMYtoread) # requires pvlib > 0.7.0 #, coerce_year=2001
        # rename different field parameters to match DNI, DHI, DryBulb, Wspd
        #pvlib uses -1hr offset that needs to be un-done. Why did they do this?
        myTMY3.index = myTMY3.index+pd.Timedelta(hours=1)   
    else:
        raise Exception('Incorrect extension for TMYtoread. Either .csv (TMY3) .epw or None')
    
    myTMY3.rename(columns={'dni':'DNI', 'ghi':'GHI',
                           'dhi':'DHI',
                           'temp_air':'DryBulb',
                           'wind_speed':'Wspd',
                           'albedo': 'Alb'}, inplace=True)
        
    return myTMY3, meta

def fixintervalTMY(myTMY3, meta):
    '''
    If data is passed in TMY3 format but has a interval smaller than 1 HR, this 
    function fixes the timestamps from the already imported TMY3 data with 
    readInputTMY. It assume there is a column labeld 'Time (HH:MM)' in myTMY3
    '''
    import pandas as pd
    
    myTMY3['Datetime'] = pd.to_datetime(myTMY3['Date (MM/DD/YYYY)'] + ' ' + myTMY3['Time (HH:MM)'])
    myTMY3 = myTMY3.set_index('Datetime').tz_localize(int(meta['TZ'] * 3600))

    return myTMY3, meta

def getEPW(lat=None, lon=None, GetAll=False, path = None):
    """
    Subroutine to download nearest epw files to latitude and longitude provided,
    into the directory EPWs
    based on github/aahoo.
    
    .. warning::
        verify=false is required to operate within NREL's network.
        to avoid annoying warnings, insecurerequestwarning is disabled
        currently this function is not working within NREL's network.  annoying!
    
    Parameters
    ----------
    lat : decimal 
        Used to find closest EPW file.
    lon : decimal 
        Longitude value to find closest EPW file.
    GetAll : boolean 
        Download all available files. Note that no epw file will be loaded into memory
    
    
    """

    import requests, re
    from requests.packages.urllib3.exceptions import InsecureRequestWarning
    requests.packages.urllib3.disable_warnings(InsecureRequestWarning)
    hdr = {'User-Agent' : "Magic Browser",
           'Accept': 'text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8'
           }

    def _setPath(path):
            """
            setPath - move path and working directory
        
            """
            path = os.path.abspath(path)
        
            print('path = '+ path)
            try:
                os.chdir(path)
            except:
                print("Error on Path passed")
        
            # check for path in the new Radiance directory:
            def _checkPath(path):  # create the file structure if it doesn't exist
                if not os.path.exists(path):
                    os.makedirs(path)
                    print('Making path: '+path)
        
            _checkPath('EPWs')
            
    if path is None:
        _setPath(os.getcwd())
    else:
        _setPath(path)
      
    path_to_save = os.path.join('EPWs') # create a directory and write the name of directory here
    if not os.path.exists(path_to_save):
        os.makedirs(path_to_save)

    def _returnEPWnames():
        ''' return a dataframe with the name, lat, lon, url of available files'''
        r = requests.get('https://github.com/NREL/EnergyPlus/raw/develop/weather/master.geojson', verify=False)
        data = r.json() #metadata for available files
        #download lat/lon and url details for each .epw file into a dataframe
        df = pd.DataFrame({'url':[], 'lat':[], 'lon':[], 'name':[]})
        for location in data['features']:
            match = re.search(r'href=[\'"]?([^\'" >]+)', location['properties']['epw'])
            if match:
                url = match.group(1)
                name = url[url.rfind('/') + 1:]
                lontemp = location['geometry']['coordinates'][0]
                lattemp = location['geometry']['coordinates'][1]
                dftemp = pd.DataFrame({'url':[url], 'lat':[lattemp], 'lon':[lontemp], 'name':[name]})
                #df = df.append(dftemp, ignore_index=True)
                df = pd.concat([df, dftemp], ignore_index=True)
        return df

    def _findClosestEPW(lat, lon, df):
        #locate the record with the nearest lat/lon
        errorvec = np.sqrt(np.square(df.lat - lat) + np.square(df.lon - lon))
        index = errorvec.idxmin()
        url = df['url'][index]
        name = df['name'][index]
        return url, name

    def _downloadEPWfile(url, path_to_save, name):
        r = requests.get(url, verify=False, headers=hdr)
        if r.ok:
            filename = os.path.join(path_to_save, name)
            # py2 and 3 compatible: binary write, encode text first
            with open(filename, 'wb') as f:
                f.write(r.text.encode('ascii', 'ignore'))
            print(' ... OK!')
        else:
            print(' connection error status code: %s' %(r.status_code))
            r.raise_for_status()

    # Get the list of EPW filenames and lat/lon
    df = _returnEPWnames()

    # find the closest EPW file to the given lat/lon
    if (lat is not None) & (lon is not None) & (GetAll is False):
        url, name = _findClosestEPW(lat, lon, df)

        # download the EPW file to the local drive.
        print('Getting weather file: ' + name)
        _downloadEPWfile(url, path_to_save, name)
        #self.epwfile = os.path.join('EPWs', name)
        epwfile = os.path.join('EPWs', name)

    elif GetAll is True:
        if input('Downloading ALL EPW files available. OK? [y/n]') == 'y':
            # get all of the EPW files
            for index, row in df.iterrows():
                print('Getting weather file: ' + row['name'])
                _downloadEPWfile(row['url'], path_to_save, row['name'])
        #self.epwfile = None
        epwfile = None
    else:
        print('Nothing returned. Proper usage: epwfile = getEPW(lat,lon)')
        #self.epwfile = None
        epwfile = None
        
    #return self.epwfile
    return epwfile

def simulate(myTMY3, meta, writefiletitle=None, tilt=0, sazm=180, 
             clearance_height=None, hub_height = None, 
             pitch=None, rowType='interior', transFactor=0.01, sensorsy=6, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=None,  
             tracking=False, backtrack=True, limit_angle=45,
             calculatePVMismatch=False, cellsnum= 72, 
             portraitorlandscape='landscape', bififactor = 1.0,
             calculateBilInterpol=False, BilInterpolParams=None,
             deltastyle='TMY3', agriPV=False):

        '''
      
        Description
        -----------
        Main function to run the bifacialvf routines 
    
        Parameters
        ---------- 
        myTMY3 (pd.DataFrame): A pandas DataaFrame containing for each timestep columns:
            DNI, DHI, it can also have DryBulb, Wspd, zenith, azimuth,
        meta (dict): A dictionary conatining keys: 'latitude', 'longitude', 'TZ', 'Name'
        writefiletitle:  name of output file
        tilt:    tilt angle in degrees.  Not used for tracking
        sazm:    surface azimuth orientation in degrees east of north. For tracking this is the tracker axis orientation
        C:       normalized ground clearance.  For trackers, this is the module height at zero tilt
        pitch:     row-to-row normalized distance.  = 1/GCR
        transFactor:   PV module transmission fraction.  Default 1% (0.01)
        sensorsy:      Number of points along the module chord to return irradiance values.  Default 6 (1-up landscape module)
        limit_angle:     1-axis tracking maximum limits of rotation
        tracking, backtrack:  boolean to enable 1-axis tracking and pvlib backtracking algorithm, respectively
        albedo:     If a value is passed, that value will be used for all the simulations.
                    If None is passed (or albedo argument is not passed), program will search the 
                    TMY file for the "Albe (unitless)" column and use those values

        New Parameters: 
        # Dictionary input example:
        # calculateBilInterpol = {'interpolA':0.005, 'IVArray':None, 'beta_voc_all':None, 'm_all':None, 'bee_all':None}

        
        Returns
        -------
        none
        '''    

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

        D = pitch - math.cos(tilt / 180.0 * math.pi)

        if writefiletitle == None:
            writefiletitle = "data/Output/TEST.csv"
        

        noRows, noCols = myTMY3.shape
        lat = meta['latitude']; lng = meta['longitude']; tz = meta['TZ']
        try:
            name = meta['Name'] #TMY3
        except KeyError:  
            name = meta['city'] #EPW
        
        ## infer the data frequency in minutes
        dataInterval = (myTMY3.index[1]-myTMY3.index[0]).total_seconds()/60
    
        if not (('azimuth' in myTMY3) and ('zenith' in myTMY3) and ('elevation' in myTMY3)):
            solpos, sunup = sunrisecorrectedsunposition(myTMY3, meta, deltastyle = deltastyle)
            myTMY3['zenith'] = np.radians(solpos['zenith'])
            myTMY3['azimuth'] = np.radians(solpos['azimuth'])
            myTMY3['elevation']=np.radians(solpos['elevation'])
        
        
        if tracking == True:        
                        
            if not (('trackingdata_surface_tilt' in myTMY3) and ('trackingdata_surface_azimuth' in myTMY3)):
                gcr=1/pitch  
                trackingdata = pvlib.tracking.singleaxis(np.degrees(myTMY3['zenith']), 
                                                         np.degrees(myTMY3['azimuth']),
                                                         axis_tilt, axis_azimuth, 
                                                         limit_angle, backtrack, gcr)
                
                trackingdata['surface_tilt'] = trackingdata.surface_tilt.fillna(stowingangle)
                myTMY3['trackingdata_surface_tilt'] = trackingdata['surface_tilt']         
                myTMY3['trackingdata_surface_azimuth'] = trackingdata['surface_azimuth']      
            
            [myTMY3['C'], myTMY3['D']] = trackingBFvaluescalculator(myTMY3['trackingdata_surface_tilt'], hub_height, pitch)
                
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

        ## Distance between rows for no shading on Dec 21 at 9 am
        print( " ")
        print( "********* ")
        print( "Running Simulation for TMY3: ")
        print( "Location:  ", name)
        print( "Lat: ", lat, " Long: ", lng, " Tz ", tz)
        print( "Parameters: tilt: ", tilt, "  Sazm: ", sazm, "   ", 
              heightlabel, ": ", C, "  Pitch: ", pitch, "  Row type: ", rowType, 
              "  Albedo: ", albedo)
        print( "Saving into", writefiletitle)
        print( " ")
        print( " ")
        
        DD = rowSpacing(tilt, sazm, lat, lng, tz, 9, 0.0);          ## Distance between rows for no shading on Dec 21 at 9 am
        print( "Distance between rows for no shading on Dec 21 at 9 am solar time = ", DD)
        print( "Actual distance between rows = ", D  )
        print( " ")
    
        if tracking==False:        
            ## Sky configuration factors are the same for all times, only based on geometry and row type
            [rearSkyConfigFactors, frontSkyConfigFactors] = getSkyConfigurationFactors(rowType, tilt, C, D)       ## Sky configuration factors are the same for all times, only based on geometry and row type
    
     
        ## Create WriteFile and write labels at this time
        
        #check that the save directory exists, unless it's in root
        savedirectory = os.path.dirname(writefiletitle)
        if ( (not os.path.exists(savedirectory)) and (savedirectory != '')):
            os.makedirs(savedirectory)
        

                    
        with open (writefiletitle,'w') as csvfile:
            sw = csv.writer(csvfile, delimiter=',', quotechar='|', 
                            quoting=csv.QUOTE_MINIMAL, lineterminator='\n')
            # Write Simulation Parameters (from setup file)
            
            
            if tracking==False and backtrack==True:
                print("Warning: tracking=False, but backtracking=True. ",
                      "Setting backtracking=False because it doesn't make ",
                      "sense to backtrack on fixed tilt systems.")
                backtrack = False
            outputheader=['Latitude(deg)','Longitude(deg)', 'Time Zone','Tilt(deg)', 
                         'PV Azimuth(deg)',heightlabel, 'Pitch', 'RowType(first interior last single)',
                         'TransmissionFactor(open area fraction)','sensorsy(# hor rows in panel)', 
                         'PVfrontSurface(glass or ARglass)', 'PVbackSurface(glass or ARglass)',
                         'Albedo',  'Tracking', 'backtracking']
            outputheadervars=[lat, lng, tz, tilt, sazm, clearance_height, pitch, 
                              rowType, transFactor, sensorsy, PVfrontSurface,
                             PVbackSurface, albedo, tracking, backtrack]
            
                                            
            sw.writerow(outputheader)
            sw.writerow(outputheadervars)
            
            # Write Results names"
            allrowfronts=[]
            allrowbacks=[]
            for k in range(0, sensorsy):
                allrowfronts.append("No_"+str(k+1)+"_RowFrontGTI")
                allrowbacks.append("No_"+str(k+1)+"_RowBackGTI")      
            outputtitles=['date', 'DNI', 'DHI', 
                          'albedo', 'decHRs', 'ghi', 'inc', 'zen', 'azm', 'pvFrontSH', 
                         'aveFrontGroundGHI', 'GTIfrontBroadBand', 'pvBackSH', 
                         'aveBackGroundGHI', 'GTIbackBroadBand', 'maxShadow', 'Tamb', 'VWind']
            outputtitles+=allrowfronts
            outputtitles+=allrowbacks
            if tracking == True:
                print( " ***** IMPORTANT --> THIS SIMULATION Has Tracking Activated")
                print( "Backtracking Option is set to: ", backtrack)
                outputtitles+=['tilt']
                outputtitles+=['sazm']
                outputtitles+=['height']
                outputtitles+=['D']

            if agriPV:
                print("Saving Ground Irradiance Values for AgriPV Analysis. ")
                outputtitles+=['Ground Irradiance Values']
            sw.writerow(outputtitles)
            
            ## Loop through file.  TODO: replace this with for loop.
            rl = 0
            
            
            for rl in tqdm(range(noRows)):
                
                myTimestamp=myTMY3.index[rl]
                hour = myTimestamp.hour
                minute = myTimestamp.minute
                dni = myTMY3.DNI.iloc[rl]#get_value(rl,5,"False")
                dhi = myTMY3.DHI.iloc[rl]#get_value(rl,8,"False")
                if 'DryBulb' in myTMY3: Tamb=myTMY3.DryBulb.iloc[rl]
                else: Tamb=0	            
                if 'Wspd' in myTMY3: VWind = myTMY3.Wspd.iloc[rl]	           
                else: VWind=0
                

                
                if useTMYalbedo:
                    albedo = myTMY3.Alb[rl]
                                                              
                zen = myTMY3['zenith'].iloc[rl]
                azm = myTMY3['azimuth'].iloc[rl]
                elv = myTMY3['elevation'].iloc[rl]
    
                if (zen < 0.5 * math.pi):    # If daylight hours
                
                    # a. CALCULATE THE IRRADIANCE DISTRIBUTION ON THE GROUND 
                    #********************************************************
                    #double[] rearGroundGHI = new double[100], frontGroundGHI = new double[100]
                    # For global horizontal irradiance for each of 100 ground segments, to the rear and front of front of row edge         
                    # Determine where on the ground the direct beam is shaded for a sun elevation and azimuth
                    #int[] rearGroundSH = new int[100], frontGroundSH = new int[100]
                    # Front and rear row-to-row spacing divided into 100 segments, (later becomes 1 if direct beam is shaded, 0 if not shaded)
                    #double pvFrontSH = 0.0, pvBackSH = 0.0, maxShadow    
                    # Initialize fraction of PV module front and back surfaces that are shaded to zero (not shaded), and maximum shadow projected from front of row.
                    
                    # TRACKING ROUTINE CALULATING GETSKYCONFIGURATION FACTORS
                    if tracking == True:                                   
                        tilt = myTMY3['trackingdata_surface_tilt'].iloc[rl]
                        sazm = myTMY3['trackingdata_surface_azimuth'].iloc[rl]
                        C = myTMY3['C'].iloc[rl]                        
                        D = myTMY3['D'].iloc[rl]
                        
                        [rearSkyConfigFactors, frontSkyConfigFactors] = getSkyConfigurationFactors(rowType, tilt, C, D)       ## Sky configuration factors are the same for all times, only based on geometry and row type

                    rearGroundGHI=[]
                    frontGroundGHI=[]
                    pvFrontSH, pvBackSH, maxShadow, rearGroundSH, frontGroundSH = getGroundShadeFactors (rowType, tilt, C, D, elv, azm, sazm)
            
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

     
        print( "Finished")
        
        return

def skycomposition_method(myTMY3, spectral_file_path, meta, writefiletitle=None, tilt=0, sazm=180, 
             clearance_height=None, hub_height = None, 
             pitch=None, rowType='interior', transFactor=0.01, sensorsy=6, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=None,  
             tracking=False, backtrack=True, limit_angle=45,
             calculatePVMismatch=False, cellsnum= 72, 
             portraitorlandscape='landscape', bififactor = 1.0,
             calculateBilInterpol=False, BilInterpolParams=None,
             deltastyle='TMY3', agriPV=False):


        '''
      
        Description
        -----------
        Main function to run the bifacialvf routines 
    
        Parameters
        ---------- 
        tmyFile (string): Baseline weather file used to modify input weather data for sky composition method
            DNI, DHI, it can also have DryBulb, Wspd, zenith, azimuth,
        meta (dict): A dictionary conatining keys: 'latitude', 'longitude', 'TZ', 'Name'
        writefiletitle:  name of output file
        tilt:    tilt angle in degrees.  Not used for tracking
        sazm:    surface azimuth orientation in degrees east of north. For tracking this is the tracker axis orientation
        C:       normalized ground clearance.  For trackers, this is the module height at zero tilt
        pitch:     row-to-row normalized distance.  = 1/GCR
        transFactor:   PV module transmission fraction.  Default 1% (0.01)
        sensorsy:      Number of points along the module chord to return irradiance values.  Default 6 (1-up landscape module)
        limit_angle:     1-axis tracking maximum limits of rotation
        tracking, backtrack:  boolean to enable 1-axis tracking and pvlib backtracking algorithm, respectively
        albedo:     If a value is passed, that value will be used for all the simulations.
                    If None is passed (or albedo argument is not passed), program will search the 
                    TMY file for the "Albe (unitless)" column and use those values

        New Parameters: 
        # Dictionary input example:
        # calculateBilInterpol = {'interpolA':0.005, 'IVArray':None, 'beta_voc_all':None, 'm_all':None, 'bee_all':None}

        
        Returns
        -------
        none
        '''    

        if writefiletitle == None:
            writefiletitle = "data/Output/baseline.csv"

        #Baseline
        write_file_baseline = writefiletitle[:-4] + "_baseline.csv"
        simulate(myTMY3=myTMY3, meta=meta, writefiletitle=write_file_baseline, tilt=tilt, sazm=sazm, 
             clearance_height=clearance_height, hub_height = hub_height, 
             pitch=pitch, rowType='interior', transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=albedo,  
             tracking=tracking, backtrack=backtrack, limit_angle=limit_angle,
             calculatePVMismatch=calculatePVMismatch, cellsnum= cellsnum, 
             portraitorlandscape=portraitorlandscape, bififactor = bififactor,
             calculateBilInterpol=calculateBilInterpol, BilInterpolParams=BilInterpolParams,
             deltastyle=deltastyle, agriPV=agriPV)
        
        (data, metadata) = loadVFresults(write_file_baseline)
        composite_data = data

        #C - DNI & Alb to 0 - Only thing available is DHI from other sources
        myTMY3_DNIALB0 = myTMY3.copy()
        myTMY3_DNIALB0.loc[:, "DNI"] = 0
        myTMY3_DNIALB0.loc[:, "Alb"] = 0
        write_file_DNIALB0 = writefiletitle[:-4] + "_DNIALB0.csv"
        myTMY3_DNIALB0.to_csv('DNIALB0_weather.csv')
        simulate(myTMY3=myTMY3_DNIALB0, meta=meta, writefiletitle=write_file_DNIALB0, tilt=tilt, sazm=sazm, 
             clearance_height=clearance_height, hub_height = hub_height, 
             pitch=pitch, rowType='interior', transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=0,  
             tracking=tracking, backtrack=backtrack, limit_angle=limit_angle,
             calculatePVMismatch=calculatePVMismatch, cellsnum= cellsnum, 
             portraitorlandscape=portraitorlandscape, bififactor = bififactor,
             calculateBilInterpol=calculateBilInterpol, BilInterpolParams=BilInterpolParams,
             deltastyle=deltastyle, agriPV=agriPV)
        (data, metadata) = loadVFresults(write_file_DNIALB0)
        composite_data["Avg_RowBackGTI_DHIDirect"] = 0
        for x in np.arange(1, sensorsy+1, 1):
            #data.loc[:, "No_"+str(x)+"_RowBackGTI"] = data.loc[:, "No_"+str(x)+"_RowBackGTI"] - composite_data.loc[:, "No_"+str(x)+"_RowBackGTI_DNI0"]
            old_string = "No_"+str(x)+"_RowBackGTI" 
            new_string = "No_"+str(x) + "_RowBackGTI_DHIDirect"
            composite_data[new_string] = data[old_string]
            composite_data["Avg_RowBackGTI_DHIDirect"] = composite_data["Avg_RowBackGTI_DHIDirect"] + (composite_data[new_string]/sensorsy)
        #composite_data["Avg_RowBackGTI_DHIDirect"] = composite_data["Avg_RowBackGTI_DHIDirect"]
        #data.rename(dict)
        #composite_data = pd.concat([composite_data, data])
        
        #D - DHI & Alb to 0 - only thing available is DNI from other sources
        myTMY3_DHIALB0 = myTMY3.copy()
        myTMY3_DHIALB0.loc[:, "DHI"] = 0
        myTMY3_DHIALB0.loc[:, "Alb"] = 0
        myTMY3_DHIALB0.to_csv('DHIALB0_weather.csv')
        write_file_DHIALB0 = writefiletitle[:-4] + "_DHIALB0.csv"
        simulate(myTMY3=myTMY3_DHIALB0, meta=meta, writefiletitle=write_file_DHIALB0, tilt=tilt, sazm=sazm, 
             clearance_height=clearance_height, hub_height = hub_height, 
             pitch=pitch, rowType='interior', transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=0,  
             tracking=tracking, backtrack=backtrack, limit_angle=limit_angle,
             calculatePVMismatch=calculatePVMismatch, cellsnum= cellsnum, 
             portraitorlandscape=portraitorlandscape, bififactor = bififactor,
             calculateBilInterpol=calculateBilInterpol, BilInterpolParams=BilInterpolParams,
             deltastyle=deltastyle, agriPV=agriPV)
        
        (data, metadata) = loadVFresults(write_file_DHIALB0)
        composite_data["Avg_RowBackGTI_DNIDirect"] = 0
        for x in np.arange(1, sensorsy+1, 1):
            #data.loc[:, "No_"+str(x)+"_RowBackGTI"] = data.loc[:, "No_"+str(x)+"_RowBackGTI"] - composite_data.loc[:, "No_"+str(x)+"_RowBackGTI_DHI0"]
            old_string = "No_"+str(x)+"_RowBackGTI"
            new_string = "No_"+str(x) + "_RowBackGTI_DNIDirect"
            composite_data[new_string] = data[old_string]
            composite_data["Avg_RowBackGTI_DNIDirect"] = composite_data["Avg_RowBackGTI_DNIDirect"] + (composite_data[new_string]/sensorsy)
        #composite_data["Avg_RowBackGTI_DNIDirect"] = composite_data["Avg_RowBackGTI_DNIDirect"] 
        #data.rename(dict)
        #composite_data = pd.concat([composite_data, data])
        
        #A - DHI to 0 - DNI from ground + other sources, subtract other sources
        myTMY3_DHI0 = myTMY3.copy()
        myTMY3_DHI0.loc[:,'DHI'] = 0
        myTMY3_DHI0.to_csv('DHI0_weather.csv')
        write_file_DHI0 = writefiletitle[:-4] + "_DHI0.csv"
        dict = {}
        simulate(myTMY3=myTMY3_DHI0, meta=meta, writefiletitle=write_file_DHI0, tilt=tilt, sazm=sazm, 
             clearance_height=clearance_height, hub_height = hub_height, 
             pitch=pitch, rowType='interior', transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=albedo,  
             tracking=tracking, backtrack=backtrack, limit_angle=limit_angle,
             calculatePVMismatch=calculatePVMismatch, cellsnum= cellsnum, 
             portraitorlandscape=portraitorlandscape, bififactor = bififactor,
             calculateBilInterpol=calculateBilInterpol, BilInterpolParams=BilInterpolParams,
             deltastyle=deltastyle, agriPV=agriPV)
        (data, metadata) = loadVFresults(write_file_DHI0)
        composite_data["Avg_RowBackGTI_DNIReflected"] = 0
        for x in np.arange(1, sensorsy+1, 1):
            old_string = "No_"+str(x)+"_RowBackGTI"
            new_string = "No_" + str(x) + "_RowBackGTI_DNIReflected"
            composite_data[new_string] = data[old_string] - composite_data["No_"+ str(x) + "_RowBackGTI_DNIDirect"]
            composite_data["Avg_RowBackGTI_DNIReflected"] = composite_data["Avg_RowBackGTI_DNIReflected"] + (composite_data[new_string])/sensorsy

        #data.rename(dict)
        #composite_data = pd.concat([composite_data, data])
        
        #B - DNI to 0, DHI from ground and other sources, subtract other sources
        myTMY3_DNI0 = myTMY3.copy()
        myTMY3_DNI0.loc[:, "DNI"] = 0
        myTMY3_DNI0.to_csv('DNI0_weather.csv')
        write_file_DNI0 = writefiletitle[:-4] + "_DNI0.csv"
        simulate(myTMY3=myTMY3_DNI0, meta=meta, writefiletitle=write_file_DNI0, tilt=tilt, sazm=sazm, 
             clearance_height=clearance_height, hub_height = hub_height, 
             pitch=pitch, rowType='interior', transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface='glass', PVbackSurface='glass', albedo=albedo,  
             tracking=tracking, backtrack=backtrack, limit_angle=limit_angle,
             calculatePVMismatch=calculatePVMismatch, cellsnum= cellsnum, 
             portraitorlandscape=portraitorlandscape, bififactor = bififactor,
             calculateBilInterpol=calculateBilInterpol, BilInterpolParams=BilInterpolParams,
             deltastyle=deltastyle, agriPV=agriPV)
        (data, metadata) = loadVFresults(write_file_DNI0)
        composite_data["Avg_RowBackGTI_DHIReflected"] = 0
        for x in np.arange(1, sensorsy+1, 1):
            old_string = "No_"+str(x)+ "_RowBackGTI"
            new_string = "No_"+ str(x) + "_RowBackGTI_DHIReflected"
            composite_data[new_string] = data[old_string] - composite_data["No_"+ str(x) + "_RowBackGTI_DHIDirect"]
            composite_data["Avg_RowBackGTI_DHIReflected"] = composite_data["Avg_RowBackGTI_DHIReflected"] + (composite_data[new_string])/sensorsy
        #data.rename(dict)
        #composite_data = pd.concat([composite_data, data])

        min_lambda = 280
        max_lambda = 400
        lambda_iter = 5
        lambda_range = np.arange(min_lambda, max_lambda, lambda_iter)
        #spectral_file_path = ''
        composite_data['Sum_DNI'] = 0
        composite_data['Sum_DHI'] = 0
        composite_data['Sum_DNI_ALB'] = 0
        composite_data['Sum_DHI_ALB'] = 0
        for file in os.listdir(spectral_file_path):
            spec_tmy = pd.read_csv('data/spectral_tmys/' + file, skiprows=1)
            composite_data['Sum_DNI'] = composite_data['Sum_DNI'] + spec_tmy['DNI']
            composite_data['Sum_DHI'] = composite_data['Sum_DHI'] + spec_tmy['DHI']
            composite_data['Sum_DNI_ALB'] = composite_data['Sum_DNI_ALB'] + spec_tmy['DNI'] * spec_tmy['ALB']
            composite_data['Sum_DHI_ALB'] = composite_data['Sum_DHI_ALB'] + spec_tmy['DHI'] * spec_tmy['ALB']
            num_from_str = int(file[-8:-4])
            if num_from_str in lambda_range:
                composite_data['DNI_' + file[-8:-4]] = spec_tmy['DNI']
                composite_data['DHI_' + file[-8:-4]] = spec_tmy['DHI']
                composite_data['ALB_' + file[-8:-4]] = spec_tmy['ALB']

        
        for x in lambda_range:
            composite_data['G_rear_DNI_' + str(x)] = composite_data['Avg_RowBackGTI_DNIDirect'] / composite_data['Sum_DNI'] * composite_data['DNI_0' + str(x)]
            composite_data['G_rear_DHI_' + str(x)] = composite_data['Avg_RowBackGTI_DHIDirect'] / composite_data['Sum_DHI'] * composite_data['DHI_0' + str(x)]
            composite_data['G_rear_DNI_reflected_' + str(x)] = composite_data['Avg_RowBackGTI_DNIReflected'] / composite_data['Sum_DNI_ALB'] * composite_data['DNI_0' + str(x)] * composite_data['ALB_0' + str(x)]
            composite_data['G_rear_DHI_reflected_' + str(x)] = composite_data['Avg_RowBackGTI_DHIReflected'] / composite_data['Sum_DHI_ALB'] * composite_data['DHI_0' + str(x)] * composite_data['ALB_0' + str(x)]
            composite_data['G_rear_' + str(x)] = composite_data['G_rear_DNI_' + str(x)] + composite_data['G_rear_DHI_' + str(x)] + composite_data['G_rear_DNI_reflected_' + str(x)] + composite_data['G_rear_DHI_reflected_' + str(x)]

        
        composite_data.to_csv(writefiletitle[:-4] + "_composite.csv")

        return composite_data

        
if __name__ == "__main__":    

    # IO Files
    TMYtoread="data/724010TYA.csv"   # VA Richmond
    writefiletitle="data/Output/Test_RICHMOND_1.0.csv"
    composition_file = "data/Output/Test_sky_composition.csv"

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

    # read input
    myTMY3, meta = readInputTMY(TMYtoread)
    deltastyle = 'TMY3'
    # Function
    simulate(myTMY3, meta, writefiletitle=writefiletitle, 
             tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
             rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
             albedo=albedo, tracking=tracking, backtrack=backtrack, 
             limit_angle=limit_angle, calculatePVMismatch=calculatePVMismatch,
             cellsnum = cellsnum, bififactor=bififactor,
             calculateBilInterpol=calculateBilInterpol,
             portraitorlandscape=portraitorlandscape, deltastyle=deltastyle)
    
    skycomposition_method(myTMY3, meta, writefiletitle=composition_file, 
             tilt=tilt, sazm=sazm, pitch=pitch, clearance_height=clearance_height, 
             rowType=rowType, transFactor=transFactor, sensorsy=sensorsy, 
             PVfrontSurface=PVfrontSurface, PVbackSurface=PVbackSurface, 
             albedo=albedo, tracking=tracking, backtrack=backtrack, 
             limit_angle=limit_angle, calculatePVMismatch=calculatePVMismatch,
             cellsnum = cellsnum, bififactor=bififactor,
             calculateBilInterpol=calculateBilInterpol,
             portraitorlandscape=portraitorlandscape, deltastyle=deltastyle)
                                        
    #Load the results from the resultfile
    from loadVFresults import loadVFresults
    (data, metadata) = loadVFresults(writefiletitle)
    #print data.keys()
    # calculate average front and back global tilted irradiance across the module chord
    data['GTIFrontavg'] = data[['No_1_RowFrontGTI', 'No_2_RowFrontGTI','No_3_RowFrontGTI','No_4_RowFrontGTI','No_5_RowFrontGTI','No_6_RowFrontGTI']].mean(axis=1)
    data['GTIBackavg'] = data[['No_1_RowBackGTI', 'No_2_RowBackGTI','No_3_RowBackGTI','No_4_RowBackGTI','No_5_RowBackGTI','No_6_RowBackGTI']].mean(axis=1)
    
    # Print the annual bifacial ratio.
    frontIrrSum = data['GTIFrontavg'].sum()
    backIrrSum = data['GTIBackavg'].sum()
    print('\n The bifacial ratio for this run is: {:.1f}%'.format(backIrrSum/frontIrrSum*100))
    #print("--- %s seconds ---" % (time.time() - start_time))
