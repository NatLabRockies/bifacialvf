# -*- coding: utf-8 -*-
"""
ViewFactor module - VF calculation helper files for bifacial-viewfactor
    
@author Bill Marion
@translated to python by sayala 06/09/17

"""

# ensure python3 compatible division and printing
from __future__ import division, print_function, absolute_import
import math
import numpy as np
from bifacialvf.sun import solarPos, sunIncident, perezComp, aOIcorrection
import logging

# TODO: set level or add formatters if more advanced logging required
LOGGER = logging.getLogger(__name__)  # only used to raise errors
DTOR = math.pi / 180.0  # Factor for converting from degrees to radians


def getBackSurfaceIrradiances(rowType, maxShadow, PVbackSurface, beta, sazm,
                              dni, dhi, C, D, albedo, zen, azm, cellRows,
                              pvBackSH, rearGroundGHI, frontGroundGHI,
                              frontReflected, offset=0):      
    """
    This method calculates the AOI corrected irradiance on the back of the PV
    module/panel. 11/19/2015

    Added rowType and other changes to distinguish between types of rows.
    4/19/2016

    Added input of offset of reference cell from PV module back (in PV panel
    slope lengths) for modeling Sara's reference cell measurements, should be
    set to zero for PV module cell irradiances.

    Added while loop so projected Xs aren't too negative causing array index
    problems (<0) 12/13/2016::

        while (projectedX1 < -100.0 || projectedX2 < -100.0):
            # Offset so array indexes are >= -100.0  12/13/2016

            projectedX1 += 100.0;
            projectedX2 += 100.0;

    Parameters
    ----------
    rowType : str
        Type of row: "first", "interior", "last", or "single" 
    maxShadow
        Maximum shadow length projected to the front(-) or rear (+) from the
        front of the module 
    PVbackSurface
        PV module back surface material type, either "glass" or "ARglass"
    beta
        Tilt from horizontal of the PV modules/panels (deg) (for front surface)
    sazm
        Surface azimuth of PV panels (deg) (for front surface)
    dni
        Direct normal irradiance (W/m2)
    dhi
        Diffuse horizontal irradiance (W/m2)
    C
        Ground clearance of PV panel (in PV panel slope lengths)            
    D
        Horizontal distance between rows of PV panels (in PV panel slope
        lengths) 
    albedo
        Ground albedo
    zen
        Sun zenith (in radians)
    azm
        Sun azimuth (in radians)
    pvBackSH
        Decimal fraction of the back surface of the PV panel that is shaded,
        0.0 to 1.0
    rearGroundGHI : array of size [100]
        Global horizontal irradiance for each of 100 ground segments (W/m2)
    frontGroundGHI : array of size [100]
        Global horizontal irradiance for each of 100 ground segments (W/m2)
    frontReflected : array of size [cellRows]
        Irradiance reflected from the front of the PV module/panel (W/m2) in
        the row behind the one of interest
    offset
        Offset of reference cell from PV module back (in PV panel slope
        lengths), set to zero for PV module cell irradiances

    Returns
    -------
    backGTI : array of size [cellRows]
        AOI corrected irradiance on back side of PV module/panel, one for each
        cell row (W/m2)
    aveGroundGHI : numeric
        Average GHI on ground under PV array

    Notes
    -----
    1-degree hemispherical segment AOI correction factor for glass (index=0)
    and ARglass (index=1)
    """
    backGTI = []

    SegAOIcor = [
        [0.057563, 0.128570, 0.199651, 0.265024, 0.324661, 0.378968, 0.428391, 0.473670, 0.514788, 0.552454, 
         0.586857, 0.618484, 0.647076, 0.673762, 0.698029, 0.720118, 0.740726, 0.759671, 0.776946, 0.792833, 
         0.807374, 0.821010, 0.833534, 0.845241, 0.855524, 0.865562, 0.874567, 0.882831, 0.890769, 0.897939, 
         0.904373, 0.910646, 0.916297, 0.921589, 0.926512, 0.930906, 0.935179, 0.939074, 0.942627, 0.946009, 
         0.949096, 0.952030, 0.954555, 0.957157, 0.959669, 0.961500, 0.963481, 0.965353, 0.967387, 0.968580, 
         0.970311, 0.971567, 0.972948, 0.974114, 0.975264, 0.976287, 0.977213, 0.978142, 0.979057, 0.979662, 
         0.980460, 0.981100, 0.981771, 0.982459, 0.982837, 0.983199, 0.983956, 0.984156, 0.984682, 0.985026, 
         0.985364, 0.985645, 0.985954, 0.986241, 0.986484, 0.986686, 0.986895, 0.987043, 0.987287, 0.987388, 
         0.987541, 0.987669, 0.987755, 0.987877, 0.987903, 0.987996, 0.988022, 0.988091, 0.988104, 0.988114, 
         0.988114, 0.988104, 0.988091, 0.988022, 0.987996, 0.987903, 0.987877, 0.987755, 0.987669, 0.987541, 
         0.987388, 0.987287, 0.987043, 0.986895, 0.986686, 0.986484, 0.986240, 0.985954, 0.985645, 0.985364, 
         0.985020, 0.984676, 0.984156, 0.983956, 0.983199, 0.982837, 0.982459, 0.981771, 0.981100, 0.980460, 
         0.979662, 0.979057, 0.978142, 0.977213, 0.976287, 0.975264, 0.974114, 0.972947, 0.971567, 0.970311, 
         0.968580, 0.967387, 0.965353, 0.963481, 0.961501, 0.959671, 0.957157, 0.954555, 0.952030, 0.949096, 
         0.946009, 0.942627, 0.939074, 0.935179, 0.930906, 0.926512, 0.921589, 0.916297, 0.910646, 0.904373, 
         0.897939, 0.890769, 0.882831, 0.874567, 0.865562, 0.855524, 0.845241, 0.833534, 0.821010, 0.807374, 
         0.792833, 0.776946, 0.759671, 0.740726, 0.720118, 0.698029, 0.673762, 0.647076, 0.618484, 0.586857, 
         0.552454, 0.514788, 0.473670, 0.428391, 0.378968, 0.324661, 0.265024, 0.199651, 0.128570, 0.057563],
        [0.062742, 0.139913, 0.216842, 0.287226, 0.351055, 0.408796, 0.460966, 0.508397, 0.551116, 0.589915,
         0.625035, 0.657029, 0.685667, 0.712150, 0.735991, 0.757467, 0.777313, 0.795374, 0.811669, 0.826496, 
         0.839932, 0.852416, 0.863766, 0.874277, 0.883399, 0.892242, 0.900084, 0.907216, 0.914023, 0.920103, 
         0.925504, 0.930744, 0.935424, 0.939752, 0.943788, 0.947313, 0.950768, 0.953860, 0.956675, 0.959339, 
         0.961755, 0.964039, 0.965984, 0.967994, 0.969968, 0.971283, 0.972800, 0.974223, 0.975784, 0.976647, 
         0.977953, 0.978887, 0.979922, 0.980773, 0.981637, 0.982386, 0.983068, 0.983759, 0.984436, 0.984855, 
         0.985453, 0.985916, 0.986417, 0.986934, 0.987182, 0.987435, 0.988022, 0.988146, 0.988537, 0.988792, 
         0.989043, 0.989235, 0.989470, 0.989681, 0.989857, 0.990006, 0.990159, 0.990263, 0.990455, 0.990515, 
         0.990636, 0.990731, 0.990787, 0.990884, 0.990900, 0.990971, 0.990986, 0.991042, 0.991048, 0.991057, 
         0.991057, 0.991048, 0.991042, 0.990986, 0.990971, 0.990900, 0.990884, 0.990787, 0.990731, 0.990636, 
         0.990515, 0.990455, 0.990263, 0.990159, 0.990006, 0.989857, 0.989681, 0.989470, 0.989235, 0.989043, 
         0.988787, 0.988532, 0.988146, 0.988022, 0.987435, 0.987182, 0.986934, 0.986417, 0.985916, 0.985453, 
         0.984855, 0.984436, 0.983759, 0.983068, 0.982386, 0.981637, 0.980773, 0.979920, 0.978887, 0.977953, 
         0.976647, 0.975784, 0.974223, 0.972800, 0.971284, 0.969970, 0.967994, 0.965984, 0.964039, 0.961755, 
         0.959339, 0.956675, 0.953860, 0.950768, 0.947313, 0.943788, 0.939752, 0.935424, 0.930744, 0.925504, 
         0.920103, 0.914023, 0.907216, 0.900084, 0.892242, 0.883399, 0.874277, 0.863766, 0.852416, 0.839932, 
         0.826496, 0.811669, 0.795374, 0.777313, 0.757467, 0.735991, 0.712150, 0.685667, 0.657029, 0.625035, 
         0.589915, 0.551116, 0.508397, 0.460966, 0.408796, 0.351055, 0.287226, 0.216842, 0.139913, 0.062742]]

    # Tilt from horizontal of the PV modules/panels, in radians
    beta = beta * DTOR
    sazm = sazm * DTOR  # Surface azimuth of PV module/panels, in radians

    # 1. Calculate and assign various paramters to be used for modeling
    #    irradiances

    # For calling PerezComp to break diffuse into components for zero tilt
    # (horizontal)
    iso_dif = 0.0; circ_dif = 0.0; horiz_dif = 0.0; grd_dif = 0.0; beam = 0.0
    # Call to get iso_dif for horizontal surface
    ghi, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(
        dni, dhi, albedo, zen, 0.0, zen)

    # Isotropic irradiance from sky on horizontal surface, used later for
    # determining isotropic sky component
    iso_sky_dif = iso_dif

    # For calling PerezComp to break diffuse into components for 90 degree tilt
    # (vertical)
    inc, tiltr, sazmr = sunIncident(0, 90.0, 180.0, 45.0, zen, azm)

    # Call to get horiz_dif for vertical surface
    vti, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(
        dni, dhi, albedo, inc, tiltr, zen)

    # Horizon diffuse irradiance on a vertical surface, used later for
    # determining horizon brightening irradiance component
    F2DHI = horiz_dif

    index = -99
    n2 = -99.9
    if (PVbackSurface == "glass"):
        # Index to use with 1-degree hemispherical segment AOI correction
        # factor array
        index = 0
        n2 = 1.526  # Index of refraction for glass
    elif (PVbackSurface == "ARglass"):
        # Index to use with 1-degree hemispherical segment AOI correction
        # factor array
        index = 1
        n2 = 1.300  # Index of refraction for ARglass
    else:
        raise Exception(
            "Incorrect text input for PVbackSurface."
            " Must be glass or ARglass.")

    # Reflectance at normal incidence, Duffie and Beckman p217
    Ro = math.pow((n2 - 1.0) / (n2 + 1.0), 2.0)

    # Average GHI on ground under PV array for cases when x projection exceed
    # 2*rtr
    aveGroundGHI = 0.0
    for i in range(0,100):
        aveGroundGHI += rearGroundGHI[i] / 100.0

    # Calculate x,y coordinates of bottom and top edges of PV row in back of desired PV row so that portions of sky and ground viewed by the 
    # PV cell may be determined. Origin of x-y axis is the ground pobelow the lower front edge of the desired PV row. The row in back of 
    # the desired row is in the positive x direction.

    h = math.sin(beta);          # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    x1 = math.cos(beta);         # Horizontal distance from front of panel to rear of panel (in PV panel slope lengths)
    rtr = D + x1;                # Row-to-row distance (in PV panel slope lengths)
    PbotX = rtr;                 # x value for poon bottom egde of PV module/panel of row in back of (in PV panel slope lengths)
    PbotY = C;                   # y value for poon bottom egde of PV module/panel of row in back of (in PV panel slope lengths)
    PtopX = rtr + x1;            # x value for poon top egde of PV module/panel of row in back of (in PV panel slope lengths)
    PtopY = h + C;               # y value for poon top egde of PV module/panel of row in back of (in PV panel slope lengths)

    # 2. Calculate diffuse and direct component irradiances for each cell row 
    for i in range (0, cellRows):
    
        # Calculate diffuse irradiances and reflected amounts for each cell row over it's field of view of 180 degrees, 
        # beginning with the angle providing the upper most view of the sky (j=0)
        #PcellX = x1 * (i + 0.5) / ((double)cellRows);                    # x value for location of PV cell
        #PcellY = C + h * (i + 0.5) / ((double)cellRows);                 # y value for location of PV cell
        PcellX = x1 * (i + 0.5) / (cellRows) + offset * math.sin(beta);    # x value for location of PV cell with OFFSET FOR SARA REFERENCE CELLS     4/26/2016
        PcellY = C + h * (i + 0.5) / (cellRows) - offset * math.cos(beta); # y value for location of PV cell with OFFSET FOR SARA REFERENCE CELLS     4/26/2016
        elvUP = math.atan((PtopY - PcellY) / (PtopX - PcellX));          # Elevation angle up from PV cell to top of PV module/panel, radians
        elvDOWN = math.atan((PcellY - PbotY) / (PbotX - PcellX));        # Elevation angle down from PV cell to bottom of PV module/panel, radians
        if (rowType == "last" or rowType == "single"):                         # 4/19/16 No array to the rear for these cases
        
            elvUP = 0.0;
            elvDOWN = 0.0;
        
        #Console.WriteLine("ElvUp = 0", elvUP / DTOR);
        #if (i == 0)
        #    Console.WriteLine("ElvDown = 0", elvDOWN / DTOR);

        #123
        #iStopIso = Convert.ToInt32((beta - elvUP) / DTOR);        # Last whole degree in arc range that sees sky, first is 0
        #Console.WriteLine("iStopIso = 0", iStopIso);
        #iHorBright = Convert.ToInt32(max(0.0, 6.0 - elvUP / DTOR));    # Number of whole degrees for which horizon brightening occurs
        #iStartGrd = Convert.ToInt32((beta + elvDOWN) / DTOR);               # First whole degree in arc range that sees ground, last is 180

        iStopIso = int(round((beta - elvUP) / DTOR));        # Last whole degree in arc range that sees sky, first is 0
        #Console.WriteLine("iStopIso = 0", iStopIso);
        iHorBright = int(round(max(0.0, 6.0 - elvUP / DTOR)));    # Number of whole degrees for which horizon brightening occurs
        iStartGrd = int(round((beta + elvDOWN) / DTOR));               # First whole degree in arc range that sees ground, last is 180

        backGTI.append(0.0)                                                      # Initialtize front GTI

        for j in range (0, iStopIso):                                      # Add sky diffuse component and horizon brightening if present
        
            backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * iso_sky_dif;                               # Sky radiation
    #          backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * iso_sky_dif;                               # Sky radiation

            if ((iStopIso - j) <= iHorBright):                                   # Add horizon brightening term if seen
            
                backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * F2DHI / 0.052264;  # 0.052246 = 0.5 * [cos(84) - cos(90)]
            #backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * F2DHI / 0.052264;  # 0.052246 = 0.5 * [cos(84) - cos(90)]
            
        

        if (rowType == "interior" or rowType == "first"):                          # 4/19/16 Only add reflections from PV modules for these cases
        

            for j in range (iStopIso, iStartGrd):      #j = iStopIso; j < iStartGrd; j++)                              # Add relections from PV module front surfaces
            
                L = (PbotX - PcellX) / math.cos(elvDOWN);                    # Diagonal distance from cell to bottom of module in row behind
                startAlpha = -(j - iStopIso) * DTOR + elvUP + elvDOWN;
                stopAlpha = -(j + 1 - iStopIso) * DTOR + elvUP + elvDOWN;
                m = L * math.sin(startAlpha);
                theta = math.pi - elvDOWN - (math.pi / 2.0 - startAlpha) - beta;
                projectedX2 = m / math.cos(theta);                           # Projected distance on sloped PV module
                m = L * math.sin(stopAlpha);
                theta = math.pi - elvDOWN - (math.pi / 2.0 - stopAlpha) - beta;
                projectedX1 = m / math.cos(theta);                           # Projected distance on sloped PV module
                projectedX1 = max(0.0, projectedX1);
                #Console.WriteLine("j= 0 projected X1 = 1,6:0.000  projected X2 = 2,6:0.000", j, projectedX1, projectedX2);

                PVreflectedIrr = 0.0;                                        # Irradiance from PV module front cover reflections
                deltaCell = 1.0 / cellRows;                          # Length of cell in sloped direction in module/panel units (dimensionless)
                for k in range (0, cellRows):                                  # Determine which cells in behind row are seen, and their reflected irradiance
                
                    cellBot = k * deltaCell;                                 # Position of bottom of cell along PV module/panel
                    cellTop = (k + 1) * deltaCell;                           # Position of top of cell along PV module/panel
                    cellLengthSeen = 0.0;                                    # Length of cell seen for this row, start with zero
                    if (cellBot >= projectedX1 and cellTop <= projectedX2):
                        cellLengthSeen = cellTop - cellBot;                         # Sees the whole cell
                    elif (cellBot <= projectedX1 and cellTop >= projectedX2):
                        cellLengthSeen = projectedX2 - projectedX1;                 # Sees portion in the middle of cell
                    elif (cellBot >= projectedX1 and projectedX2 > cellBot and cellTop >= projectedX2):
                        cellLengthSeen = projectedX2 - cellBot;                     # Sees bottom of cell
                    elif (cellBot <= projectedX1 and projectedX1 < cellTop and cellTop <= projectedX2):
                        cellLengthSeen = cellTop - projectedX1;                     # Sees top of cell
                    #Console.WriteLine("cell= 0 cellBot = 1,5:0.00 cellTop = 2,5:0.00  Cell length seen = 3,5:0.00", k, cellBot, cellTop, cellLengthSeen);
                    PVreflectedIrr += cellLengthSeen * frontReflected[k];           # Add reflected radiation for this PV cell, if seen, weight by cell length seen
                
                PVreflectedIrr /= projectedX2 - projectedX1;                        # Reflected irradiance from PV modules (W/m2)
                backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * PVreflectedIrr;     # Radiation reflected from PV module surfaces onto back surface of module
            
            # End of adding reflections from PV module surfaces
        #Console.WriteLine("");
        #if (i == 0)
        #Console.WriteLine("iStartGrd = 0", iStartGrd);
        for j in range (iStartGrd, 180):                                  # Add ground reflected component
        
            startElvDown = (j - iStartGrd) * DTOR + elvDOWN;             # Start and ending down elevations for this j loop 
            stopElvDown = (j + 1 - iStartGrd) * DTOR + elvDOWN;
            projectedX2 = PcellX + np.float64(PcellY) / math.tan(startElvDown);      # Projection of ElvDown to ground in +x direction (X1 and X2 opposite nomenclature for front irradiance method)
            projectedX1 = PcellX + PcellY / math.tan(stopElvDown);
            actualGroundGHI = 0.0;                                       # Actuall ground GHI from summing array values
            #if (i == 0)
            #    Console.WriteLine("j= 0 projected X1 = 1,6:0.0", j, 100 * projectedX1 / rtr);
            if (abs(projectedX1 - projectedX2) > 0.99 * rtr):
            
                if (rowType == "last" or rowType == "single"):                  # 4/19/16 No array to rear for these cases
                
                    actualGroundGHI = ghi;                                      # Use total value if projection approximates the rtr
                
                else:
                    actualGroundGHI = aveGroundGHI;                                 # Use average value if projection approximates the rtr                        
            
            else:
            
                projectedX1 = 100.0 * projectedX1 / rtr;                        # Normalize projections and multiply by 100
                projectedX2 = 100.0 * projectedX2 / rtr;
                #Console.WriteLine("projectedX1 = 0 projectedX2 = 1", projectedX1, projectedX2);

                if ((rowType == "last" or rowType == "single") and (abs(projectedX1) > 99.0 or abs(projectedX2) > 99.0)):    #4/19/2016
                
                    actualGroundGHI = ghi;                                      # Use total value if projection > rtr for "last" or "single"
                
                else:
                
                    while (projectedX1 >= 100.0 or projectedX2 >= 100.0):            # Offset so array indexes are less than 100
                    
                        projectedX1 -= 100.0;
                        projectedX2 -= 100.0;
                    
                    while (projectedX1 < -100.0 or projectedX2 < -100.0):            # Offset so array indexes are >= -100.0  12/13/2016
                    
                        projectedX1 += 100.0;
                        projectedX2 += 100.0;
                    

                    #Console.WriteLine("projectedX1 = 0 projectedX2 = 1", projectedX1, projectedX2);
                    index1 = (int)(projectedX1 + 100.0) - 100;                  # Determine indexes for use with rearGroundGHI array and frontGroundGHI array(truncates values)
                    index2 = (int)(projectedX2 + 100.0) - 100;                  # (int)(1.9) = 1 and (int)(-1.9) = -1; (int)(1.9+100) - 100 = 1 and (int)(-1.9+100) - 100 = -2
                    #Console.WriteLine("index1=0 index2=1", index1, index2);
                    if (index1 == index2):
                    
                        if (index1 < 0):
                            actualGroundGHI = frontGroundGHI[index1 + 100];
                        #actualGroundGHI = 0.0;
                        else:
                            actualGroundGHI = rearGroundGHI[index1];                        # x projections in same groundGHI element THIS SEEMS TO ADD HICCUP 4/26/2016 ***************************
                        #actualGroundGHI = 0.0;
                    
                    else:
                    
                        for k in range (index1, index2+1):  #for (k = index1; k <= index2; k++)                      # Sum the irradiances on the ground if projections are in different groundGHI elements
                        
                            if (k == index1):
                            
                                if (k < 0):
                                    actualGroundGHI += frontGroundGHI[k + 100] * (k + 1.0 - projectedX1);
                                else:
                                    actualGroundGHI += rearGroundGHI[k] * (k + 1.0 - projectedX1);
                            
                            elif (k == index2):
                            
                                if (k < 0):
                                    actualGroundGHI += frontGroundGHI[k + 100] * (projectedX2 - k);
                                else:
                                    actualGroundGHI += rearGroundGHI[k] * (projectedX2 - k);
                            
                            else:
                            
                                if (k < 0):
                                    actualGroundGHI += frontGroundGHI[k + 100];
                                else:
                                    actualGroundGHI += rearGroundGHI[k];
                            
                        
                        actualGroundGHI /= projectedX2 - projectedX1;                # Irradiance on ground in the 1 degree field of view
                    
                

                #if (i == 0)
                #    Console.WriteLine("j=0 index1=1 index2=2 projectX1=3,5:0.0 projectX2=4,5:0.0 actualGrdGHI=5,6:0.0", j, index1, index2, projectedX1, projectedX2, actualGroundGHI);
                # End of if looping to determine actualGroundGHI

            
            backGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * actualGroundGHI * albedo;     # Add ground reflected component

            #Console.WriteLine("actualGroundGHI = 0,6:0.0 inputGHI = 1,6:0.0 aveArrayGroundGHI = 2,6:0.0", actualGroundGHI, dhi + dni * math.cos(zen), aveGroundGHI);
            
            # End of j loop for adding ground reflected componenet 

        # Calculate and add direct and circumsolar irradiance components
        inc, tiltr, sazmr = sunIncident(0, 180-beta / DTOR, sazm / DTOR - 180, 45.0, zen, azm)  # For calling PerezComp to break diffuse into components for downward facing tilt
        
        gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen) # Call to get components for the tilt

        cellShade = pvBackSH * cellRows - i;
        if (cellShade > 1.0):    # Fully shaded if > 1, no shade if < 0, otherwise fractionally shaded
            cellShade = 1.0;
        elif (cellShade < 0.0):
            cellShade = 0.0;

        if (cellShade < 1.0 and inc < math.pi / 2.0):  # Cell not shaded entirely and inc < 90 deg
        
            cor = aOIcorrection(n2, inc);                # Get AOI correction for beam and circumsolar
            backGTI[i] += (1.0 - cellShade) * (beam + circ_dif) * cor; # Add beam and circumsolar radiation
        
        # End of for i = 0; i < cellRows loop
    
    return backGTI, aveGroundGHI;
    # End of GetBackSurfaceIrradiances


def getFrontSurfaceIrradiances(rowType, maxShadow, PVfrontSurface, beta, sazm,
                               dni, dhi, C, D, albedo, zen, azm, cellRows,
                               pvFrontSH, frontGroundGHI):      
    """
    This method calculates the AOI corrected irradiance on the front of the PV
    module/panel and the irradiance reflected from the the front of the PV
    module/panel. 11/12/2015
    
    Added row type and MaxShadow and changed code to accommodate 4/19/2015
    
    Parameters
    ----------
    rowType : str
        Type of row: "first", "interior", "last", or "single" 
    maxShadow
        Maximum shadow length projected to the front (-) or rear (+) from the
        front of the module row (in PV panel slope lengths), only used for
        `rowTypes` other than "interior"
    PVfrontSurface
        PV module front surface material type, either "glass" or "ARglass"
    beta
        Tilt from horizontal of the PV modules/panels (deg)
    sazm
        Surface azimuth of PV panels (deg)
    dni
        Direct normal irradiance (W/m2)
    dhi
        Diffuse horizontal irradiance (W/m2)
    C
        Ground clearance of PV panel (in PV panel slope lengths)            
    D
        Horizontal distance between rows of PV panels (in PV panel slope
        lengths) 
    albedo
        Ground albedo
    zen
        Sun zenith (in radians)
    azm
        Sun azimuth (in radians)
    pvFrontSH
        Decimal fraction of the front surface of the PV panel that is shaded,
        0.0 to 1.0
    froutGroundGHI : array of size [100]
        Global horizontal irradiance for each of 100 ground segments in front
        of the module row
    
    Returns
    -------
    frontGTI : array of size [cellRows]
        AOI corrected irradiance on front side of PV module/panel, one for each
        cell row (W/m2)
    frontReflected : array of size [cellRows]
        Irradiance reflected from the front of the PV module/panel (W/m2)
    aveGroundGHI : numeric
        Average GHI on the ground (includes effects of shading by array) from
        the array frontGroundGHI[100]

    Notes
    -----
    1-degree hemispherical segment AOI correction factor for glass (index=0)
    and ARglass (index=1). Creates a list containing 5 lists, each of 8 items,
    all set to 0
    """
    frontGTI = []
    frontReflected = []
    #w, h = 2, 180;
    #SegAOIcor = [[0 for x in range(w)] for y in range(h)] 
    SegAOIcor = ([[0.057563, 0.128570, 0.199651, 0.265024, 0.324661, 0.378968, 0.428391, 0.473670, 0.514788, 0.552454, 
            0.586857, 0.618484, 0.647076, 0.673762, 0.698029, 0.720118, 0.740726, 0.759671, 0.776946, 0.792833, 
            0.807374, 0.821010, 0.833534, 0.845241, 0.855524, 0.865562, 0.874567, 0.882831, 0.890769, 0.897939, 
            0.904373, 0.910646, 0.916297, 0.921589, 0.926512, 0.930906, 0.935179, 0.939074, 0.942627, 0.946009, 
            0.949096, 0.952030, 0.954555, 0.957157, 0.959669, 0.961500, 0.963481, 0.965353, 0.967387, 0.968580, 
            0.970311, 0.971567, 0.972948, 0.974114, 0.975264, 0.976287, 0.977213, 0.978142, 0.979057, 0.979662, 
            0.980460, 0.981100, 0.981771, 0.982459, 0.982837, 0.983199, 0.983956, 0.984156, 0.984682, 0.985026, 
            0.985364, 0.985645, 0.985954, 0.986241, 0.986484, 0.986686, 0.986895, 0.987043, 0.987287, 0.987388, 
            0.987541, 0.987669, 0.987755, 0.987877, 0.987903, 0.987996, 0.988022, 0.988091, 0.988104, 0.988114, 
            0.988114, 0.988104, 0.988091, 0.988022, 0.987996, 0.987903, 0.987877, 0.987755, 0.987669, 0.987541, 
            0.987388, 0.987287, 0.987043, 0.986895, 0.986686, 0.986484, 0.986240, 0.985954, 0.985645, 0.985364, 
            0.985020, 0.984676, 0.984156, 0.983956, 0.983199, 0.982837, 0.982459, 0.981771, 0.981100, 0.980460, 
            0.979662, 0.979057, 0.978142, 0.977213, 0.976287, 0.975264, 0.974114, 0.972947, 0.971567, 0.970311, 
            0.968580, 0.967387, 0.965353, 0.963481, 0.961501, 0.959671, 0.957157, 0.954555, 0.952030, 0.949096, 
            0.946009, 0.942627, 0.939074, 0.935179, 0.930906, 0.926512, 0.921589, 0.916297, 0.910646, 0.904373, 
            0.897939, 0.890769, 0.882831, 0.874567, 0.865562, 0.855524, 0.845241, 0.833534, 0.821010, 0.807374, 
            0.792833, 0.776946, 0.759671, 0.740726, 0.720118, 0.698029, 0.673762, 0.647076, 0.618484, 0.586857, 
            0.552454, 0.514788, 0.473670, 0.428391, 0.378968, 0.324661, 0.265024, 0.199651, 0.128570, 0.057563],
            [0.062742, 0.139913, 0.216842, 0.287226, 0.351055, 0.408796, 0.460966, 0.508397, 0.551116, 0.589915,
            0.625035, 0.657029, 0.685667, 0.712150, 0.735991, 0.757467, 0.777313, 0.795374, 0.811669, 0.826496, 
            0.839932, 0.852416, 0.863766, 0.874277, 0.883399, 0.892242, 0.900084, 0.907216, 0.914023, 0.920103, 
            0.925504, 0.930744, 0.935424, 0.939752, 0.943788, 0.947313, 0.950768, 0.953860, 0.956675, 0.959339, 
            0.961755, 0.964039, 0.965984, 0.967994, 0.969968, 0.971283, 0.972800, 0.974223, 0.975784, 0.976647, 
            0.977953, 0.978887, 0.979922, 0.980773, 0.981637, 0.982386, 0.983068, 0.983759, 0.984436, 0.984855, 
            0.985453, 0.985916, 0.986417, 0.986934, 0.987182, 0.987435, 0.988022, 0.988146, 0.988537, 0.988792, 
            0.989043, 0.989235, 0.989470, 0.989681, 0.989857, 0.990006, 0.990159, 0.990263, 0.990455, 0.990515, 
            0.990636, 0.990731, 0.990787, 0.990884, 0.990900, 0.990971, 0.990986, 0.991042, 0.991048, 0.991057, 
            0.991057, 0.991048, 0.991042, 0.990986, 0.990971, 0.990900, 0.990884, 0.990787, 0.990731, 0.990636, 
            0.990515, 0.990455, 0.990263, 0.990159, 0.990006, 0.989857, 0.989681, 0.989470, 0.989235, 0.989043, 
            0.988787, 0.988532, 0.988146, 0.988022, 0.987435, 0.987182, 0.986934, 0.986417, 0.985916, 0.985453, 
            0.984855, 0.984436, 0.983759, 0.983068, 0.982386, 0.981637, 0.980773, 0.979920, 0.978887, 0.977953, 
            0.976647, 0.975784, 0.974223, 0.972800, 0.971284, 0.969970, 0.967994, 0.965984, 0.964039, 0.961755, 
            0.959339, 0.956675, 0.953860, 0.950768, 0.947313, 0.943788, 0.939752, 0.935424, 0.930744, 0.925504, 
            0.920103, 0.914023, 0.907216, 0.900084, 0.892242, 0.883399, 0.874277, 0.863766, 0.852416, 0.839932, 
            0.826496, 0.811669, 0.795374, 0.777313, 0.757467, 0.735991, 0.712150, 0.685667, 0.657029, 0.625035, 
            0.589915, 0.551116, 0.508397, 0.460966, 0.408796, 0.351055, 0.287226, 0.216842, 0.139913, 0.062742]]);

    beta = beta * DTOR                 # Tilt from horizontal of the PV modules/panels, in radians
    sazm = sazm * DTOR                 # Surface azimuth of PV module/panels, in radians

    # 1. Calculate and assign various paramters to be used for modeling irradiances
    iso_dif = 0.0; circ_dif = 0.0; horiz_dif = 0.0; grd_dif = 0.0; beam = 0.0;   # For calling PerezComp to break diffuse into components for zero tilt (horizontal)                           
    ghi, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, zen, 0.0, zen) # Call to get iso_dif for horizontal surface
    #            print "PEREZCOMP1 = "
    #            print "ghi = ", ghi
    #            print "iso_dif = ", iso_dif
    #            print "circ_dif = ", circ_dif
    #            print "horiz_dif = ", horiz_dif
    #            print "grd_dif = ", grd_dif
    #            print "beam = ", beam

    iso_sky_dif = iso_dif;       # Isotropic irradiance from sky on horizontal surface, used later for determining isotropic sky component

    inc, tiltr, sazmr = sunIncident(0, 90.0, 180.0, 45.0, zen, azm) # For calling PerezComp to break diffuse into components for 90 degree tilt (vertical)
    #            print "sunIncident 1."
    #            print "inc = ", inc
    #            print "tiltr = ", tiltr
    #            print "sazmr = ", sazmr

    vti, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen) # Call to get horiz_dif for vertical surface
    #            print "PEREZCOMP1 = "
    #            print "vti = ", vti
    #            print "iso_dif = ", iso_dif
    #            print "circ_dif = ", circ_dif
    #            print "horiz_dif = ", horiz_dif
    #            print "grd_dif = ", grd_dif
    #            print "beam = ", beam

    F2DHI = horiz_dif;           # Horizon diffuse irradiance on a vertical surface, used later for determining horizon brightening irradiance component

    index = -99;
    n2 = -99.9;
    if (PVfrontSurface == "glass"):
    
        index = 0;          # Index to use with 1-degree hemispherical segment AOI correction factor array
        n2 = 1.526;         # Index of refraction for glass
    
    elif (PVfrontSurface == "ARglass"):
    
        index = 1;          # Index to use with 1-degree hemispherical segment AOI correction factor array
        n2 = 1.300;         # Index of refraction for ARglass
    
    else:
        raise Exception("Incorrect text input for PVfrontSurface. Must be glass or ARglass.")
    
    Ro = math.pow((n2 - 1.0) / (n2 + 1.0), 2.0);     # Reflectance at normal incidence, Duffie and Beckman p217

    aveGroundGHI = 0.0;          # Average GHI on ground under PV array for cases when x projection exceed 2*rtr
    for i in range (0,100):
        aveGroundGHI += frontGroundGHI[i] / 100.0;

    # Calculate x,y coordinates of bottom and top edges of PV row in front of desired PV row so that portions of sky and ground viewed by the 
    # PV cell may be determined. Origin of x-y axis is the ground pobelow the lower front edge of the desired PV row. The row in front of 
    # the desired row is in the negative x direction.

    h = math.sin(beta);          # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    x1 = math.cos(beta);         # Horizontal distance from front of panel to rear of panel (in PV panel slope lengths)
    rtr = D + x1;                # Row-to-row distance (in PV panel slope lengths)
    PbotX = -rtr;                # x value for poon bottom egde of PV module/panel of row in front of (in PV panel slope lengths)
    PbotY = C;                   # y value for poon bottom egde of PV module/panel of row in front of (in PV panel slope lengths)
    PtopX = -D;                  # x value for poon top egde of PV module/panel of row in front of (in PV panel slope lengths)
    PtopY = h + C;               # y value for poon top egde of PV module/panel of row in front of (in PV panel slope lengths)

    # 2. Calculate diffuse and direct component irradiances for each cell row 
    
    
    for i in range (0, cellRows):
    
        # Calculate diffuse irradiances and reflected amounts for each cell row over it's field of view of 180 degrees, 
        # beginning with the angle providing the upper most view of the sky (j=0)
        PcellX = x1 * (i + 0.5) / (cellRows);                    # x value for location of PV cell
        PcellY = C + h * (i + 0.5) / (cellRows);                 # y value for location of PV cell
        elvUP = math.atan((PtopY - PcellY) / (PcellX - PtopX));          # Elevation angle up from PV cell to top of PV module/panel, radians
        elvDOWN = math.atan((PcellY - PbotY) / (PcellX - PbotX));        # Elevation angle down from PV cell to bottom of PV module/panel, radians
        if (rowType == "first" or rowType == "single"):                         # 4/19/16 No array in front for these cases
        
            elvUP = 0.0;
            elvDOWN = 0.0;
        
        #Console.WriteLine("ElvUp = 0", elvUP / DTOR);
        #if (i == 0)
        #    Console.WriteLine("ElvDown = 0", elvDOWN / DTOR);

        if math.isnan(beta):
            print( "Beta is Nan")
        if math.isnan(elvUP):
            print( "elvUP is Nan")
        if math.isnan((math.pi - beta - elvUP) / DTOR):
            print( "division is Nan")
        
        
        iStopIso = int(round(np.float64((math.pi - beta - elvUP)) / DTOR)) # Last whole degree in arc range that sees sky, first is 0
        #Console.WriteLine("iStopIso = 0", iStopIso);
        iHorBright = int(round(max(0.0, 6.0 - elvUP / DTOR)));    # Number of whole degrees for which horizon brightening occurs
        iStartGrd = int(round((math.pi - beta + elvDOWN) / DTOR));     # First whole degree in arc range that sees ground, last is 180
    #                print "iStopIso = ", iStopIso
    #                print "iHorBright = ", iHorBright
    #                print "iStartGrd = ", iStartGrd

        frontGTI.append(0.0)                                                      # Initialtize front GTI
        frontReflected.append(0.0);                                                # Initialize reflected amount from front

        for j in range (0, iStopIso):                                        # Add sky diffuse component and horizon brightening if present
        #for (j = 0; j < iStopIso; j++)                                      
        
            frontGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * iso_sky_dif;                               # Sky radiation
            frontReflected[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * iso_sky_dif * (1.0 - SegAOIcor[index][j] * (1.0 - Ro));    # Reflected radiation from module
            if ((iStopIso - j) <= iHorBright):                                   # Add horizon brightening term if seen
            
                frontGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * F2DHI / 0.052264;  # 0.052246 = 0.5 * [cos(84) - cos(90)]
                frontReflected[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * (F2DHI / 0.052264) * (1.0 - SegAOIcor[index][j] * (1.0 - Ro));    # Reflected radiation from module
            
        

        #if (i == 0)
        #    Console.WriteLine("iStartGrd = 0", iStartGrd);
        for j in range (iStartGrd, 180):                                     # Add ground reflected component
        #(j = iStartGrd; j < 180; j++)                                   
            startElvDown = (j - iStartGrd) * DTOR + elvDOWN;             # Start and ending down elevations for this j loop 
            stopElvDown = (j + 1 - iStartGrd) * DTOR + elvDOWN;
            projectedX1 = PcellX - np.float64(PcellY) / math.tan(startElvDown);      # Projection of ElvDown to ground in -x direction
            projectedX2 = PcellX - PcellY / math.tan(stopElvDown);
            actualGroundGHI = 0.0;                                       # Actuall ground GHI from summing array values
            #if (i == 0)
            #    Console.WriteLine("j= 0 projected X1 = 1,6:0.0", j, 100 * projectedX1 / rtr);
            if (abs(projectedX1 - projectedX2) > 0.99 * rtr):
            
                if (rowType == "first" or rowType == "single"):                  # 4/19/16 No array in front for these cases
                
                    actualGroundGHI = ghi;                                      # Use total value if projection approximates the rtr
                
                else:
                    actualGroundGHI = aveGroundGHI;                                 # Use average value if projection approximates the rtr                        
            
            else:
            
                projectedX1 = 100.0 * projectedX1 / rtr;                        # Normalize projections and multiply by 100
                projectedX2 = 100.0 * projectedX2 / rtr;
                if ((rowType == "first" or rowType == "single") and (abs(projectedX1) > rtr or abs(projectedX2) > rtr)):    #4/19/2016
                
                    actualGroundGHI = ghi;                                      # Use total value if projection > rtr for "first" or "single"
                
                else:
                
                    while (projectedX1 < 0.0 or projectedX2 < 0.0):                  # Offset so array indexes are positive
                    
                        projectedX1 += 100.0;
                        projectedX2 += 100.0;
                    
                    index1 = int(projectedX1);                                  # Determine indexes for use with groundGHI array (truncates values)
                    index2 = int(projectedX2);
                    if (index1 == index2):
                    
                        actualGroundGHI = frontGroundGHI[index1];                        # x projections in same groundGHI element
                    
                    else:
                    
                        for k in range (index1, index2+1):                   # Sum the irradiances on the ground if projections are in different groundGHI elements
                        #for (k = index1; k <= index2; k++)
                        
                            #Console.WriteLine("index1=0 index2=1", index1,index2);
                            if (k == index1):
                                actualGroundGHI += frontGroundGHI[k] * (k + 1.0 - projectedX1);
                            elif (k == index2):
                            
                                if (k < 100):
                                    actualGroundGHI += frontGroundGHI[k] * (projectedX2 - k);
                                else:
                                    actualGroundGHI += frontGroundGHI[k - 100] * (projectedX2 - k);
                            
                            else:
                            
                                if (k < 100):
                                    actualGroundGHI += frontGroundGHI[k];
                                else:
                                    actualGroundGHI += frontGroundGHI[k - 100];
                            
                        
                        actualGroundGHI /= projectedX2 - projectedX1;                # Irradiance on ground in the 1 degree field of view
                    
                
                #if (i == 0)
                #    Console.WriteLine("j=0 index1=1 index2=2 projectX1=3,5:0.0 projectX2=4,5:0.0 actualGrdGHI=5,6:0.0", j, index1, index2, projectedX1, projectedX2, actualGroundGHI);
            
            frontGTI[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * SegAOIcor[index][j] * actualGroundGHI * albedo;     # Add ground reflected component
            frontReflected[i] += 0.5 * (math.cos(j * DTOR) - math.cos((j + 1) * DTOR)) * actualGroundGHI * albedo * (1.0 - SegAOIcor[index][j] * (1.0 - Ro));    # Reflected ground radiation from module
            #Console.WriteLine("actualGroundGHI = 0,6:0.0 inputGHI = 1,6:0.0 aveArrayGroundGHI = 2,6:0.0", actualGroundGHI, dhi + dni * math.cos(zen), aveGroundGHI);
            # End of j loop for adding ground reflected componenet 

        # Calculate and add direct and circumsolar irradiance components
        inc, tiltr, sazmr = sunIncident(0, beta / DTOR, sazm / DTOR, 45.0, zen, azm) # For calling PerezComp to break diffuse into components for 90 degree tilt (vertical)
    #                print "sunIncident 2."
    #                print "inc = ", inc
    #                print "tiltr = ", tiltr
    #                print "sazmr = ", sazmr

    #                print " INCIDENT REALY NEEDED for AOI ", inc
        gtiAllpc, iso_dif, circ_dif, horiz_dif, grd_dif, beam = perezComp(dni, dhi, albedo, inc, tiltr, zen) # Call to get components for the tilt
    #                print "PEREZCOMP 2 = "
    #                print "gtiAllpc = ", vti
    #                print "iso_dif = ", iso_dif
    #                print "circ_dif = ", circ_dif
    #                print "horiz_dif = ", horiz_dif
    #                print "grd_dif = ", grd_dif
    #                print "beam = ", beam

        cellShade = pvFrontSH * cellRows - i;
        if (cellShade > 1.0):    # Fully shaded if > 1, no shade if < 0, otherwise fractionally shaded
            cellShade = 1.0;
        elif (cellShade < 0.0):
            cellShade = 0.0;

        if (cellShade < 1.0 and inc < math.pi / 2.0):  # Cell not shaded entirely and inc < 90 deg
        
            cor = aOIcorrection(n2, inc);                # Get AOI correction for beam and circumsolar
            frontGTI[i] += (1.0 - cellShade) * (beam + circ_dif) * cor; # Add beam and circumsolar radiation
            #frontReflected[i] += (1.0 - cellShade) * (beam + circ_dif) * (1.0 - cor * (1.0 - Ro));    # Reflected beam and circumsolar radiation from module
        
        # End of for i = 0; i < cellRows loop
    return aveGroundGHI, frontGTI, frontReflected;
    # End of GetFrontSurfaceIrradiances

    
def getGroundShadeFactors(rowType, beta, C, D, elv, azm, sazm):
    """
    This method determines if the ground is shaded from direct beam radiation
    for points on the ground from the leading edge of one row of PV panels to
    the leading edge of the next row of PV panels behind it. This row-to-row
    dimension is divided into 100 ground segments and a ground shade factor is
    returned for each ground segment, with values of 1 for shaded segments and
    values of 0 for non shaded segments. The fractional amounts of shading of
    the front and back surfaces of the PV panel are also returned. 8/20/2015

    4/18/2016 - Modified to account for different row types. Because the ground
    factors may now be different depending on row, they are calculated for the
    row-to-row dimension to the rear of the leading module edge and to the
    front of the leading edge. Also returned is the maximum shadow length
    projected to the front or rear from the front of the module row

    Parameters
    ----------
    rowType : str
        "first", "interior", "last", or "single" 
    beta
        Tilt from horizontal of the PV modules/panels (deg)
    C
        Ground clearance of PV panel (in PV panel slope lengths)            
    D
        Horizontal distance between rows of PV panels (in PV panel slope
        lengths) 
    elv
        Sun elevation (in radians)
    azm
        Sun azimuth (in radians)
    sazm
        Surface azimuth of PV panels (deg)

    Returns
    -------
    pvFrontSH : numeric
        Decimal fraction of the front surface of the PV panel that is shaded,
        0.0 to 1.0
    pvBackSH : numeric
        Decimal fraction of the back surface of the PV panel that is shaded,
        0.0 to 1.0
    rearGroundSH : array of size [100]
        Ground shade factors for ground segments to the rear, 0 = not shaded,
        1 = shaded
    frontGroundSH : array of size [100]
        Ground shade factors for ground segments to the front, 0 = not shaded,
        1 = shaded
    maxShadow : numeric
        Maximum shadow length projected to the front(-) or rear (+) from the
        front of the module row (in PV panel slope lengths), only used later
        for rowTypes other than "interior"
    """
    rearGroundSH = []
    frontGroundSH = []
    
    beta = beta * DTOR  # Tilt from horizontal of the PV modules/panels, in radians
    sazm = sazm * DTOR  # Surface azimuth of PV module/pamels, in radians

    h = math.sin(beta);          # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    x1 = math.cos(beta);         # Horizontal distance from front of panel to rear of panel (in PV panel slope lengths)
    rtr = D + x1;                # Row-to-row distance (in PV panel slope lengths)

    # Divide the row-to-row spacing into 100 intervals for calculating ground shade factors
    delta = rtr / 100.0;
    x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals

    Lh = (h / math.tan(elv)) * math.cos(sazm - azm); # Horizontal length of shadow perpindicular to row from top of module to bottom of module
    Lhc = ((h + C) / math.tan(elv)) * math.cos(sazm - azm); # Horizontal length of shadow perpindicular to row from top of module to ground level
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

    elif (rowType == "first"):
    
        if (Lh > 0.0):   # Sun is on front side of PV module
        
            pvFrontSH = 0.0;
            pvBackSH = 1.0;
            ss1 = Lc;           # Ground shaded from shadow of lower edge
            se1 = x1 + Lhc;     # to shadow of upper edge                                        
            # End of if sun on front side of PV module

        elif (Lh < -(rtr + x1)):  # Back side of PV module partially shaded from row to rear, front completely shaded, ground completely shaded
        
            pvFrontSH = 1.0;
            pvBackSH = (Lh + rtr + x1) / (Lh + x1);
            ss1 = -rtr;      # Ground shaded from -rtr to rtr
            se1 = rtr;

            # End of if back side of PV module partially shaded, front completely shaded, ground completely shaded

        else:   # Shadow to frontside of row, either front or back might be shaded, depending on tilt and other factors     
        
            if (Lc < Lhc + x1):
            
                pvFrontSH = 0.0;
                pvBackSH = 1.0;
                ss1 = Lc;         # Shadow starts at Lc
                se1 = Lhc + x1;   # Shadow ends here
            
            else:
            
                pvFrontSH = 1.0;
                pvBackSH = 0.0;
                ss1 = Lhc + x1;  # Shadow starts at Lhc + x1
                se1 = Lc;        # Shadow ends here
            

            # End of shadow to front of row 

        delta = rtr / 100.0;
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals            
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                rearGroundSH.append(1)         # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                rearGroundSH.append(0)         # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        

        x = -rtr - delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals for front interval           
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                frontGroundSH.append(1)        # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                frontGroundSH.append(0)        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        
        # End of if row type == "first"

    elif (rowType == "last"):
    
        if (Lh > D): # Front side of PV module partially shaded, back completely shaded, ground completely shaded
        
            pvFrontSH = (Lh - D) / (Lh + x1);
            pvBackSH = 1.0;
            ss1 = -rtr;      # Ground shaded from -rtr to rtr
            se1 = rtr;
        
        else:   # Shadow to frontside of row, either front or back might be shaded, depending on tilt and other factors     
        
            if (Lc < Lhc + x1):
            
                pvFrontSH = 0.0;
                pvBackSH = 1.0;
                ss1 = Lc;         # Shadow starts at Lc
                se1 = Lhc + x1;   # Shadow ends here
            
            else:
            
                pvFrontSH = 1.0;
                pvBackSH = 0.0;
                ss1 = Lhc + x1;  # Shadow starts at Lhc + x1
                se1 = Lc;        # Shadow ends here
            

            # End of shadow to front of row 

        delta = rtr / 100.0;
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals            
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                rearGroundSH.append(1);        # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                rearGroundSH.append(0);        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        

        x = -rtr - delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals for front interval           
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                frontGroundSH.append(1);        # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                frontGroundSH.append(0);        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        

        # End of if row type == "last"

    elif (rowType == "single"):
    
        if (Lh > 0.0):   # Shadow to the rear
        
            pvFrontSH = 0.0;
            pvBackSH = 1.0;
            ss1 = Lc;           # Ground shaded from shadow of lower edge
            se1 = x1 + Lhc;     # to shadow of upper edge                                        
            # End of if sun on front side of PV module

        else:   # Shadow to frontside of row, either front or back might be shaded, depending on tilt and other factors     
        
            if (Lc < Lhc + x1):
            
                pvFrontSH = 0.0;
                pvBackSH = 1.0;
                ss1 = Lc;         # Shadow starts at Lc
                se1 = Lhc + x1;   # Shadow ends here
            
            else:
            
                pvFrontSH = 1.0;
                pvBackSH = 0.0;
                ss1 = Lhc + x1;  # Shadow starts at Lhc + x1
                se1 = Lc;        # Shadow ends here
            

            # End of shadow to front of row 

        delta = rtr / 100.0;
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals            
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                rearGroundSH.append(1);        # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                rearGroundSH.append(0);        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        

        x = -rtr - delta / 2.0;    # Initialize horizontal dimension x to provide midpoof intervals for front interval           
        for i in range(0,100):
        
            x += delta;
            if (x >= ss1 and x < se1):
                frontGroundSH.append(1);        # x within a shaded interval, set groundSH to 1 to indicate shaded
            else:
                frontGroundSH.append(0);        # x not within a shaded interval, set groundSH to 0 to indicated not shaded, i.e. sunny
        

        # End of if row type == "single"
    else:
        print ("ERROR: Incorrect row type not passed to function GetGroundShadedFactors ");

    if (abs(ss1) > abs(se1)):      # Maximum shadow length projected from the front of the PV module row
        maxShadow = ss1;
    else:
        maxShadow = se1;

    #Console.WriteLine("elv = 0,6:0.00  azm = 1,6:0.00  sazm = 2,6:0.00", elv * 180.0 / math.pi, azm * 180.0 / math.pi, sazm * 180.0 / math.pi);
    #Console.WriteLine("ss1 = 0,6:0.0000 se1 = 1,6:0.0000 ss2 = 2,6:0.0000 se2 = 3,6:0.0000     rtr = 4,6:0.000", ss1, se1, ss2, se2, rtr);
    #Console.WriteLine("pvFrontSH = 0,6:0.00 pvBackSH = 1,6:0.00", pvFrontSH, pvBackSH);

    # End of GetGroundShadedFactors
    
    #print "rearGroundSH", rearGroundSH[0]
    return pvFrontSH, pvBackSH, maxShadow, rearGroundSH, frontGroundSH;
    # End of getGroundShadeFactors


def getSkyConfigurationFactors(rowType, beta, C, D):
    """
    This method determines the sky configuration factors for points on the
    ground from the leading edge of one row of PV panels to the leading edge of
    the next row of PV panels behind it. This row-to-row dimension is divided
    into 100 ground segments and a sky configuration factor is returned for
    each ground segment. The sky configuration factor represents the fraction
    of the isotropic diffuse sky radiation (unobstructed) that is present on
    the ground when partially obstructed by the rows of PV panels. The
    equations follow that on pages in the notebook dated 8/12/2015. 8/20/2015            

    4/15/2016 Modifed for calculations other than just the interior rows. Row
    type is identified with the string `rowType`, with the possilbe values:

    * first = first row of the array
    * interior = interior row of array
    * last = last row of the array
    * single = a single row array

    Because the sky configuration factors may now be different depending on
    row, they are calculated for the row-to-row dimension to the rear of the
    leading module edge and to the front of the leading edge.

    Parameters
    ----------
    rowType : str
        "first", "interior", "last", or "single"
    beta : float
        Tilt from horizontal of the PV modules/panels (deg)
    C : float
        Ground clearance of PV panel (in PV module/panel slope lengths)            
    D : float
        Horizontal distance between rows of PV panels (in PV module/panel slope
        lengths)

    Returns
    -------
    rearSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)
    frontSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)

    Notes
    -----
    The horizontal distance between rows, `D`, is from the back edge of one row
    to the front edge of the next, and it is not the row-to-row spacing.
    """

    rearSkyConfigFactors = []
    frontSkyConfigFactors = []
    
    # Tilt from horizontal of the PV modules/panels, in radians
    beta = beta * DTOR
    # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    h = math.sin(beta)
    # Horizontal distance from front of panel to rear of panel (in PV panel
    # slope lengths)
    x1 = math.cos(beta)
    rtr = D + x1  # Row-to-row distance (in PV panel slope lengths)

    # Forced fix for case of C = 0
    # FIXME: for some reason the Config Factors go from 1 to 2 and not 0 to 1.
    # TODO: investigate why this is happening in the code.
    if C==0:
        C=0.0000000001
        
    if C < 0:
        LOGGER.error(
            "Height is below ground level. Function GetSkyConfigurationFactors"
            " will continue but results might be unreliable")
        
    # Divide the row-to-row spacing into 100 intervals and calculate
    # configuration factors
    delta = rtr / 100.0

    if (rowType == "interior"):
        # Initialize horizontal dimension x to provide midpoint of intervals
        x = -delta / 2.0

        for i in range(0,100):        
            x += delta
            # all dimensions are relative to module height
            # in other words, consider a module height of 1 arb. unit
            #
            #   <--rtr=x1+D--><--rtr=x1+D--><--rtr=x1+D-->
            #  :\\     row   :\\     row   :\\     row   :\\
            #  : \\`    0    : \\     1    : \\     2   /: \\
            #  h  \\  `      h  \\         h  \\      /  h  \\
            #  :   \\    `   :   \\        :   \\   /    :   \\ pi-beta
            #  :.x1.\\...D..`:.x1.\\...D...:.x1.\\/...D..:.x1.\\...
            #  :               `   <------x-----/|
            #  C                  `           /
            #  :              angA   `      /  angB
            #  +------------------------`-/---------------------
            #                            x
            # use ATAN2: 4-quadrant tangent instead of ATAN
            # check 2 rows away
            angA = math.atan2(h + C, (2.0 * rtr + x1 - x))
            angB = math.atan2(C, (2.0 * rtr - x))
            beta1 = max(angA, angB)
            # check 1 rows away
            angA = math.atan2(h + C, (rtr + x1 - x))
            angB = math.atan2(C, (rtr - x))
            beta2 = min(angA, angB)
            # check 0 rows away
            beta3 = max(angA, angB)
            beta4 = math.atan2(h + C, (x1 - x))
            beta5 = math.atan2(C, (-x))
            beta6 = math.atan2(h + C, (-D - x))
            sky1 =0; sky2 =0; sky3 =0
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2))
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4))
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6))
            skyAll = sky1 + sky2 + sky3

            # Save as arrays of values, same for both to the rear and front
            rearSkyConfigFactors.append(skyAll)
            frontSkyConfigFactors.append(skyAll)        
        # End of if "interior"

    elif (rowType == "first"):
    
        # RearSkyConfigFactors don't have a row in front, calculation of sky3
        # changed, beta6 = 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);    # Save as arrays of values                   
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  frontSkyConfigFactors don't have a row in front, calculation of sky3 included as part of revised sky2,
        #  beta 4 set to 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));

            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "first"

    elif (rowType == "last"):
    
        #  RearSkyConfigFactors don't have a row to the rear, combine sky1 into sky 2, set beta 3 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have beta1 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            frontSkyConfigFactors.append(skyAll);      # Save as arrays of values,
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "last" row

    elif (rowType == "single"):
    
                    
        #  RearSkyConfigFactors don't have a row to the rear ir front, combine sky1 into sky 2, set beta 3 = 0.0,
        #  for sky3, beta6 = 180.0.
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values                    
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        

        #  FrontSkyConfigFactors have only a row to the rear, combine sky3 into sky2, set beta1 = 0, beta4 = 180
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):
        
            x += delta;
            angA = math.atan((h + C) / (2.0 * rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (rtr + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (rtr - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
        
        # End of if "single"
    else:
        print("ERROR: Incorrect row type not passed to function GetSkyConfigurationFactors ");

    return rearSkyConfigFactors, frontSkyConfigFactors;
# End of GetSkyConfigurationFactors


def rowSpacing(beta, sazm, lat, lng, tz, hour, minute):
    """
    This method determines the horizontal distance D between rows of PV panels
    (in PV module/panel slope lengths) for no shading on December 21 (north
    hemisphere) June 21 (south hemisphere) for a module tilt angle beta and
    surface azimuth sazm, and a given latitude, longitude, and time zone and
    for the time passed to the method (typically 9 am).

    (Ref: the row-to-row spacing is then ``D + cos(beta)``)
    8/21/2015

    Parameters
    ----------
    beta : double
        Tilt from horizontal of the PV modules/panels (deg)
    sazm : double
        Surface azimuth of the PV modules/panels (deg)
    lat : double
        Site latitude (deg)
    lng : double
        Site longitude (deg)
    tz : double
        Time zone (hrs)
    hour : int
        hour for no shading criteria
    minute: double
        minute for no shading

    Returns
    -------
    D : numeric
        Horizontal distance between rows of PV panels (in PV panel slope
        lengths)
    """     
    beta = beta * DTOR  # Tilt from horizontal of the PV modules/panels, in radians
    sazm = sazm * DTOR  # Surface azimuth of PV module/pamels, in radians
    if lat >= 0:
        [azm, zen, elv, dec, sunrise, sunset, Eo, tst] = solarPos (2014, 12, 21, hour, minute, lat, lng, tz)
    else:
        [azm, zen, elv, dec, sunrise, sunset, Eo, tst] = solarPos (2014, 6, 21, hour, minute, lat, lng, tz)
    tst = 8.877  ##DLL Forced value
    minute -= 60.0 * (tst - hour);      # Adjust minute so sun position is calculated for a tst equal to the
      # time passed to the function

    if lat >= 0:
        [azm, zen, elv, dec, sunrise, sunset, Eo, tst] = solarPos(2014, 12, 21, hour, minute, lat, lng, tz)
    else:
        [azm, zen, elv, dec, sunrise, sunset, Eo, tst] = solarPos(2014, 6, 21, hour, minute, lat, lng, tz)
      
    # Console.WriteLine("tst = {0} azm = {1} elv = {2}", tst, azm * 180.0 / Math.PI, elv * 180.0 / Math.PI);
    D = math.cos(sazm - azm) * math.sin(beta) / math.tan(elv)
    return D
# End of RowSpacing  


def trackingBFvaluescalculator(beta, hub_height, r2r):
    '''
    1-axis tracking helper file

    Parameters
    ----------
    beta : series of floats, or float 
        Tilt from horizontal of the PV modules/panels, in radians
    hub_height : float
        tracker hub height
    r2r : float
        Row-to-row distance (in PV panel slope lengths)

    Returns
    -------
    C : float
        ground clearance of PV panel
    D : float
        row-to-row distance (each in PV panel slope lengths)
    '''

    
    beta = beta * DTOR  # Tilt from horizontal of the PV modules/panels, in radians

    try:
        x1 = beta.apply(math.cos)   # Horizontal distance from front of panel to rear of panel (in PV panel slope lengths)
        D = r2r - x1                # Calculates D DistanceBetweenRows(panel slope lengths)
        hm = beta.apply(math.sin)*0.5        # vertical distance from bottom of panel to top of panel (in PV panel slope lengths)
        #C = 0.5+Cv-hm               # Ground clearance of PV panel (in PV panel slope lengths). 
        C = hub_height - hm          #Adding a 0.5 for half a panel slope length, since it is assumed the panel is rotating around its middle axis 
    except: # in case only one value is passed
        x1 = math.cos(beta);         # Horizontal distance from front of panel to rear of panel (in PV panel slope lengths)
        #rtr = D + x1;                # Row-to-row distance (in PV panel slope lengths)
        D = r2r - x1;                # Calculates D DistanceBetweenRows(panel slope lengths)
        hm = 0.5*math.sin(beta);        # vertical distance from bottom of panel to top of panel (in PV panel slope lengths)
        #C = 0.5+Cv-hm               # Ground clearance of PV panel (in PV panel slope lengths). 
        C = hub_height - hm          #Adding a 0.5 for half a panel slope length, since it is assumed the panel is rotating around its middle axis 

    return C, D



def getSkyConfigurationFactors2(rowType, beta, C, D, pitch):
    """
    Converting inputs to float or series!
    
    This method determines the sky configuration factors for points on the
    ground from the leading edge of one row of PV panels to the leading edge of
    the next row of PV panels behind it. This row-to-row dimension is divided
    into 100 ground segments and a sky configuration factor is returned for
    each ground segment. The sky configuration factor represents the fraction
    of the isotropic diffuse sky radiation (unobstructed) that is present on
    the ground when partially obstructed by the rows of PV panels. The
    equations follow that on pages in the notebook dated 8/12/2015. 8/20/2015            
    4/15/2016 Modifed for calculations other than just the interior rows. Row
    type is identified with the string `rowType`, with the possilbe values:
    * first = first row of the array
    * interior = interior row of array
    * last = last row of the array
    * single = a single row array
    Because the sky configuration factors may now be different depending on
    row, they are calculated for the row-to-row dimension to the rear of the
    leading module edge and to the front of the leading edge.
    Parameters
    ----------
    rowType : str
        "first", "interior", "last", or "single"
    beta : float
        Tilt from horizontal of the PV modules/panels (deg)
    C : float
        Ground clearance of PV panel (in PV module/panel slope lengths)            
    D : float
        Horizontal distance between rows of PV panels (in PV module/panel slope
        lengths)
    Returns
    -------
    rearSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)
    frontSkyConfigFactors : array of size [100]
        Sky configuration factors to rear of leading PVmodule edge (decimal
        fraction)
    Notes
    -----
    The horizontal distance between rows, `D`, is from the back edge of one row
    to the front edge of the next, and it is not the row-to-row spacing.
    """

    rearSkyConfigFactors = []
    frontSkyConfigFactors = []

    # Tilt from horizontal of the PV modules/panels, in radians
    beta = beta * DTOR
    # Vertical height of sloped PV panel (in PV panel slope lengths)                    
    h = beta.apply(math.sin)
    # Horizontal distance from front of panel to rear of panel (in PV panel
    # slope lengths)
    x1 = beta.apply(math.cos)

    # Forced fix for case of C = 0
    # FIXME: for some reason the Config Factors go from 1 to 2 and not 0 to 1.
    # TODO: investigate why this is happening in the code. #Sil: Did we fix this?        
    C[C==0] = 0.0000000001

    if not C[C<0].empty:    
        LOGGER.error("ERROR: Clearance height is below ground level."
                     " Function GetSkyConfigurationFactors "
                     " will continue but results might be unreliable")
        C[C<0] = 0.0000000001

    # Divide the row-to-row spacing into 100 intervals and calculate
    # configuration factors
    delta = pitch / 101.0
    x=np.linspace(delta, pitch-delta, num=100, endpoint=True)
    df_x1 = pd.DataFrame(zip(*[x1 for i in range(100)]))

    if (rowType == "interior"):
    # Initialize horizontal dimension x to provide midpoint of intervals

        #  <--rtr=x1+D--><--rtr=x1+D--><--rtr=x1+D-->
        # |\            |\            |\            |\ 
        # | \ `         | \           | \          /| \
        # h  \   `      h  \          h  \       /  h  \
        # |   \     `   |   \         |   \    /    |   \
        # |_x1_\____D__`|_x1_\____D___|_x1_\_/_D____|_x1_\_
        # |               `   <------x-----/|
        # C                  `           /
        # |              angA   `      /  angB
        # *------------------------`-/---------------------
        #                          x
        # use ATAN2: 4-quadrant tangent instead of ATAN
        # check 2 rows away

        #angA = math.atan2(h + C, (2.0 * pitch + x1 - x))
        # atan -- y/x .. (y, x)
        clear = h+C
        df_clear = pd.DataFrame(zip(*[clear for i in range(100)]))
        df = df_x1.subtract(x)  #x1-x
        df = df+pitch*2.0
        df_yx = df_clear/df
        angA = df_yx.applymap(math.atan)

        #angB = math.atan2(C, (2.0 * pitch - x))
        df_y = pd.DataFrame(zip(*[C for i in range(100)]))
        df = 2.0*pitch - x
        df_yx = df_y/df
        angB = df_yx.applymap(math.atan)

        beta1 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

        ## check 1 rows away
        #angA = math.atan2(h + C, (pitch + x1 - x))
        df = df_x1.subtract(x)  #x1-x
        df = df+pitch
        df_yx = df_clear/df
        angA = df_yx.applymap(math.atan)

        #angB = math.atan2(C, (pitch - x))
        df = pitch - x
        df_yx = df_y/df
        angB = df_yx.applymap(math.atan)

        # beta2 = min(angA, angB)
        beta2 = angA.where(angA < angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

        #beta3 = max(angA, angB)
        beta3 = angA.where(angA > angB, angB).fillna(angA) # selects the max of each element-wise between both dataframes.

        # check 0 rows away
        #beta4 = math.atan2(h + C, (x1 - x))
        df_yx = df_clear / df_x1.subtract(x)  #x1-x
        beta4 = df_yx.applymap(math.atan)

        #beta5 = math.atan2(C, (-x))
        df_yx = df_y / -x
        beta5 = df_yx.applymap(math.atan) 

        #beta6 = math.atan2(h + C, (-D - x))
        df_D = pd.DataFrame(zip(*[D for i in range(100)]))
        df_D = -df_D-x
        df_yx = df_clear / df_D
        beta6 = df_yx.applymap(math.atan)

        #SIL CONTINUE FROM HERE
        sky1 =0; sky2 =0; sky3 =0
        if (beta2 > beta1):
            sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2))
        if (beta4 > beta3):
            sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4))
        if (beta6 > beta5):
            sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6))
        skyAll = sky1 + sky2 + sky3

        # Save as arrays of values, same for both to the rear and front
        rearSkyConfigFactors.append(skyAll)
        frontSkyConfigFactors.append(skyAll)        
    # End of if "interior"

    elif (rowType == "first"):

        # RearSkyConfigFactors don't have a row in front, calculation of sky3
        # changed, beta6 = 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);    # Save as arrays of values                   
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);


        #  frontSkyConfigFactors don't have a row in front, calculation of sky3 included as part of revised sky2,
        #  beta 4 set to 180 degrees
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));

            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);

        # End of if "first"

    elif (rowType == "last"):

        #  RearSkyConfigFactors don't have a row to the rear, combine sky1 into sky 2, set beta 3 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);


        #  FrontSkyConfigFactors have beta1 = 0.0
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.atan((h + C) / (-D - x));
            if (beta6 < 0.0):
                beta6 += math.pi;

            sky1 = 0.0; sky2 = 0.0; sky3 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky1 + sky2 + sky3;
            frontSkyConfigFactors.append(skyAll);      # Save as arrays of values,
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);

        # End of if "last" row

    elif (rowType == "single"):


        #  RearSkyConfigFactors don't have a row to the rear ir front, combine sky1 into sky 2, set beta 3 = 0.0,
        #  for sky3, beta6 = 180.0.
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;

            beta3 = 0.0;

            beta4 = math.atan((h + C) / (x1 - x));
            if (beta4 < 0.0):
                beta4 += math.pi;

            beta5 = math.atan(C / (-x));
            if (beta5 < 0.0):
                beta5 += math.pi;

            beta6 = math.pi;

            sky2 = 0.0; sky3 = 0.0;
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            if (beta6 > beta5):
                sky3 = 0.5 * (math.cos(beta5) - math.cos(beta6));
            skyAll = sky2 + sky3;
            rearSkyConfigFactors.append(skyAll);     # Save as arrays of values                    
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);


        #  FrontSkyConfigFactors have only a row to the rear, combine sky3 into sky2, set beta1 = 0, beta4 = 180
        x = -delta / 2.0;    # Initialize horizontal dimension x to provide midpoint of intervals

        for i in range(0,100):

            x += delta;
            angA = math.atan((h + C) / (2.0 * pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (2.0 * pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta1 = max(angA, angB);
            beta1 = 0.0;

            angA = math.atan((h + C) / (pitch + x1 - x));
            if (angA < 0.0):
                angA += math.pi;
            angB = math.atan(C / (pitch - x));
            if (angB < 0.0):
                angB += math.pi;
            beta2 = min(angA, angB);

            beta3 = max(angA, angB);

            beta4 = math.pi;

            sky1 = 0.0; sky2 = 0.0;
            if (beta2 > beta1):
                sky1 = 0.5 * (math.cos(beta1) - math.cos(beta2));
            if (beta4 > beta3):
                sky2 = 0.5 * (math.cos(beta3) - math.cos(beta4));
            skyAll = sky1 + sky2;
            frontSkyConfigFactors.append(skyAll);  # Save as arrays of values
            #Console.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);
            #sw.WriteLine("0,5:0.000,1,5:0.000,2,5:0.000,3,5:0.000,4,5:0.000", x, sky1, sky2, sky3, skyAll);

        # End of if "single"
    else:
        print("ERROR: Incorrect row type not passed to function GetSkyConfigurationFactors ");

    return rearSkyConfigFactors, frontSkyConfigFactors;

