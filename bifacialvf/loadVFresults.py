#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 14:18:48 2017


"""
from __future__ import division, print_function # ensure python3 compatible division and printing
import io
try:
    from urllib2 import urlopen
except ImportError:
    from urllib.request import urlopen

import pandas as pd



def loadVFresults(filename=None):
    '''
    Read a VF CSV Result  file in to a pandas dataframe.

    Parameters
    ----------
    filename : None or string
        If None, attempts to use a Tkinter file browser. A string can be
        a relative file path, absolute file path, or url.

    Returns
    -------
    Tuple of the form (data, metadata).

    data : DataFrame
        A pandas dataframe with the columns described in the table
        below. For more detailed descriptions of each component, please
        consult the TMY3 User's Manual ([1]), especially tables 1-1
        through 1-6.

    metadata : dict
        The site metadata available in the file.

    Notes
    -----

    

    '''

    if filename is None:
        try:
            filename = _interactive_load()
        except:
            raise Exception('Interactive load failed. Tkinter not supported ' +
                            'on this system. Try installing X-Quartz and ' +
                            'reloading')

    
    try:
        csvdata = open(filename, 'r')
    except IOError:
        response = urlopen(filename)
        csvdata = io.StringIO(response.read().decode(errors='ignore'))
    
    # read in file metadata
    headerline = csvdata.readline() # read headers labels, read values
    head1 = headerline.strip().split(',')
    meta = dict(zip(head1, csvdata.readline().rstrip('\n').split(",")))
    
 
    
    import numpy as np
    data = pd.read_csv(filename, skiprows=1, header=1)

    if 'Ground Irradiance Values' in data.columns:
        data['Ground Irradiance Values'] = data['Ground Irradiance Values'].apply(
            lambda s: np.fromstring(s.strip('[]').replace(',', ' '), sep=' ') if isinstance(s, str) else s)

    return data, meta
#
#class Program:
#
#    data, meta = loadVFresults('data//output//LOADTEST.csv')    
#    print meta
def _interactive_load():
    import Tkinter
    from tkFileDialog import askopenfilename
    Tkinter.Tk().withdraw()  # Start interactive file input
    return askopenfilename()