# -*- coding: utf-8 -*-
"""
Created on Tue Jun 15 13:17:20 2021

@author: snowy
"""

import re
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
#First we need to make all Headers into txt folders
TimeH = []
TimeM = []
TimeS = []
RA_D1 = []
DEC_D1 = []
COMMENTS = []

ABBA_pattern = 'a', 'b', 'b', 'a'
n = 0
No_list = np.arange(87) + 41
warnings.filterwarnings('ignore', 'Header block contains null bytes instead of spaces for padding',)
warnings.filterwarnings('ignore', 'non-ASCII characters are present in the FITS file header',)
for x in range(87):
    xi = No_list[x]
    if xi < 10:
        filename = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s000' + str(xi) + '.fits'
    elif xi >= 10 and xi < 100: 
        filename = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s00' + str(xi) + '.fits'
    else:
        filename = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s0' + str(xi) + '.fits'
    hdu = fits.open(filename, ignore_missing_end=True)
    hdr = hdu[0].header
    Time = str(hdr['UTC'])
    RA = hdr['RA']
    RA_D1.append(RA)
    DEC = hdr['DEC']
    DEC_D1.append(DEC)
    Time = re.findall(r'\d+', Time)
    TimeH.append(int(Time[0]))
    TimeM.append(int(Time[1]))
    TimeS.append(int(Time[2]))
    
Deg_RA = []
Deg_DEC = []
for a in range(84):
    Full_Ra = RA_D1[a]
    Deg_RA.append(Full_Ra)
    Full_DEC = DEC_D1[a]
    Deg_DEC.append(Full_DEC)

plt.figure()
plt.plot(np.arange(84), Deg_RA)
plt.ylabel('Right Ascension in Degrees')
plt.xlabel('ABBAABBA set No. in Kyle Folder')

plt.figure()
plt.plot(np.arange(84), Deg_DEC)
plt.ylabel('Declination in Degrees')
plt.xlabel('ABBAABBA set No. in Kyle Folder')

plt.figure()
plt.plot(Deg_RA, Deg_DEC, 'bo')
plt.ylabel('Declination in Degrees')
plt.xlabel('Right Ascension in Degrees')