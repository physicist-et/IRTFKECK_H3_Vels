# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 02:07:18 2021

@author: snowy
"""

#Adjusting with shift the Q2 - > Q3,1 order for iSHELL data
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits #import the relevant directories to read a fits file from your directory and plot it
import h3ppy
from ACRE_tss import acre
from HenrikTheory import fit, wave

fit_IR = fit
wave_IR = wave

#%%
Wave_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.4.fits.gz'
Wave_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.5.fits.gz'
Wave_info6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.6.fits.gz'
Wave_info7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.7.fits.gz'

hdu = fits.open(Wave_info4)
hdr = hdu[0].header
#print(hdr)

Wavelength_O4 = fits.getdata(Wave_info4, ext=0)
Wavelength_O5 = fits.getdata(Wave_info5, ext=0)
Wavelength_O6 = fits.getdata(Wave_info6, ext=0)
Wavelength_O7 = fits.getdata(Wave_info7, ext=0)

#%%

image_file1 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.01.fits.gz' #imports a single fits file from the cwd + the file path in the commas
image_file2 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.02.fits.gz'
image_file3 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.03.fits.gz'
image_file4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.04.fits.gz'
image_file5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.05.fits.gz'
image_file6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.06.fits.gz'
image_file7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.07.fits.gz'
image_file8 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.08.fits.gz'
image_file9 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.014/uranus.014.order.09.fits.gz'
image_data1 = fits.getdata(image_file1, ext=0) #Focuses specifically on the array data contained within the fits file

image_data3 = fits.getdata(image_file3, ext=0)
image_data4 = fits.getdata(image_file4, ext=0) #Focuses specifically on the array data contained within the fits file
image_data5 = fits.getdata(image_file5, ext=0)
image_data6 = fits.getdata(image_file6, ext=0)
image_data7 = fits.getdata(image_file7, ext=0)
image_data_alt = np.concatenate((image_data3, image_data4, image_data5, image_data6, image_data7), axis=0)
#Single_IRTF_Data_Alt = acre(image_data_alt, width = 100, verbose = False)
image_data_alt[image_data_alt > 0.025] = np.nan
image_data_alt[image_data_alt < -0.025] = np.nan

plt.figure()
plt.imshow(np.fliplr(np.flipud(image_data_alt)), cmap='gist_gray')
plt.xlabel(r'Spectral pixel position (Pixel No.) ', fontsize=16)
#label_x = 0, 100, 200, 300, 400, 500, 600, 700
#plt.xticks(label_x, ('3.9506', '3.9564', '3.9622', '3.9680', '3.9738', '3.9796', '3.9854', '3.9919'), fontsize=15)
plt.ylabel('Spatial position across Slit (Pixel No.)', fontsize=16)
#plt.axvline(833, color='red', ls='-.', alpha=0.125)
plt.annotate(r'Q(1,0${^-}$)', (862, 395), color='red', fontsize=12)
#plt.axvline(688, color='blue', ls='-.', alpha=0.125)
plt.annotate(r'Q(3,0${^-}$)', (717, 285), color='blue', fontsize=12)
#plt.axhline(40, color='r', ls='--')
plt.title(r'Single AB of Uranus with ISHELL on $11^{th}$ October 2016', fontsize=23)
cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
cbar.set_label(r'Intensity ($Wm^{-2}sr^{-1}$)', fontsize=20)

#%%
image_data2 = fits.getdata(image_file2, ext=0)
image_data3 = fits.getdata(image_file3, ext=0)
image_data4 = fits.getdata(image_file4, ext=0)
image_data5 = fits.getdata(image_file5, ext=0)
image_data6 = fits.getdata(image_file6, ext=0)
image_data7 = fits.getdata(image_file7, ext=0)
image_data8 = fits.getdata(image_file8, ext=0)
image_data9 = fits.getdata(image_file9, ext=0)

#First create lists in which the arrays of Keck data will go into
IRTF_DataT = np.concatenate((image_data1, image_data2, image_data3, image_data4, image_data5, image_data6, image_data7, image_data8, image_data9))
Q_IRTF_Data = np.concatenate((image_data5, image_data6, image_data7), axis=0)

#plt.figure()
#plt.imshow(Q_IRTF_Data, cmap='gist_gray')

No_list = 14, 17, 18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41, 42, 45, 46, 49, 50, 53, 54, 57, 58, 61, 62, 65, 66, 69, 70, 73, 74, 77, 78, 81, 82, 85, 86, 89, 90, 93, 94, 97, 98, 101, 102
s = 0 #This will be used to create the ABBA pattern
All_IRTF_Data = []
for n in range(45): #We use this list to create a list which holds all the data from Order19
    if n < 43:
        for x in range(9):
            if x == 0:
                xx = x + 1
                image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.0' + str(No_list[n]) + '/uranus.0' + str(No_list[n]) + '.order.0' + str(xx) + '.fits.gz'
                image_datai = fits.getdata(image_filei, ext=0)
                IRTF_Data = image_datai
            else:
                xx = x + 1
                image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.0' + str(No_list[n]) + '/uranus.0' + str(No_list[n]) + '.order.0' + str(xx) + '.fits.gz'
                image_datai = fits.getdata(image_filei, ext=0)
                IRTF_Data = np.concatenate((IRTF_Data, image_datai))
        All_IRTF_Data.append(IRTF_Data)
    else:
        for x in range(9):
            if x == 0:
                xx = x + 1
                image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.' + str(No_list[n]) + '/uranus.' + str(No_list[n]) + '.order.0' + str(xx) + '.fits.gz'
                image_datai = fits.getdata(image_filei, ext=0)
                IRTF_Data = image_datai
            else:
                xx = x + 1
                image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus/uranus.' + str(No_list[n]) + '/uranus.' + str(No_list[n]) + '.order.0' + str(xx) + '.fits.gz'
                image_datai = fits.getdata(image_filei, ext=0)
                IRTF_Data = np.concatenate((IRTF_Data, image_datai))
        All_IRTF_Data.append(IRTF_Data)
        
#%% Now to call in the sky lines for Q2 - > Q3,1
from lmfit import Model
from types import SimpleNamespace

#First lets fit the Q1, Q2, Q3, Q3,1 and Q3,2 if possible (Only Q3,1 possible)
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

gmodel = Model(gauss_fit)
Sky_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.4.fits.gz'
Sky_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.5.fits.gz'
Sky_info6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.6.fits.gz'
Sky_info7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.7.fits.gz'

sky_data5 = fits.getdata(Sky_info5, ext=0)
rough_sky_data5 = sky_data5[8:85,:]
rough_sky_data5 = acre(rough_sky_data5, width = 10, verbose = False)
Sky_lines5 = np.zeros(2048)

for i in range(77):
    Sky_lines5 += rough_sky_data5[i,:]

plt.figure()
plt.imshow(sky_data5, cmap='gist_gray', vmin = 0, vmax = 10e4)

X = np.arange(60)
Mid_point1 = []
A0s = []

#Now to go through all 80 pixels and see if the medium position has changed for the first Sky line
for i in range(77):
    sky_line_O5 = rough_sky_data5[i,:]
    a2_guess = 2.5
    a1_guess = 31
    Data = sky_line_O5[1085:1145]
    resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = a2_guess, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
    pSky51 = SimpleNamespace(**resultSky5.best_values)
    Mid_point1.append(pSky51.a1)
    A0s.append(pSky51.a0)

#%% Order 2020 also showed odd reduction attempts so instead lets find out what the skew is and attempt to move the pixels into position
#To do that we need to measure the gaps and see how they consistantly move, it maybe worth doing a Gaussian? 

plt.figure()
plt.plot(sky_data5[50])

Mid_point2 = []
X = np.arange(30)
A0s = []

sky_line_O5 = rough_sky_data5[50,:]
dataiiloc = np.where(sky_line_O5 == np.max(sky_line_O5[1455:1485]))
dataiiloc = dataiiloc[0][0]
Data = sky_line_O5[1455:1485]
a1_guess = dataiiloc - 1455
resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = 1.5, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
pSky51 = SimpleNamespace(**resultSky5.best_values)
Mid_point2.append(pSky51.a1)
A0s.append(pSky51.a0)

#%%

Region1 = rough_sky_data5
RegionA = Region1[:, 1570:1630]
#RegionA = acre(RegionA, width = 3, verbose = False)
RegionB = Region1[:, 1440:1500]
#RegionB = acre(RegionB, width = 3, verbose = False)
RegionC = Region1[:, 1950:2010]
#RegionC = acre(RegionC, width = 3, verbose = False)
RegionD = Region1[:, 1890:1950]
#RegionD = acre(RegionD, width = 3, verbose = False)

Regions = RegionA + RegionB


plt.figure()
plt.imshow(RegionB, cmap='gist_gray')
plt.vlines(20, -0.5, 76.5, color='b')
plt.vlines(40, -0.5, 76.5, color='r')

plt.figure()
plt.imshow(RegionA, cmap='gist_gray')
plt.vlines(20, -0.5, 76.5, color='b')
plt.vlines(40, -0.5, 76.5, color='r')

plt.figure()
plt.imshow(RegionD, cmap='gist_gray')
plt.vlines(20, -0.5, 76.5, color='b')
plt.vlines(40, -0.5, 76.5, color='r')

plt.figure()
plt.imshow(RegionC, cmap='gist_gray')
plt.vlines(20, -0.5, 76.5, color='b')
plt.vlines(40, -0.5, 76.5, color='r')
#%% Now find 2/3 lines on the left side to show how Mid point varies

X = np.arange(60)
Mid_point2 = []

#Line 1 = 735
for i in range(77):
    sky_line_O5 = rough_sky_data5[i,:]
    dataiiloc = np.where(sky_line_O5 == np.max(sky_line_O5[705:765]))
    dataiiloc = dataiiloc[0][0] - 705
    a2_guess = 2.5
    a1_guess = dataiiloc
    Data = sky_line_O5[705:765]
    resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = a2_guess, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
    pSky51 = SimpleNamespace(**resultSky5.best_values)
    Mid_point2.append(pSky51.a1)

#Line 2 = 85
Mid_point3 = []

for i in range(77):
    sky_line_O5 = rough_sky_data5[i,:]
    dataiiloc = np.where(sky_line_O5 == np.max(sky_line_O5[55:115]))
    dataiiloc = dataiiloc[0][0] - 55
    a2_guess = 2.5
    a1_guess = dataiiloc
    a2_guess = 2.5
    a1_guess = 31
    Data = sky_line_O5[55:115]
    resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = a2_guess, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
    pSky51 = SimpleNamespace(**resultSky5.best_values)
    Mid_point3.append(pSky51.a1)

#Line 3 = 395
Mid_point4 = []

for i in range(77):
    sky_line_O5 = rough_sky_data5[i,:]
    dataiiloc = np.where(sky_line_O5 == np.max(sky_line_O5[365:425]))
    dataiiloc = dataiiloc[0][0] - 365
    a2_guess = 2.5
    a1_guess = dataiiloc
    a2_guess = 2.5
    a1_guess = 31
    Data = sky_line_O5[365:425]
    resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = a2_guess, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
    pSky51 = SimpleNamespace(**resultSky5.best_values)
    Mid_point4.append(pSky51.a1)

#And a fouth line  = 712
Mid_point5 = []

for i in range(77):
    sky_line_O5 = rough_sky_data5[i,:]
    dataiiloc = np.where(sky_line_O5 == np.max(sky_line_O5[682:742]))
    dataiiloc = dataiiloc[0][0] - 682
    a2_guess = 2.5
    a1_guess = dataiiloc
    a2_guess = 2.5
    a1_guess = 31
    Data = sky_line_O5[682:742]
    resultSky5 = gmodel.fit(Data, x=X, a0 = 100000, a1 = a1_guess, a2 = a2_guess, a3 = 2662504, a4 = 1992.35482, a5 = 46.3771111)
    pSky51 = SimpleNamespace(**resultSky5.best_values)
    Mid_point5.append(pSky51.a1)

#Now to see how they all plot out

plt.figure()
plt.plot(np.arange(77), Mid_point2, color='r')
plt.plot(np.arange(77), Mid_point3, color='b')
plt.plot(np.arange(77), Mid_point4, color='g')
plt.plot(np.arange(77), Mid_point5, color='m')

#Now to find the relationship in these lines to find out if they are 
#%%
plt.figure()
plt.plot(Mid_point3, np.arange(77), color='b')
plt.xlabel('Middle Point of the Gaussian from Pixel No. 55', fontsize=20)
plt.ylabel('Slit No. Pixel Position (+8)', fontsize=20)

plt.figure()
plt.plot(Mid_point4, np.arange(77), color='g')
plt.xlabel('Middle Point of the Gaussian from Pixel No. 365', fontsize=20)
plt.ylabel('Slit No. Pixel Position (+8)', fontsize=20)

plt.figure()
plt.plot(Mid_point5, np.arange(77), color='m')
plt.xlabel('Middle Point of the Gaussian from Pixel No. 682', fontsize=20)
plt.ylabel('Slit No. Pixel Position (+8)', fontsize=20)

plt.figure()
plt.plot(Mid_point2, np.arange(77), color='r')
plt.xlabel('Middle Point of the Gaussian from Pixel No. 705', fontsize=20)
plt.ylabel('Slit No. Pixel Position (+8)', fontsize=20)
