# -*- coding: utf-8 -*-
"""
Created on Wed Mar 24 16:21:40 2021

@author: Emma 'Amelia'
"""
#Attempt to pull the correct IRTF data to find Q lines

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits #import the relevant directories to read a fits file from your directory and plot it
import h3ppy
from ACRE_tss import acre
# from HenrikTheory import fit, wave
from numpy import inf
from types import SimpleNamespace
import math
import scipy

# fit_IR = fit
# wave_IR = wave

#%%
Wave_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.4.fits.gz'
Wave_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.5.fits.gz'
Wave_info6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.6.fits.gz'
Wave_info7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.7.fits.gz'

hdu = fits.open(Wave_info4)
hdr = hdu[0].header
#print(hdr)

Wavelength_O4 = np.load('Wavelength_O4.npy')
Wavelength_O5 = np.load('Wavelength_O5.npy')
Wavelength_O6 = fits.getdata(Wave_info6, ext=0)
Wavelength_O7 = fits.getdata(Wave_info7, ext=0)

#%%

image_file1 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.01.fits.gz' #imports a single fits file from the cwd + the file path in the commas
image_file2 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.02.fits.gz'
image_file3 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.03.fits.gz'
image_file4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.04.fits.gz'
image_file5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.05.fits.gz'
image_file6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.06.fits.gz'
image_file7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.07.fits.gz'
image_file8 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.08.fits.gz'
image_file9 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.014/uranus.014.order.09.fits.gz'
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

#To shift data correctly we add in the popt values by hand
A_average = 0
B_average = 0
TIMES = []

No_list = 14, 17, 18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41, 42, 45, 46, 49, 50, 53, 54, 57, 58, 61, 62, 65, 66, 69, 70, 73, 74, 77, 78, 81, 82, 85, 86, 89, 90, 93, 94, 97, 98, 101, 102
s = 0 #This will be used to create the ABBA pattern
All_IRTF_Data_Q1 = []
for n in range(45): #We use this list to create a list which holds all the data from Order19
    if n < 43:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.0' + str(No_list[n]) + '/uranus.0' + str(No_list[n]) + '.order.04.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        hdu = fits.open(image_filei)
        hdr = hdu[0].header
        TIMES.append(hdr['TIME_OBS'])
        IRTF_Data = image_datai
        IRTF_Data[IRTF_Data == -inf] = 0
        IRTF_Data[IRTF_Data == inf] = 0
        IRTF_Data[np.isnan(IRTF_Data)] = 0
        All_IRTF_Data_Q1.append(IRTF_Data)
    else:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.' + str(No_list[n]) + '/uranus.' + str(No_list[n]) + '.order.04.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        IRTF_Data = image_datai
        hdu = fits.open(image_filei)
        hdr = hdu[0].header
        TIMES.append(hdr['TIME_OBS'])
        IRTF_Data[IRTF_Data == -inf] = 0
        IRTF_Data[IRTF_Data == inf] = 0
        IRTF_Data[np.isnan(IRTF_Data)] = 0
        All_IRTF_Data_Q1.append(IRTF_Data)

All_IRTF_Data_Q3 = []
for n in range(45): #We use this list to create a list which holds all the data from Order19
    if n < 43:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.0' + str(No_list[n]) + '/uranus.0' + str(No_list[n]) + '.order.05.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        IRTF_Data = image_datai
        IRTF_Data[IRTF_Data == -inf] = 0
        IRTF_Data[IRTF_Data == inf] = 0
        IRTF_Data[np.isnan(IRTF_Data)] = 0
        All_IRTF_Data_Q3.append(IRTF_Data)
    else:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.' + str(No_list[n]) + '/uranus.' + str(No_list[n]) + '.order.05.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        IRTF_Data = image_datai
        IRTF_Data[IRTF_Data == -inf] = 0
        IRTF_Data[IRTF_Data == inf] = 0
        IRTF_Data[np.isnan(IRTF_Data)] = 0
        All_IRTF_Data_Q3.append(IRTF_Data)    

#%% Now we need to bin the data into appropriate boxes
IRTF_Data_Q1 = []
Approx_Mid = (89-7)/2 + 7
S_point = All_IRTF_Data_Q1[0]
import scipy.ndimage as f #This will allow the equivalent of fshift-ing in IDL
corrective_x = np.load('Corrective_Array_O4.npy')

for n in range(45):
    if n == 0:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q1[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = corrective_x[i][1216] #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q1 = shift_Spec_row
    else:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q1[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = corrective_x[i][1216] #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q1 = np.dstack([IRTF_Data_Q1, shift_Spec_row])
        
for n in range(45):
    if n == 0:
        Total_IRTF_Data = IRTF_Data_Q1[:,:,0]
    else:
        Total_IRTF_Data = np.add(Total_IRTF_Data, IRTF_Data_Q1[:,:,n])
        
# plt.figure()
# plt.imshow(Total_IRTF_Data, cmap='gist_gray', vmax = 100, vmin = -100)
# plt.colorbar()
        
corrective_x = np.load('Corrective_Array_O5.npy')
IRTF_Data_Q3 = []
Approx_Mid = (90-6)/2 + 6
S_point = All_IRTF_Data_Q3[0]

for n in range(45):
    if n == 0:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q3[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = corrective_x[i][1366] #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q3 = shift_Spec_row
    else:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q3[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = corrective_x[i][1366]  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q3 = np.dstack((IRTF_Data_Q3, shift_Spec_row))
        
for n in range(45):
    if n == 0:
        Total_IRTF_DataQ3 = IRTF_Data_Q3[:,:,0]
    else:
        Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, IRTF_Data_Q3[:,:,n])
        
IRTF_Data_Q31 = []
Approx_Mid = (90-6)/2 + 6
S_point = All_IRTF_Data_Q3[0]

for n in range(45):
    if n == 0:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q3[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][1445:1465]) #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q31 = shift_Spec_row
    else:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q3[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][1445:1465])  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q31 = np.dstack((IRTF_Data_Q31, shift_Spec_row))
        
for n in range(45):
    if n == 0:
        Total_IRTF_DataQ31 = IRTF_Data_Q31[:,:,0]
    else:
        Total_IRTF_DataQ31 = np.add(Total_IRTF_DataQ31, IRTF_Data_Q31[:,:,n])
        
# plt.figure()
# plt.imshow(Total_IRTF_DataQ3, cmap='gist_gray', vmax = 2*10**-1, vmin=-2*10**-1) #NEED TO DO A Q3,1 corrected spectra!!!! #CHECK BOTH REDUCTION RESIDUALS TO SEE IF THE OFFSET MIGHT EXPLAIN THE FINAL OFFSET ISSUE!
        
#%%
h3p = h3ppy.h3p(line_list_file = 'C:/Users/snowy/OneDrive/Documents/Python work/h3p_line_list_neale_1996_subset.txt')
wave2 = h3p.wavegen(np.nanmin(Wavelength_O4), np.nanmax(Wavelength_O4), 2048)
Total_IRTF_Data_Alt = Total_IRTF_Data/90
Total_IRTF_Data_Alt = acre(Total_IRTF_Data_Alt, width = 15, verbose = False)
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt > 0.025] = 0
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt < -0.025] = 0

plt.figure()
plt.imshow(Total_IRTF_Data, cmap='gist_gray', vmax = 0.25, vmin = -0.25)
plt.xlabel(r'Wavelength ' + '($\mu$m)', fontsize = 15)
#label_x = 0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000
#plt.xticks(label_x, (str(wave2[0]), str(wave2[200]), str(wave2[400]), str(wave2[600]), str(wave2[800]), str(wave2[1000]), str(wave2[1200]), str(wave2[1400]), str(wave2[1600]), str(wave2[1800]), str(wave2[2000])), fontsize=10)
plt.ylabel('Spatial position across Slit (Pixel No.)', fontsize=16)
# =============================================================================
# plt.axvline(1214, color='red', ls='-.', alpha=0.125)
# plt.axvline(1359, color='blue', ls='-.', alpha=0.125)
# plt.axvline(1456, color='green', ls='-.', alpha=0.125)
# plt.annotate(r'Q(1,0${^-}$)', (1225, 620), color='red', fontsize=12)
# plt.annotate(r'Q(3,0${^-}$)', (1371, 530), color='blue', fontsize=12)
# plt.annotate(r'Q(3,1${^-}$)', (1476, 530), color='green', fontsize=12)
# =============================================================================
#plt.axhline(40, color='r', ls='--')
plt.title(r'Total 5.5hr observations of Uranus with ISHELL on $11^{th}$ October 2016', fontsize=23)
#plt.Circle(, radius=12.049815498154981, color ='r', lw=2, ls='--', fill=False)
cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
cbar.set_label(r'Intensity ($Wm^{-2}sr^{-1}$)', fontsize=20)

wave3 = h3p.wavegen(np.nanmin(Wavelength_O5), np.nanmax(Wavelength_O5), 2048)
Total_IRTF_Data_AltQ3 = Total_IRTF_DataQ3/45
Total_IRTF_Data_Alt = acre(Total_IRTF_Data_Alt, width = 15, verbose = False)
Total_IRTF_Data_AltQ3 = acre(Total_IRTF_Data_AltQ3, width = 15, verbose = False)
# Total_IRTF_Data_AltQ3[Total_IRTF_Data_AltQ3 > 0.025] = 0
# Total_IRTF_Data_AltQ3[Total_IRTF_Data_AltQ3 < -0.025] = 0

plt.figure()
plt.imshow(Total_IRTF_DataQ3, cmap='gist_gray', vmax = 0.25, vmin = -0.25)
plt.xlabel(r'Wavelength ' + '($\mu$m)', fontsize = 15)
#label_x = 0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000
#plt.xticks(label_x, (str(wave2[0]), str(wave2[200]), str(wave2[400]), str(wave2[600]), str(wave2[800]), str(wave2[1000]), str(wave2[1200]), str(wave2[1400]), str(wave2[1600]), str(wave2[1800]), str(wave2[2000])), fontsize=10)
plt.ylabel('Spatial position across Slit (Pixel No.)', fontsize=16)
# =============================================================================
# plt.axvline(1214, color='red', ls='-.', alpha=0.125)
# plt.axvline(1359, color='blue', ls='-.', alpha=0.125)
# plt.axvline(1456, color='green', ls='-.', alpha=0.125)
# plt.annotate(r'Q(1,0${^-}$)', (1225, 620), color='red', fontsize=12)
# plt.annotate(r'Q(3,0${^-}$)', (1371, 530), color='blue', fontsize=12)
# plt.annotate(r'Q(3,1${^-}$)', (1476, 530), color='green', fontsize=12)
# =============================================================================
#plt.axhline(40, color='r', ls='--')
plt.title(r'Total 5.5hr observations of Uranus with ISHELL on $11^{th}$ October 2016', fontsize=23)
cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
cbar.set_label(r'Intensity ($Wm^{-2}sr^{-1}$)', fontsize=20)

#%% There's a slight misalignment so we're going to move the positive data up and along to match up 
# First lets find the negative and positive line middle positions and then include this into the translation
from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Pos_Q1Tn = []
Pos_Q1Tp = []
Err_Q1Tn = []
Err_Q1Tp = []

#For Q1 Neg (Need to clean up a little)
for m in range(27):
    dataI = Total_IRTF_Data
    datai = dataI[13+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmin(datai[1205:1225]))
    a1_pointP = a1_pointP[0][0] - 1110
    ex_resultQ1P = gmodel.fit(datai[1110:1310], x=np.arange(200), a0=np.nanmin(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q1Tn.append(ex_pQ1P.a1 + 1110)
    Err_Q1Tn.append(ex_eQ1P[1])
    datai = dataI[54+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    a1_pointP = a1_pointP[0][0] - 1110
    ex_resultQ1P = gmodel.fit(datai[1110:1310], x=np.arange(200), a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q1Tp.append(ex_pQ1P.a1 + 1110)
    Err_Q1Tp.append(ex_eQ1P[1])

Olap = []
Err_lap = []

for m in range(27):
    Olap.append(Pos_Q1Tp[m] - Pos_Q1Tn[m])
    Err_lap.append(np.sqrt((Err_Q1Tp[m])**2+(Err_Q1Tn[m])**2))
    
Pos_Q3Tn = []
Pos_Q3Tp = []
Err_Q3Tn = []
Err_Q3Tp = []

#For Q3 Neg (Need to clean up a little)
for m in range(25):       
    dataI = Total_IRTF_DataQ3
    datai = dataI[14+m, :]
    datai[datai > 0] = 0
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmin(datai[1350:1370]))
    a1_pointP = a1_pointP[0][0] - 1250
    ex_resultQ1P = gmodel.fit(datai[1250:1450], x=np.arange(200), a0=np.nanmin(datai[1350:1370]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q3Tn.append(ex_pQ1P.a1 + 1250)
    Err_Q3Tn.append(ex_eQ1P[1])
    dataI = Total_IRTF_DataQ3
    datai = dataI[55+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1350:1370]))
    a1_pointP = a1_pointP[0][0] - 1250
    ex_resultQ1P = gmodel.fit(datai[1250:1450], x=np.arange(200), a0=np.nanmax(datai[1350:1370]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q3Tp.append(ex_pQ1P.a1 + 1250)
    Err_Q3Tp.append(ex_eQ1P[1])

Olapv3 = []
Err_lapv3 = []

for m in range(25):
    Olapv3.append(Pos_Q3Tp[m] - Pos_Q3Tn[m])
    Err_lapv3.append(np.sqrt((Err_Q3Tp[m])**2+(Err_Q3Tn[m])**2))
    
#%% For Q3,1
Pos_Q31Tn = []
Pos_Q31Tp = []
Err_Q31Tn = []
Err_Q31Tp = []

#For Q3 Neg (Need to clean up a little)
for m in range(25):
    dataI = Total_IRTF_DataQ31
    datai = dataI[15+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmin(datai[1445:1465]))
    a1_pointP = a1_pointP[0][0] - 1350
    ex_resultQ1P = gmodel.fit(datai[1350:1550], x=np.arange(200), a0=np.nanmin(datai[1445:1465]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q31Tn.append(ex_pQ1P.a1 + 1350)
    Err_Q31Tn.append(ex_eQ1P[1])
    datai = dataI[56+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1445:1465]))
    a1_pointP = a1_pointP[0][0] - 1350
    ex_resultQ1P = gmodel.fit(datai[1350:1550], x=np.arange(200), a0=np.nanmax(datai[1445:1465]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q31Tp.append(ex_pQ1P.a1 + 1350)
    Err_Q31Tp.append(ex_eQ1P[1])

Olapv31 = []
Err_lapv31 = []

for m in range(25):
    Olapv31.append(Pos_Q31Tp[m] - Pos_Q31Tn[m])
    Err_lapv31.append(np.sqrt((Err_Q31Tp[m])**2+(Err_Q31Tn[m])**2))

#%% Now plot a single line of data and make a line graph
#Total_IRTF_Data = acre(Total_IRTF_Data[330:551,:], width = 5, verbose = True)
T_IRTFData = Total_IRTF_Data_Alt
h3p = h3ppy.h3p(line_list_file = 'C:/Users/snowy/OneDrive/Documents/Python work/h3p_line_list_neale_1996_subset.txt')
wave2 = Wavelength_O4

#%% Here use the advantage if the star seperation to help here
T_IRTFData = Total_IRTF_Data_Alt

Order_4_Offset = np.ones(2048)*41.2337 #These are in microns
Order_5_Offset = np.ones(2048)*41.6508

#Lets start shifting the positive line onto the negative line looping through all of x
Shape_Order4 = np.shape(IRTF_Data_Q1[:,:,0])
Combined_Q1 = []

for x in range(45):
    IRTFData = IRTF_Data_Q1[:,:,x]
    if x == 35:
        IRTFData = acre(IRTFData, width = 15, verbose = False)
    for n in range(Shape_Order4[1]):
        Spec_col = IRTFData[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[Order_4_Offset[n]], mode='wrap')
        if n == 0:
            Order_4_shift = shifted_Spec_col
        else:
            Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))
    Order_4_shiftv1 = np.flipud(np.rot90(Order_4_shift))
    #plt.imshow(Order_4_shiftv1, cmap='gist_gray', vmax = 0.025, vmin = -0.025)
    for n in range(Shape_Order4[0]):
        Spec_col = Order_4_shiftv1[n,:]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[np.nanmean(Olap)], mode='wrap')
        if n == 0:
            Order_4_shift = shifted_Spec_col
        else:
            Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))
    Order_4_shift = Order_4_shift
    Combined_Q1.append(Order_4_shift)
    # plt.imshow(Combined_Q1[x], cmap='gist_gray', vmax = 0.025, vmin = -0.025)
#%%
#Lets start shifting the positive line onto the negative line looping through all of x
Shape_Order5 = np.shape(IRTF_Data_Q3[:,:,0])
Combined_Q3 = []

for x in range(45):
    IRTFData = IRTF_Data_Q3[:,:,x]
    #IRTFData = acre(IRTFData, width = 15, verbose = False)
    for n in range(Shape_Order5[1]):
        Spec_col = IRTFData[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[Order_5_Offset[n]], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shiftv1 = np.flipud(np.rot90(Order_5_shift))
    for n in range(Shape_Order5[0]):
        Spec_col = Order_5_shiftv1[n,:]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[np.nanmean(Olapv3)], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shift = Order_5_shift
    Combined_Q3.append(Order_5_shift)

#%%
Shape_Order5A = np.shape(IRTF_Data_Q31[:,:,0])
Combined_Q31 = []

for x in range(45):
    IRTFData = IRTF_Data_Q31[:,:,x]
    #IRTFData = acre(IRTFData, width = 15, verbose = False)
    for n in range(Shape_Order5[1]):
        Spec_col = IRTFData[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[Order_5_Offset[n]], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shiftv1 = np.flipud(np.rot90(Order_5_shift))
    for n in range(Shape_Order5[0]):
        Spec_col = Order_5_shiftv1[n,:]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[np.nanmean(Olapv31)], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shift = Order_5_shift
    Combined_Q31.append(Order_5_shift)

#%%
import matplotlib.ticker as tkr
AB_Combined_Order4 = []
AB_Combined_Order5 = []
AB_Combined_Order5B = []

for x in range(45):
    IRTFData_Q1 = -1*Combined_Q1[x] + IRTF_Data_Q1[:,:,x]
    # plt.figure()
    # plt.imshow(IRTFData_Q1, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    AB_Combined_Order4.append(IRTFData_Q1)

for x in range(45):
    IRTFData_Q3 = -1*Combined_Q3[x] + IRTF_Data_Q3[:,:,x]
    AB_Combined_Order5.append(IRTFData_Q3)
    
for x in range(45):
    IRTFData_Q31 = -1*Combined_Q31[x] + IRTF_Data_Q31[:,:,x]
    AB_Combined_Order5B.append(IRTFData_Q31)    

Q1_Mean = np.nanmean(AB_Combined_Order4[16:24], axis = 0)
Q3_Mean = np.nanmean(AB_Combined_Order5[16:24], axis = 0)

# fig, (ax1, ax2) = plt.subplots(1, 2)
# fig.suptitle('a) Averaged $H_{3}^{+}$ emission spectrum for Keck II NIRSPEC in October 2014', fontsize = 30)
# ax1.imshow(Q1_Mean[40:100, 1116:1316], cmap='gist_heat', vmax = 0.02, vmin = -0.02)
# ax2.imshow(Q3_Mean[:, 1260:1555], cmap='gist_heat', vmax = 0.02, vmin = -0.02)

plt.figure()
plt.imshow(np.flipud(Q1_Mean[50:90, 716:1716]), cmap='gist_heat', vmax = 0.02, vmin = 0)
plt.title('d) Averaged $H_{3}^{+}$ emission spectrum for IRTF iSHELL in October 2016', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.3f'))
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Intensity ${\mu}Wm^{-2}sr^{-1}$', fontsize = 25)
plt.xlabel('Wavelength ($\mu$m)', fontsize = 25)
plt.ylabel('Spatial axis (Pixels)', fontsize = 25)
plt.yticks(fontsize = 20)
plt.xticks(ticks = (0, 999), labels=("{:.3f}".format(round(Wavelength_O4[716], 2)), "{:.3f}".format(round(Wavelength_O4[1715], 2))), fontsize = 20)
plt.vlines((488, 507), ymin = 0, ymax = 60, color='w', lw = 3, ls = 'dashed')
plt.ylim(0, 40)
plt.text(469, 4, '$Q(1,0^{-})$', fontsize = 13, color = 'w', rotation = 90)
# plt.text(345.5, 10, '$Q(2,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
# plt.text(594.5, 10, '$Q(3,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
# plt.text(754.5, 10, '$Q(3,2^{-})$', fontsize = 15, color = 'w', rotation = 90)
# plt.text(620.5, 10, '$Q(3,1^{-})$', fontsize = 15, color = 'w', rotation = 90)

plt.figure()
plt.imshow(np.flipud(Q3_Mean[50:90, 908:1908]), cmap='gist_heat', vmax = 0.02, vmin = 0)
plt.title('e) Averaged $H_{3}^{+}$ emission spectrum for IRTF iSHELL in October 2016', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.3f'))
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Intensity ${\mu}Wm^{-2}sr^{-1}$', fontsize = 25)
plt.xlabel('Wavelength ($\mu$m)', fontsize = 25)
plt.ylabel('Spatial axis (pixels)', fontsize = 25)
plt.yticks(fontsize = 20)
plt.xticks(ticks = (0, 999), labels=("{:.3f}".format(round(Wavelength_O5[908], 2)), "{:.3f}".format(round(Wavelength_O5[1907], 2))), fontsize = 20)
plt.vlines((439, 464), ymin = 0, ymax = 60, color='w', lw = 3, ls = 'dashed')
plt.ylim(0, 40)
plt.text(421, 4, '$Q(3,0^{-})$', fontsize = 13, color = 'w', rotation = 90)
plt.text(521, 4, '$Q(3,1^{-})$', fontsize = 13, color = 'w', rotation = 90)
# plt.text(754.5, 10, '$Q(3,2^{-})$', fontsize = 15, color = 'w', rotation = 90)
# plt.text(620.5, 10, '$Q(3,1^{-})$', fontsize = 15, color = 'w', rotation = 90)

plt.figure()
plt.imshow(np.flipud(Q3_Mean[50:90, 0:1000])*1000, cmap='gist_heat', vmax = 20, vmin = 0)
plt.title('e) Averaged $H_{3}^{+}$ emission spectrum for IRTF iSHELL in October 2016', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.1f'), orientation='horizontal')
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_xlabel('Intensity $nWm^{-2}sr^{-1}$', fontsize = 25)
cbar.ax.xaxis.set_label_position('bottom')
#cbar.ax.set_xticks(fontsize = 20)
plt.xlabel('Wavelength ($\mu$m)', fontsize = 25)
plt.ylabel('Spatial axis (pixels)', fontsize = 25)
plt.yticks(fontsize = 20)
plt.xticks(ticks = (0, 999), labels=("{:.1f}".format(round(Wavelength_O5[908], 2)), "{:.1f}".format(round(Wavelength_O5[1907], 2))), fontsize = 20)
#plt.vlines((439, 464), ymin = 0, ymax = 60, color='w', lw = 3, ls = 'dashed')
plt.ylim(0, 40)
plt.text(442, 4, '$Q(2,0^{-})$', fontsize = 13, color = 'w', rotation = 90)
#plt.text(521, 7, '$Q(3,1^{-})$', fontsize = 15, color = 'w', rotation = 90)

#%% Lets attempt to fit H3+ onto this (which didn' work so lets just focus Q1 and Q3 seperately so we use the Keck example from h3ppy)
from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

ex_A0_Q1 = []
ex_Err_A0_Q1 = []
A0_Q3 = []
ex_A0_Q3 = []
Err_A0_Q3 = []
ex_A2_Q1 = []
ex_Err_A2_Q1 = []
A2_Q3 = []
ex_A2_Q3 = []
Err_A2_Q3 = []
ex_FWHMQ1 = []
FWHMQ3 = []
ex_FWHMQ3 = []
ex_INTSQ1 = []
INTSQ3 = []
ex_INTSQ3 = []
ex_Err_A0_Q3 = []
ex_Err_A2_Q3 = []
ex_A1_Q1 = []
ex_Err_A1_Q1 = []
ex_A1_Q3 = []
ex_Err_A1_Q3 = []
Pos_Q1 = []
Pos_Q3 = []

#Lets make a special case
ex_A0_Q18 = []
ex_Err_A0_Q18 = []
A0_Q38 = []
ex_A0_Q38 = []
Err_A0_Q38 = []
ex_A2_Q18 = []
ex_Err_A2_Q18 = []
A2_Q38 = []
ex_A2_Q38 = []
Err_A2_Q38 = []
ex_FWHMQ18 = []
FWHMQ38 = []
ex_FWHMQ38 = []
ex_INTSQ18 = []
INTSQ38 = []
ex_INTSQ38 = []
ex_Err_A0_Q38 = []
ex_Err_A2_Q38 = []
ex_A1_Q18 = []
ex_Err_A1_Q18 = []
ex_A1_Q38 = []
ex_Err_A1_Q38 = []

XX = np.arange(200) # the wavelengths either side of Q1 are 3.9529269748606297 and 3.953021233729671 over 3 pixels either side of Q1
o = 0
Q1_Data_Capture = []
Q3_Data_Capture = []

for x in range(11):
    dataI = (AB_Combined_Order4[o]+AB_Combined_Order4[o+1]+AB_Combined_Order4[o+2]+AB_Combined_Order4[o+3])/8
    #dataI[dataI > 0.020] = 0.0025
    Q1_Data_Capture.append(dataI[49:89, :])
    plt.figure()
    plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    o += 4
    for m in range(22):
        datai = dataI[57+m, :]
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
        ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
        A0 = ex_pQ1P.a0
        ex_A0_Q1.append(A0)
        A2 = ex_pQ1P.a2
        ex_A2_Q1.append((A2)*1.6199707031249897e-05)
        err_A0 = ex_eQ1P[0]
        ex_Err_A0_Q1.append(err_A0)
        err_A2 = ex_eQ1P[2]
        ex_Err_A2_Q1.append(err_A2*1.6199707031249897e-05)
        ex_A1_Q1.append(ex_pQ1P.a1)
        ex_Err_A1_Q1.append(ex_eQ1P[1]*1.6199707031249897e-05)
        ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQ1.append(ex_FWHM)
        ex_INTSQ1.append(A0*ex_FWHM*(10**6))
        #Pos_Q1.append(ex_pQ1P.a1 + 1110)
        # if x == 9:
        #     plt.figure()
        #     plt.plot(XX, datai[1110:1310], 'bo')   
        #     plt.plot(XX, ex_resultQ1P.best_fit, '--')
        #     print(ex_pQ1P.a1)
      
Q1_IntErr = []
# #Slit is 15 so

for o in range(242):
    Q1IntErr = ex_INTSQ1[o]*np.sqrt((ex_Err_A0_Q1[o]/ex_A0_Q1[o])**2 + (ex_Err_A2_Q1[o]/ex_A2_Q1[o])**2)
    Q1_IntErr.append(Q1IntErr)

Q31_Data_Capture = []
o = 0
for x in range(11):
    dataI = (AB_Combined_Order5[o]+AB_Combined_Order5[o+1]+AB_Combined_Order5[o+2]+AB_Combined_Order5[o+3])/8
    dataI[dataI < -0.004] = 0
    #dataI[dataI > 0.020] = 0.0025
    Q3_Data_Capture.append(dataI[49:89, :])
    Datai = (AB_Combined_Order5B[o]+AB_Combined_Order5B[o+1]+AB_Combined_Order5B[o+2]+AB_Combined_Order5B[o+3])/8
    Q31_Data_Capture.append(Datai[49:89, :])
    o += 4
    for m in range(22):
        datai = dataI[57+m, :]
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1350:1370]))
        a1_pointP = a1_pointP[0][0] - 1260
        ex_resultQ3P = gmodel.fit(datai[1260:1460], x=XX, a0=np.nanmax(datai[1350:1370]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQ3P = SimpleNamespace(**ex_resultQ3P.best_values)
        ex_eQ3P = np.sqrt(np.diag(ex_resultQ3P.covar))
        A0 = ex_pQ3P.a0
        ex_A0_Q3.append(A0)
        A2 = ex_pQ3P.a2
        ex_A2_Q3.append((A2)*1.6199707031249897e-05)
        err_A0 = ex_eQ3P[0]
        ex_Err_A0_Q3.append(err_A0)
        err_A2 = ex_eQ3P[2]
        ex_Err_A2_Q3.append(err_A2*1.6199707031249897e-05)
        ex_A1_Q3.append(ex_pQ3P.a1)
        ex_Err_A1_Q3.append(ex_eQ3P[1]*1.6199707031249897e-05)
        ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQ3.append(ex_FWHM)
        ex_INTSQ3.append(A0*ex_FWHM*(10**6)) #Needs work with Q3 to skip difficult pixels?
        #Pos_Q3.append(ex_pQ3P.a1 + 1260)
        # if x == 3:
        #     if m == 21:
        #         plt.figure()
        #         plt.plot(XX, datai[1260:1460], 'bo')
                # z = (XX - ex_pQ3P.a1) / ex_pQ3P.a2
                # y = ex_pQ3P.a0 * np.exp(-z**2 / ex_pQ3P.a2) + ex_pQ3P.a3 + ex_pQ3P.a4 * XX + ex_pQ3P.a5 * XX**2
                # plt.plot(XX, y, 'g')
                # print(A0*ex_FWHM*(10**6))       
      
Q3_IntErr = []
# #Slit is 15 so

for o in range(242):
    Q3IntErr = ex_INTSQ3[o]*np.sqrt((ex_Err_A0_Q3[o]/ex_A0_Q3[o])**2 + (ex_Err_A2_Q3[o]/ex_A2_Q3[o])**2)
    Q3_IntErr.append(Q3IntErr)

#%% Now we find Q1 and Q3 values across the slits Double up the sets and pixels (Start rewrite here!!)
           
#Now to plot all eleven lines
Mapping_Q1A = np.reshape(ex_INTSQ1, (11, 22))
Mapping_Q1A = Mapping_Q1A.transpose()
Mapping_Q1A[Mapping_Q1A == 0] = np.nan

Mapping_Q3A = np.reshape(ex_INTSQ3, (11, 22))
Mapping_Q3A = Mapping_Q3A.transpose()
Mapping_Q3A[Mapping_Q3A == 0] = np.nan

#%%
TIME = ['07:45', '08:14', '08:44', '09:14', '09:44', '10:13', '10:43', '11:12', '11:42', '12:11', '12:42']
TIMEv2 = ['12:42', '12:11', '11:42', '11:12', '10:43', '10:13', '09:44', '09:14', '08:44', '08:14', '07:45']
#PIXELS = ['-']

Mapping_Q1A = Mapping_Q1A.transpose()
Mapping_Q1B = np.flipud(Mapping_Q1A)
#Mapping_Q1B[Mapping_Q1B > 2.0] = 0
#Mapping_Q1B = np.where(Mapping_Q1A==Mapping_Q1A[8,0], np.nanmean(Mapping_Q1A[4:8,0]), Mapping_Q1A)
plt.figure()
plt.imshow(np.fliplr(Mapping_Q1B), cmap='nipy_spectral')
cbar = plt.colorbar()
cbar.set_label(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) +- ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr[75:154]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=17.5) 
plt.title(r'$H_{3}^{+}$ $Q(1,0^{-})$ Intensity of Uranus observations from 07:30 to 12:57 UTC on $11^{th}$ October 2016', pad = 45, fontsize = 20)
plt.xlabel('Arbitrary Longitude across Uranus (Pixels)', fontsize=17.5, labelpad=10)
plt.ylabel('UTC Time (HH:MM +- 00:15)', fontsize=17.5, labelpad=10)
plt.yticks(np.arange(11), TIMEv2)

Mapping_Q3A = Mapping_Q3A.transpose()
Mapping_Q3B = np.flipud(Mapping_Q3A)
plt.figure()
plt.imshow(np.fliplr(Mapping_Q3B), cmap='nipy_spectral')
cbar = plt.colorbar()
cbar.set_label(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) +- ' + str("{:.2f}".format(round(np.nanmean(Q3_IntErr[75:154]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=17.5) 
plt.title(r'$H_{3}^{+}$ $Q(3,0^{-})$ Intensity of Uranus observations from 07:30 to 12:57 UTC on $11^{th}$ October 2016', pad = 45, fontsize = 20)
plt.xlabel('Arbitrary Longitude across Uranus (Pixels)', fontsize=17.5, labelpad=10)
plt.ylabel('UTC Time (HH:MM +- 00:15)', fontsize=17.5, labelpad=10)
plt.yticks(np.arange(11), TIMEv2)

Mapping_Q1Err = np.reshape(Q1_IntErr, (11, 22))
Mapping_Q1Err = np.flipud(Mapping_Q1Err)
plt.figure()
plt.imshow(np.fliplr(Mapping_Q1Err), cmap='nipy_spectral')
cbar = plt.colorbar()
cbar.set_label(r'Intensity Error (${\mu}Wm^{-2}sr^{-1}$)', fontsize=17.5) 
plt.title(r'$H_{3}^{+}$ $Q(1,0^{-})$ Intensity Error of Uranus observations from 07:30 to 12:57 UTC on $11^{th}$ October 2016', pad = 45, fontsize = 20)
plt.xlabel('Arbitrary Longitude across Uranus (Pixels)', fontsize=17.5, labelpad=10)
plt.ylabel('UTC Time (HH:MM +- 00:15)', fontsize=17.5, labelpad=10)
plt.yticks(np.arange(11), TIMEv2)

Mapping_Q3Err = np.reshape(Q3_IntErr, (11, 22))
Mapping_Q3Err = np.flipud(Mapping_Q3Err)
plt.figure()
plt.imshow(np.fliplr(Mapping_Q3Err), cmap='nipy_spectral')
cbar = plt.colorbar()
cbar.set_label(r'Intensity Error (${\mu}Wm^{-2}sr^{-1}$)', fontsize=17.5) 
plt.title(r'$H_{3}^{+}$ $Q(3,0^{-})$ Intensity Error of Uranus observations from 07:30 to 12:57 UTC on $11^{th}$ October 2016', pad = 45, fontsize = 20)
plt.xlabel('Arbitrary Longitude across Uranus (Pixels)', fontsize=17.5, labelpad=10)
plt.ylabel('UTC Time (HH:MM +- 00:15)', fontsize=17.5, labelpad=10)
plt.yticks(np.arange(11), TIMEv2)

#%% Now lets approximate the Temperature from Q1 and Q3 (assuming Q3 is reduced correctly)
h = 6.63*(10**-34)
c = 3.00*(10**8)
kb = 1.38*(10**-23)
Q1Int = np.fliplr(Mapping_Q1B[:,0:23])
Q3Int = np.fliplr(Mapping_Q3B[:,0:23])
EQ1 = 2552.57
EQ3 = 2961.84
AQ1 = 128.7
AQ3 = 123.2
gQ1 = 4
gQ3 = 4
JQ1 = 3
JQ3 = 9
wQ1 = 2529.73
wQ3 = 2509.08

#Now find gamma to then find Temperature
gamma = (gQ1 * ((2 * JQ1) + 1)* h * c * wQ1 * AQ1)/(gQ3 *((2 * JQ3) +1)* h* c * wQ3 * AQ3)

Temps = []
#Now find Temperature and map it across, compare this to IDL to see what similarities fall out
for i in range(11):
    Rotat_Temp = ((EQ1 - EQ3) * 100 * (h * c/ kb))/(np.log(gamma) - np.log(Q1Int[i,:]/Q3Int[i,:]))
    Temps.append(Rotat_Temp)
        
TIMEv2 = ['12:42', '12:11', '11:42', '11:12', '10:43', '10:13', '09:44', '09:14', '08:44', '08:14', '07:45']
# Temps = Temps[:,0:20]
#Temps[Temps > 800] = np.nan

plt.figure()
plt.imshow(Temps, cmap='coolwarm', vmax = 600, vmin = 300)
cbar = plt.colorbar()
cbar.ax.tick_params(labelsize=15)
cbar.set_label(r'Temperature (K)' + u'\u00B1' + ' 167K', fontsize=17.5)
plt.title(r'Temperatures of $H_{3}^{+}$ from Uranus observations from 07:30 to 12:57 UTC on $11^{th}$ October 2016', pad = 45, fontsize = 22.5)
plt.xlabel('Arbitrary Longitude across Uranus ($^\circ$)', fontsize=17.5, labelpad=11)
plt.ylabel('Latitude across Uranus ($^\circ$)', fontsize=17.5, labelpad=11)
plt.xlabel('Arbitrary Longitude across Uranus (Pixels)', fontsize=17.5, labelpad=10)
plt.ylabel('UTC Time (HH:MM +- 00:15)', fontsize=17.5, labelpad=10)
#plt.yticks(np.arange(11), TIMEv2[0:11])

#Temperature Errors
Temperature_Err = (np.nanmean(Mapping_Q1A[:,1:22])/np.nanmean(Mapping_Q3B[:,0:22]))*np.sqrt((0.06/np.nanmean(Mapping_Q1A[:,1:22]))**2+(0.05/np.nanmean(Mapping_Q3B[:,0:22]))**2)
Temps_Err = Temperature_Err/(np.nanmean(Mapping_Q1A[:,1:22])/np.nanmean(Mapping_Q3B[:,0:22]))
Temp_Err = np.nanmean(Temps)*(Temps_Err/math.log(np.nanmean(Mapping_Q1A[:,1:22])/np.nanmean(Mapping_Q3B[:,0:22])))

print(Temp_Err) #Look into why this is so big

#%% Now we try to pair Q3 and Q1 together 

for n in range(11):
    if n == 0:
        Total_IRTF_DataQ1 = Q1_Data_Capture[0]
    else:
        Total_IRTF_DataQ1 = np.add(Total_IRTF_DataQ1, Q1_Data_Capture[n])

Scaling = np.nanmean(np.gradient(Wavelength_O5[1360:1363]))/np.nanmean(np.gradient(Wavelength_O4[1214:1217]))
Scaling2 = 41.2337/41.6508
for n in range(11):
    if n == 0:
        Total_IRTF_DataQ3 = Q3_Data_Capture[0]
        Total_IRTF_DataQ3 = f.zoom(Q3_Data_Capture[0], (Scaling2, Scaling), mode='wrap')
    else:
        Testing = f.zoom(Q3_Data_Capture[n], (Scaling2, Scaling), mode='wrap')
        Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, Testing)
        #Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, Q3_Data_Capture[n])
        
#Total_IRTF_DataQ3 = f.zoom(Total_IRTF_DataQ3, (0.9899857865875326, 1), mode='wrap')
for n in range(11):
    if n == 0:
        Total_IRTF_DataQ31 = Q31_Data_Capture[0]
    else:
        Total_IRTF_DataQ31 = np.add(Total_IRTF_DataQ31, Q31_Data_Capture[n])

    #%%
Total_IRTF_Data_Alt = acre(Total_IRTF_DataQ3, width = 2, verbose = False)
Total_IRTF_Data_Alt2 = acre(Total_IRTF_DataQ1, width = 2, verbose = False)

ex_A1_Q1T = []
ex_Err_A1_Q1T = []
Pos_Q1T = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ1
    datai = dataI[6+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    a1_pointP = a1_pointP[0][0] - 1110
    ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    ex_A1_Q1T.append(ex_pQ1P.a1)
    ex_Err_A1_Q1T.append(ex_eQ1P[1])
    Pos_Q1T.append(ex_pQ1P.a1 + 1110)

#Check if Q2 and Q31 are strong enough to include
ex_A1_Q2 = []
ex_Err_A1_Q2 = []
Pos_Q2 = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ3
    datai = dataI[6+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[470:495]))
    a1_pointP = a1_pointP[0][0] - 450
    ex_resultQ2P = gmodel.fit(datai[450:500], x=np.arange(50), a0=np.nanmax(datai[470:495]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ2P = SimpleNamespace(**ex_resultQ2P.best_values)
    ex_eQ2P = np.sqrt(np.diag(ex_resultQ2P.covar))
    Pos_Q2.append(ex_pQ2P.a1 + 450)
    
ex_A1_Q3T = []
ex_Err_A1_Q3T = []
Pos_Q3T = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ3
    datai = dataI[6+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1380:1400]))
    a1_pointP = a1_pointP[0][0] - 1340
    ex_resultQ3P = gmodel.fit(datai[1340:1440], x=np.arange(100), a0=np.nanmax(datai[1380:1440]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ3P = SimpleNamespace(**ex_resultQ3P.best_values)
    ex_eQ3P = np.sqrt(np.diag(ex_resultQ3P.covar))
    Pos_Q3T.append(ex_pQ3P.a1 + 1340)
    # plt.figure()
    # plt.plot(np.arange(100), datai[1290:1390])

# ex_A1_Q31 = []
# ex_Err_A1_Q31 = []
# Pos_Q31 = []

# o = 0
# for m in range(25):
#     dataI = Total_IRTF_DataQ31
#     datai = dataI[14+m, :]
#     gmodel = Model(gauss_fit)
#     a1_pointP = np.where(datai == np.nanmax(datai[1480:1500]))
#     a1_pointP = a1_pointP[0][0] - 1460
#     ex_resultQ31P = gmodel.fit(datai[1460:1520], x=np.arange(60), a0=np.nanmax(datai[1480:1500]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
#     ex_pQ31P = SimpleNamespace(**ex_resultQ31P.best_values)
#     ex_eQ31P = np.sqrt(np.diag(ex_resultQ31P.covar))
#     Pos_Q31.append(ex_pQ31P.a1 + 1460)
    
shift_list = []
shift_list2 = []
shift_list3 = []

for n in range(24):
    shift_list.append(Pos_Q3T[n] - Pos_Q1T[n])
    shift_list2.append(Pos_Q2[n] - Pos_Q1T[n])
    # shift_list3.append(Pos_Q31[n] - Pos_Q1T[n])
    
shift_amount1 = np.nanmean(shift_list[0:23])
shift_amount2 = np.nanmean(shift_list2[0:23])
# shift_amount3 = np.nanmean(shift_list3[0:23])

#%%
ii = 0
CombinedImage = []

for nn in range(11):
    Q1Wave = Q1_Data_Capture[nn]
    #Q3Wave = Total_IRTF_DataQ3 = zoom(Q3_Data_Capture[nn], (1.0, 0.9697589886378948), mode='constant')
    Q3Wave = f.zoom(Q3_Data_Capture[nn], (Scaling2, Scaling), mode='wrap')
    Q31Wave = Q31_Data_Capture[nn]
    for n in range(35):
        Shifted_Q3 = Q3Wave[n+3,:]
        Shifted_Q31 = Q31Wave[n+3,:] 
        shift_amount = shift_amount1 #+ 0.1533678848228073*2.5
        shift_amount2 = shift_amount2 #+ 0.1533678848228073*2.5
        # shift_amount3 = shift_amount3
        shifted_Q3 = f.shift(Shifted_Q3, shift=[-shift_amount], mode='wrap') #41.6508 for Order 5
        shifted_Q2 = f.shift(Shifted_Q3, shift=[-shift_amount2], mode='wrap')
        # shifted_Q31 = f.shift(Shifted_Q3, shift=[-shift_amount3], mode='wrap')
        if n == 0:
            Order_5_shift = np.vstack((Q3Wave[0:3, :], shifted_Q3))
            Order_5_shift1 = np.vstack((Q3Wave[0:3, :], shifted_Q2))
            # Order_5_shift2 = np.vstack((Q3Wave[0:10, :], shifted_Q31))
        elif n == 34:
            Order_5_shift = np.vstack((Order_5_shift, Q3Wave[37:41, :]))
            Order_5_shift1 = np.vstack((Order_5_shift1, Q3Wave[37:41, :]))
            # Order_5_shift2 = np.vstack((Order_5_shift2, Q3Wave[44:51, :]))
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Q3))
            Order_5_shift1 = np.vstack((Order_5_shift1, shifted_Q2))
            # Order_5_shift2 = np.vstack((Order_5_shift2, shifted_Q31))
        ii += 1
    if nn == 9:
        ComboFigure  = (Q1Wave[0:40,:] + Order_5_shift[:,0:2048])/2
        CombinedImage.append(ComboFigure)        
    else:
        ComboFigure  = (Q1Wave[0:40,:] + Order_5_shift[:,0:2048])/2
        CombinedImage.append(ComboFigure)

for n in range(11):
    if n == 0:
        Total_IRTF_Data = CombinedImage[0]
    else:
        Total_IRTF_Data = np.add(Total_IRTF_Data, CombinedImage[n])

#%%
ex_A1_Q = []
ex_Err_A1_Q = []
Q_IntErr = []

for m in range(23):
    #print(m)
    dataI = Total_IRTF_Data
    datai = dataI[6+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    a1_pointP = a1_pointP[0][0] - 1110
    ex_resultQP = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQP = SimpleNamespace(**ex_resultQP.best_values)
    ex_eQP = np.sqrt(np.diag(ex_resultQP.covar))
    ex_A1_Q.append(ex_pQP.a1+1110)
    ex_Err_A1_Q.append(ex_eQP[1])

Olap = []
Err_lap = []
for m in range(23):
    Overlap = ex_A1_Q[m] - Pos_Q1T[m]
    Olap.append(Overlap)
    Err = np.sqrt((ex_Err_A1_Q1T[m])**2+(ex_Err_A1_Q[m])**2)
    Err_lap.append(Err)
    print(Overlap, Err)
    
print(np.nanmean(Olap), np.nanmean(Err_lap))

# for o in range(242):
#     QIntErr = ex_INTSQ[o]*np.sqrt((ex_Err_A0_Q[o]/ex_A0_Q[o])**2 + (ex_Err_A2_Q[o]/ex_A2_Q[o])**2)
#     Q_IntErr.append(QIntErr)

#%% Now we add boxes
ex_A0_Q = []
ex_A2_Q = []
ex_Err_A0_Q= []
ex_Err_A2_Q = []
ex_A1_Q = []
ex_Err_A1_Q = []
ex_FWHMQ = []
ex_INTSQ = []

wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225])) #Confirm this!!!!!!!! Then go through each fitting to make sure the correct thing is happening!!
LIST = 7, 7, 7, 7, 7, 7, 7, 7, 7, 7, 7  #A rough idea of the starting point on each Image

for n in range(11):
    dataI = CombinedImage[n]
    dataI[dataI < -0.004] = 0
    # if n == 9:
    #     plt.figure()
    #     plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    #dataI[dataI > 0.01] == 0.001
    for m in range(24):
        datai = dataI[7+m, :]
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQP = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQP = SimpleNamespace(**ex_resultQP.best_values)
        ex_eQP = np.sqrt(np.diag(ex_resultQP.covar))
        A0 = ex_pQP.a0
        ex_A0_Q.append(A0)
        A2 = ex_pQP.a2
        ex_A2_Q.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQP[0]
        ex_Err_A0_Q.append(err_A0)
        err_A2 = ex_eQP[2]
        ex_Err_A2_Q.append(err_A2*wav_pixel_ratio)
        ex_A1_Q.append(ex_pQP.a1+1110)
        ex_Err_A1_Q.append(ex_eQP[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQ.append(ex_FWHM)
        ex_INTSQ.append(A0*ex_FWHM*(10**6))
        # if n == 9:
        #     plt.figure()
        #     plt.plot(XX, datai[1110:1310], 'bo')
        #     plt.plot(XX, ex_resultQP.best_fit)

Q_IntErr = []
for o in range(264):
    QIntErr = ex_INTSQ[o]*np.sqrt((ex_Err_A0_Q[o]/ex_A0_Q[o])**2 + (ex_Err_A2_Q[o]/ex_A2_Q[o])**2)
    Q_IntErr.append(QIntErr)

#%% 
A1_Q = []
#Q1 is approx from the reduction side
Q1_Pos = 1212.458357810168

o = 0
for n in range(264):
    A1 = ((ex_A1_Q[n]-Q1_Pos)*-wav_pixel_ratio)
    A1_Q.append(A1)
    if (n + 1) % 24 == 0:
        o += 24

#%%
c = 299792458 
Velocity_Q = []
lamda_H3_O4 = 3.95295 
lamda_H3_O5 = 3.98558
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846

for n in range(264):
    V = ((A1_Q[n]/lamda_H3_O4)*c)
    Velocity_Q.append(V/1000)
    
Velocity_Q_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(264):
    Err = Velocity_Q[n]*np.sqrt((ex_Err_A1_Q[n]/A1_Q[n])**2)
    Velocity_Q_Err.append(np.sqrt(Err**2))
    
Approx_Longitudes = [-88.72556327833178, -79.95553680388429, -72.36098621142304, -65.3970422375618, -58.60042435368884, -51.74045456148676, -44.6727749091867, -37.30101350396452, -29.568948625029716, -21.465373053143523, -13.03294197002669, -4.3717430760420015]
Approx_Longitudes2 = np.flip(Approx_Longitudes)*-1
Approx_Longs = np.concatenate((Approx_Longitudes, Approx_Longitudes2))

R_Uranus = 28059*1000 #m
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

#So the speed of Uranus at surface is 
Circumference_Uranus = 2*np.pi*R_Uranus
Period_Uranus = Circumference_Uranus/Time_s_Uranus #Look at scaling

#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846

Limb_velocity_pix_O4 = (Period_Uranus/1000)
Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 - (0.1*Limb_velocity_pix_O4)
Limb_velocity_pix_O4_80 = Limb_velocity_pix_O4 - (0.2*Limb_velocity_pix_O4)

Limb_velocity_pix_O4_test = []
for a in range(26):
    aa = 12.5 - a
    Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin(aa*np.pi/24))

Planet_rotation_O4 = Limb_velocity_pix_O4_test[1:25]
# Planet_rotation_O5 = np.linspace(Limb_velocity_pix_O5, -1*Limb_velocity_pix_O5, 22)

Velocity_Err_QP = []
Velocity_Err_QN = []

for n in range(264):
    Velocity_Err_QP.append(Velocity_Q[n] + Velocity_Q_Err[n])
    Velocity_Err_QN.append(Velocity_Q[n] - Velocity_Q_Err[n])
    
Ints_Err_QP = []
Ints_Err_QN = []
    
for n in range(264):
    Ints_Err_QP.append(ex_INTSQ[n] + Q_IntErr[n])
    Ints_Err_QN.append(ex_INTSQ[n] - Q_IntErr[n])    

TIME = ['07:45', '08:14', '08:44', '09:14', '09:44', '10:13', '10:43', '11:12', '11:42', '12:11', '12:42']

o = 0
for n in range(11):
    fig, ax = plt.subplots(figsize=(10,8))
    ax2 = ax.twinx()
    ax.set_title('~30 minute exposure at ' + str(TIME[n]) +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
    ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
    ax.plot(Approx_Longs, np.flip(Velocity_Q[o:o+24]), color='b', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
    #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
    ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QN[o:o+24]), np.flip(Velocity_Err_QP[o:o+24]), color='b', alpha=0.5)
    #ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
    ax2.plot(Approx_Longs, np.flip(ex_INTSQ[o:o+24]), color='r', label='IR Intensity', lw = 5)
    ax2.fill_between(Approx_Longs, np.flip(Ints_Err_QP[o:o+24]), np.flip(Ints_Err_QN[o:o+24]), color='r', alpha=0.5)
    #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
    #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
    ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
    ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
    ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
    ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
    ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    ax.set_xlim(-90, 90) #0, 3601
    ax.set_ylim(-3.5, 3.5)
    ax2.set_ylim(0, 1.0)
    ax.legend(loc='upper left', fontsize=25)
    ax2.legend(loc='upper right', fontsize=25)
    o += 24

#%% Lets make a mask for the Total figure of Q1 + Q3 to use on following stuff
Mean = np.nanmean(Total_IRTF_Data[10:31,500:1000])
Std = np.nanstd(Total_IRTF_Data[10:31,500:1000])
FilterMask = (Total_IRTF_Data > Mean + Std)
AntiMask = np.invert(FilterMask)

# plt.figure()
# plt.imshow(FilterMask)

# plt.figure()
# plt.imshow(AntiMask)

#Now we make an example of the noise to put back into the data
# MeanNoiseWave = np.nanmean(Total_IRTF_Data[10:51, 650:850], axis=0)
#MeanNoiseWave = MeanNoiseWave*AntiMask[20, 1100:1300]
# plt.figure()
# plt.plot(np.arange(100), MeanNoiseWave)
#NoiseMaker = AntiMask[10:51, 700:800]*MeanNoiseWave

#%% Now we extend out the analysis over more pixels to see if we can find more pixels even outside the disk AS = Aurora Set
A0_QAS = []
ex_A0_QAS = []
Err_A0_QAS = []
A2_QAS = []
ex_A2_QAS = []
Err_A2_QAS = []
ex_FWHMQAS = []
FWHMQAS = []
ex_INTSQAS = []
INTSQAS = []
ex_Err_A0_QAS = []
ex_Err_A2_QAS = []
ex_A1_QAS = []
ex_Err_A1_QAS = []
QAS_IntErr = []

#Lets clean the image a bit more
CleaningImage = CombinedImage[8]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
CleaningImage[31:34,1209:1211] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[9]

for n in range(0,2,1):
    if n == 0:
        dataI = CleaningImage
        plt.figure()
        plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    else:
        dataI = CleaningImage2
        plt.figure()
        plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    #dataI[dataI < -0.004] = 0
    #dataI[dataI > 0.01] = 0.001
    dataI[dataI < -0.004] = 0
    for m in range(25): #31
        datai = dataI[7+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
        a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QAS.append(ex_eQPAS[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQAS.append(ex_FWHM)
        ex_INTSQAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai, 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')
        #         print(ex_pQPAS.a1, 23-m)

for o in range(50):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(50):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[2:o+22]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)
    if n == 24:
        o = 25

Velocity_QAS = []
c = 299792458 
for n in range(50):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(50):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2 + (np.nanstd(ex_A1_QAS[o+2:o+21])/np.nanmean(ex_A1_QAS[o+2:o+21]))**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    if n == 24:
        o = 25
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(50):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(50):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

# Limb_velocity_pix_O4_test = []
# for a in range(26):
#     aa = 12.5 - a
#     Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin(aa*np.pi/24))

# Planet_rotation_O4 = Limb_velocity_pix_O4_test[1:25]

R_Uranus = 28559*1000 #m (2*28559*np.pi*np.cos(np.deg2rad(Lats)))
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

Lats = np.load('LatitudesS.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
for o in range(26):
    Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2))/ (17.24*60*60)
    Limb_velocity_pix_O4 = (Period_Uranus/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for o in range(26):
    Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2)))/(2*np.pi))
for a in range(13):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
for a in range(13):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))

#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
Planet_rotation_O4b = Planet_rotation_O4*1.1
    
# # NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# # EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
# SOUTH_SLIT = np.load('SOUTH_SLIT_CONFIG.npy', allow_pickle=True)
# Longs = np.load('LongitudesS.npy')
# South_Pos = []
# for o in range(25):
#     South_Pos.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)

# Approx_Longs = South_Pos

# fig, ax = plt.subplots(figsize=(10,8))
# ax2 = ax.twinx()
# ax.set_title('~30 minute exposure at ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
# ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# ax.plot(Approx_Longs, np.flip(Velocity_QAS[25:50]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[25:50]), np.flip(Velocity_Err_QAS_P[25:50]), color='b', alpha=0.5)
# ax2.plot(Approx_Longs, np.flip(ex_INTSQAS[25:50]), color='r', label='IR Intensity', lw = 5)
# ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N[25:50]), np.flip(INTS_Err_QAS_P[25:50]), color='r', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
# ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(-4, 4)
# ax2.set_ylim(0, 0.75)
# ax.legend(loc='lower left', fontsize=25)
# ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

# fig, ax = plt.subplots(figsize=(10,8))
# ax2 = ax.twinx()
# ax.set_title('~30 minute exposure at ' + '11:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
# ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
# #ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:25]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:25]), np.flip(Velocity_Err_QAS_P[0:25]), color='b', alpha=0.5)
# ax2.plot(Approx_Longs, np.flip(ex_INTSQAS[0:25]), color='r', label='IR Intensity', lw = 5)
# ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N[0:25]), np.flip(INTS_Err_QAS_P[0:25]), color='r', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
# ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(-4, 4)
# ax2.set_ylim(0, 0.75)
# ax.legend(loc='lower left', fontsize=25)
# ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASA = []
Vel_UranusNASA_Err = []
Int_UranusNASA = []
Int_UranusNASA_Err = []

Vel_UranusNASA.append(np.flip(Velocity_QAS[0:25]))
Vel_UranusNASA.append(np.flip(Velocity_QAS[25:50]))
Vel_UranusNASA_Err.append(np.flip(Velocity_QAS_Err[0:25]))
Vel_UranusNASA_Err.append(np.flip(Velocity_QAS_Err[25:50]))
Int_UranusNASA.append(np.flip(ex_INTSQAS[0:25]))
Int_UranusNASA.append(np.flip(ex_INTSQAS[25:50]))
Int_UranusNASA_Err.append(np.flip(QAS_IntErr[0:25]))
Int_UranusNASA_Err.append(np.flip(QAS_IntErr[25:50]))

#%% Now we check each set individually #Set 1 Check figure 1
A0_QNAS = []
ex_A0_QNAS = []
Err_A0_QNAS = []
A2_QNAS = []
ex_A2_QNAS = []
Err_A2_QNAS = []
ex_FWHMQNAS = []
FWHMQNAS = []
ex_INTSQNAS = []
INTSQNAS = []
ex_Err_A0_QNAS = []
ex_Err_A2_QNAS = []
ex_A1_QNAS = []
ex_Err_A1_QNAS = []
QNAS_IntErr = []

# CleaningImage[19:21,1213:1214] = np.nanmean(CombinedImage[0])
CleaningImage = CombinedImage[0]
#CleaningImage[28,1205:1207] = np.nanmean(CombinedImage[0])

CleaningImage2 = CombinedImage[1]

for n in range(0,2,1):
    for m in range(26): #27
        if n == 0:
            dataI = CleaningImage
            datai = dataI[6+m, 1100:1300] #6
            # datai = datai+(acre(CombinedImage[8][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            # if m == 0:
            #     plt.figure()
            #     plt.imshow(dataI, cmap='gist_gray')
        else:
            dataI = CleaningImage2
            datai = dataI[6+m, 1100:1300]
            #datai = datai+(acre(CombinedImage[9][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            # if m == 0:
            #     plt.figure()
            #     plt.imshow(dataI, cmap='gist_gray')
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
        gmodel = Model(gauss_fit)
        if m == 0:
            a1_pointP = np.where(datai == np.nanmax(datai[111:116]))
            a1_pointP = a1_pointP[0][0]            
        else:
            a1_pointP = np.where(datai == np.nanmax(datai[110:117]))
            a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[105:122]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        #datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1100)
        ex_Err_A1_QNAS.append(ex_eQPAS[1])
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 1:
        #         plt.figure()
        #         plt.plot(XX, datai, 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(52):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

#%%
A1_QNAS = []

o = 0
for n in range(52):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o:o+25]))*-wav_pixel_ratio)
    A1_QNAS.append(A1)
    if n == 25:
        o = 26

Velocity_QNAS = []
c = 299792458 
for n in range(52):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(52):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/(ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS)))**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(52):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])
    
INTS_Err_QNAS_N = []
INTS_Err_QNAS_P = []

for n in range(52):
    INTS_Err_QNAS_N.append(ex_INTSQNAS[n]-QNAS_IntErr[n])
    INTS_Err_QNAS_P.append(ex_INTSQNAS[n]+QNAS_IntErr[n])    

Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
Lats = np.load('LatitudesE.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
for o in range(27):
    Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad(Lats[o]))/ (17.24*60*60)
    Limb_velocity_pix_O4 = (Period_Uranus/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for o in range(27):
    Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad(Lats[o])))/(2*np.pi))
for a in range(14):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((13-a)*np.pi/(2*Planet_diameter[a])))
for a in range(14):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((13-a)*np.pi/(2*Planet_diameter[a])))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
    
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesE.npy')
South_Pos = []
for o in range(27):
    if o == 0:
        South_Pos.append((Longs[o] + Longs[27] - 90 - 90)/4)
    elif o == 26:
        South_Pos.append((Longs[26] + Longs[53] + 90 + 90)/4)
    else:
        South_Pos.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)

Approx_Longs = South_Pos

E_Rotation = Planet_rotation_O4
E_Longs = South_Pos

#%%
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '08:14' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:27], Planet_rotation_O4[1:27], color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:27], np.flip(Velocity_QNAS[26:52]), color='b', ls= '--', label='Combined Q(1,0), Q(3,0) and Q(3,1) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:27], np.flip(Velocity_Err_QNAS_N[26:52]), np.flip(Velocity_Err_QNAS_P[26:52]), color='b', alpha=0.5)
ax2.plot(Approx_Longs[1:27], np.flip(ex_INTSQNAS[26:52]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:27], np.flip(INTS_Err_QNAS_N[26:52]), np.flip(INTS_Err_QNAS_P[26:52]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QNAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 21), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '07:45' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:27], Planet_rotation_O4[1:27], color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:27], np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:27], np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(Approx_Longs[1:27], np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:27], np.flip(INTS_Err_QNAS_N[0:26]), np.flip(INTS_Err_QNAS_P[0:26]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QNAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 21), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASE = []
Vel_UranusNASE_Err = []
Int_UranusNASE = []
Int_UranusNASE_Err = []

Vel_UranusNASE.append(np.flip(Velocity_QNAS[0:26])) #Keep working on the velocity model as it looks like it works but south is similar to equator?
Vel_UranusNASE.append(np.flip(Velocity_QNAS[26:52]))
Vel_UranusNASE_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNASE_Err.append(np.flip(Velocity_QNAS_Err[26:52]))
Int_UranusNASE.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNASE.append(np.flip(ex_INTSQNAS[26:52]))
Int_UranusNASE_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNASE_Err.append(np.flip(QNAS_IntErr[26:52]))

#%% Set 2 #2 #3 #Look into averaging over Equator, North, Southern
A0_QNAS = []
ex_A0_QNAS = []
Err_A0_QNAS = []
A2_QNAS = []
ex_A2_QNAS = []
Err_A2_QNAS = []
ex_FWHMQNAS = []
FWHMQNAS = []
ex_INTSQNAS = []
INTSQNAS = []
ex_Err_A0_QNAS = []
ex_Err_A2_QNAS = []
ex_A1_QNAS = []
ex_Err_A1_QNAS = []
QNAS_IntErr = []

for n in range(0,2,1):
    dataI = CombinedImage[n+2]
    # dataI[dataI > 0.01] = 0.001
    for m in range(25): #31
        if n == 0:
            #datai = acre(dataI[15+m, 1100:1300], width = 2, verbose = False)
            datai = dataI[7+m, 1100:1300] # 8
            # datai = datai+(acre(CombinedImage[8][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            if m == 0:
                plt.figure()
                plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
        else:
            #datai = acre(dataI[15+m, 1100:1300], width = 2, verbose = False)
            datai = dataI[7+m, 1100:1300]
            dataI[25:27,106:111] = np.nanmean(dataI)
            #datai = datai+(acre(CombinedImage[9][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            if m == 0:
                plt.figure()
                plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[111:118]))
        a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[111:118]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1100)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if m == 0:
        #         plt.figure()
        #         plt.plot(XX, datai, 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(50):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

#%%
A1_QNAS = []

o = 0
for n in range(50):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS))*-wav_pixel_ratio)
    A1_QNAS.append(A1)
    if n == 24:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(50):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(50):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(50):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])
    
INTS_Err_QNAS_N = []
INTS_Err_QNAS_P = []

for n in range(50):
    INTS_Err_QNAS_N.append(ex_INTSQNAS[n]-QNAS_IntErr[n])
    INTS_Err_QNAS_P.append(ex_INTSQNAS[n]+QNAS_IntErr[n])  

R_Uranus = 28559*1000 #m
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
Lats = np.load('LatitudesS.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
for o in range(27):
    Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad(Lats[o]))/ (17.24*60*60)
    Limb_velocity_pix_O4 = (Period_Uranus/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for o in range(25):
    Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad(Lats[o])))/(2*np.pi))
for a in range(13):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
for a in range(13):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
    
# Longs = np.load('LongitudesS.npy')
# South_Pos = []
# for o in range(25):
#     South_Pos.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)

# Approx_Longs = South_Pos

# S_Longs = South_Pos
# S_Rotation = Planet_rotation_O4

# fig, ax = plt.subplots(figsize=(10,8))
# ax2 = ax.twinx()
# ax.set_title('~30 minute exposure at ' + '09:14' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# ax.plot(Approx_Longs, np.flip(Velocity_QNAS[25:50]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QNAS_N[25:50]), np.flip(Velocity_Err_QNAS_P[25:50]), color='b', alpha=0.5)
# ax2.plot(Approx_Longs, np.flip(ex_INTSQNAS[25:50]), color='r', label='IR Intensity', lw = 5)
# ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QNAS_N[25:50]), np.flip(INTS_Err_QNAS_P[25:50]), color='r', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# #ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
# #ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
# ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QNAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(-4, 4)
# ax2.set_ylim(0, 0.75)
# ax.legend(loc='lower left', fontsize=25)
# ax2.legend(loc='upper right', fontsize=25)

# fig, ax = plt.subplots(figsize=(10,8))
# ax2 = ax.twinx()
# ax.set_title('~30 minute exposure at ' + '08:44' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
# #ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# ax.plot(Approx_Longs, np.flip(Velocity_QNAS[0:25]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QNAS_N[0:25]), np.flip(Velocity_Err_QNAS_P[0:25]), color='b', alpha=0.5)
# ax2.plot(Approx_Longs, np.flip(ex_INTSQNAS[0:25]), color='r', label='IR Intensity', lw = 5)
# ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QNAS_N[0:25]), np.flip(INTS_Err_QNAS_P[0:25]), color='r', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# #ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
# #ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
# ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QNAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(-4, 4)
# ax2.set_ylim(0, 0.75)
# ax.legend(loc='lower left', fontsize=25)
# ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASS = []
Vel_UranusNASS_Err = []
Int_UranusNASS = []
Int_UranusNASS_Err = []

Vel_UranusNASS.append(np.flip(Velocity_QNAS[0:25]))
Vel_UranusNASS.append(np.flip(Velocity_QNAS[25:50]))
Vel_UranusNASS_Err.append(np.flip(Velocity_QNAS_Err[0:25]))
Vel_UranusNASS_Err.append(np.flip(Velocity_QNAS_Err[25:50]))
Int_UranusNASS.append(np.flip(ex_INTSQNAS[0:25]))
Int_UranusNASS.append(np.flip(ex_INTSQNAS[25:50]))
Int_UranusNASS_Err.append(np.flip(QNAS_IntErr[0:25]))
Int_UranusNASS_Err.append(np.flip(QNAS_IntErr[25:50]))

#%% Set 3 #4 #5 Take average away from the Ints to see increases and decreases!!!!
A0_QNAS = []
ex_A0_QNAS = []
Err_A0_QNAS = []
A2_QNAS = []
ex_A2_QNAS = []
Err_A2_QNAS = []
ex_FWHMQNAS = []
FWHMQNAS = []
ex_INTSQNAS = []
INTSQNAS = []
ex_Err_A0_QNAS = []
ex_Err_A2_QNAS = []
ex_A1_QNAS = []
ex_Err_A1_QNAS = []
QNAS_IntErr = []

for n in range(0,2,1):
    dataI = CombinedImage[n+4]
    for m in range(25): #25
        if n == 0:
            datai = dataI[6+m, 1100:1300] #6
            # datai = datai+(acre(CombinedImage[8][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            if m == 0:
                plt.figure()
                plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01) 
        else:
            datai = dataI[6+m, 1100:1300]
            #datai = datai+(acre(CombinedImage[9][13+m, 1100:1300], width = 2, verbose = False)*AntiMask[13+m, 1100:1300])
            if m == 0:
                plt.figure()
                plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)   
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[100:122]))
        a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[100:122]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 1:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(50):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

#%%
A1_QNAS = []

o = 0
for n in range(50):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS))*-wav_pixel_ratio)
    A1_QNAS.append(A1)
    if n == 24:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(50):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(50):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(50):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])
    
INTS_Err_QNAS_N = []
INTS_Err_QNAS_P = []

for n in range(50):
    INTS_Err_QNAS_N.append(ex_INTSQNAS[n]-QNAS_IntErr[n])
    INTS_Err_QNAS_P.append(ex_INTSQNAS[n]+QNAS_IntErr[n])  

Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
Lats = np.load('LatitudesN.npy')
Longs = np.load('LongitudesN.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
for o in range(25):
    Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+1]+Lats[o+26]+Lats[o+27])/4))/ (17.24*60*60)
    Limb_velocity_pix_O4 = (Period_Uranus/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for o in range(25):
    Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+1]+Lats[o+26]+Lats[o+27])/4)))/(2*np.pi))
for a in range(13):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
for a in range(13):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
    
Approx_Longs = []

for o in range(25):
    Approx_Longs.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)
    
N_Longs = Approx_Longs
N_Rotation = Planet_rotation_O4

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '10:13' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QNAS[25:50]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QNAS_N[25:50]), np.flip(Velocity_Err_QNAS_P[25:50]), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQNAS[25:50]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QNAS_N[25:50]), np.flip(INTS_Err_QNAS_P[25:50]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 11), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '09:44' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QNAS[0:25]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QNAS_N[0:25]), np.flip(Velocity_Err_QNAS_P[0:25]), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQNAS[0:25]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QNAS_N[0:25]), np.flip(INTS_Err_QNAS_P[0:25]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 11), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASN = []
Vel_UranusNASN_Err = []
Int_UranusNASN = []
Int_UranusNASN_Err = []

Vel_UranusNASN.append(np.flip(Velocity_QNAS[0:25]))
Vel_UranusNASN.append(np.flip(Velocity_QNAS[25:50]))
Vel_UranusNASN_Err.append(np.flip(Velocity_QNAS_Err[0:25]))
Vel_UranusNASN_Err.append(np.flip(Velocity_QNAS_Err[25:50]))
Int_UranusNASN.append(np.flip(ex_INTSQNAS[0:25]))
Int_UranusNASN.append(np.flip(ex_INTSQNAS[25:50]))
Int_UranusNASN_Err.append(np.flip(QNAS_IntErr[0:25]))
Int_UranusNASN_Err.append(np.flip(QNAS_IntErr[25:50])) #This needs work with the velocities to check if they match up!!!

#%% Set 4 #6 #7
A0_QNAS = []
ex_A0_QNAS = []
Err_A0_QNAS = []
A2_QNAS = []
ex_A2_QNAS = []
Err_A2_QNAS = []
ex_FWHMQNAS = []
FWHMQNAS = []
ex_INTSQNAS = []
INTSQNAS = []
ex_Err_A0_QNAS = []
ex_Err_A2_QNAS = []
ex_A1_QNAS = []
ex_Err_A1_QNAS = []
QNAS_IntErr = []

for n in range(0,2,1):
    dataI = CombinedImage[n+6]
    # dataI[dataI < -0.004] = 0
    # dataI[dataI > 0.01] = 0.001
    for m in range(27): #27
        datai = dataI[5+m, 1100:1300] #6
        #datai = datai+(MeanNoiseWave*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[100:122]))
        a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[100:122]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 1:
        #         plt.figure()
        #         plt.plot(XX, datai, 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(54):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

#%%
A1_QNAS = []

o = 0
for n in range(54):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS))*-wav_pixel_ratio)
    A1_QNAS.append(A1)
    if n == 25:
        o = 26

Velocity_QNAS = []
c = 299792458 
for n in range(54):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(54):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(54):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])
    
INTS_Err_QNAS_N = []
INTS_Err_QNAS_P = []

for n in range(54):
    INTS_Err_QNAS_N.append(ex_INTSQNAS[n]-QNAS_IntErr[n])
    INTS_Err_QNAS_P.append(ex_INTSQNAS[n]+QNAS_IntErr[n])  

# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy', allow_pickle=True)
#NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy', allow_pickle=True)
Equator_Pos = []
for n in range(27):
    Equator_Pos.append(EQUATOR_SLIT[3][n] - 180)
Approx_Longs = Equator_Pos

Lats = np.load('LatitudesE.npy')
Longs = np.load('LongitudesE.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
for o in range(27):
    Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+1]+Lats[o+26]+Lats[o+27])/4))/ (17.24*60*60)
    Limb_velocity_pix_O4 = (Period_Uranus/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for o in range(27):
    Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+1]+Lats[o+26]+Lats[o+27])/4)))/(2*np.pi))
for a in range(14):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((13-a)*np.pi/(2*Planet_diameter[a])))
for a in range(14):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((13-a)*np.pi/(2*Planet_diameter[a])))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:12' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:26], Planet_rotation_O4[1:26], color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:26], np.flip(Velocity_QNAS[28:53]), color='b', ls= '--', label='Combined Q(1,0), Q(3,0) and Q(3,1) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:26], np.flip(Velocity_Err_QNAS_N[28:53]), np.flip(Velocity_Err_QNAS_P[28:53]), color='b', alpha=0.5)
ax2.plot(Approx_Longs[1:26], np.flip(ex_INTSQNAS[28:53]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:26], np.flip(INTS_Err_QNAS_N[28:53]), np.flip(INTS_Err_QNAS_P[28:53]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 11), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '10:43' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:26], Planet_rotation_O4[1:26], color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:26], np.flip(Velocity_QNAS[1:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:26], np.flip(Velocity_Err_QNAS_N[1:26]), np.flip(Velocity_Err_QNAS_P[1:26]), color='b', alpha=0.5)
ax2.plot(Approx_Longs[1:26], np.flip(ex_INTSQNAS[1:26]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:26], np.flip(INTS_Err_QNAS_N[1:26]), np.flip(INTS_Err_QNAS_P[1:26]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-5, 5, 11), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-5, 5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASE.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNASE.append(np.flip(Velocity_QNAS[27:53]))
Vel_UranusNASE_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNASE_Err.append(np.flip(Velocity_QNAS_Err[27:53]))
Int_UranusNASE.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNASE.append(np.flip(ex_INTSQNAS[27:53]))
Int_UranusNASE_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNASE_Err.append(np.flip(QNAS_IntErr[27:53]))

#%% Set 5 #10
A0_QNAS = []
ex_A0_QNAS = []
Err_A0_QNAS = []
A2_QNAS = []
ex_A2_QNAS = []
Err_A2_QNAS = []
ex_FWHMQNAS = []
FWHMQNAS = []
ex_INTSQNAS = []
INTSQNAS = []
ex_Err_A0_QNAS = []
ex_Err_A2_QNAS = []
ex_A1_QNAS = []
ex_Err_A1_QNAS = []
QNAS_IntErr = []

for n in range(0,1,1):
    dataI = CombinedImage[n+10]
    # dataI[dataI < -0.004] = 0
    # dataI[dataI > 0.01] = 0.001
    for m in range(24): #24
        datai = dataI[6+m, 1100:1300] #6
        #datai = datai+(MeanNoiseWave*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[100:120]))
        a1_pointP = a1_pointP[0][0]
        ex_resultQPAS = gmodel.fit(datai, x=XX, a0=np.nanmax(datai[100:120]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*wav_pixel_ratio)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*wav_pixel_ratio)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*wav_pixel_ratio)
        ex_FWHM = A2*wav_pixel_ratio*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai, 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

QNAS_IntErr = []
for o in range(24):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

#%%
A1_QNAS = []

o = 0
for n in range(24):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS))*-wav_pixel_ratio)
    A1_QNAS.append(A1)


Velocity_QNAS = []
c = 299792458 
for n in range(24):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(24):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(24):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])
    
INTS_Err_QNAS_N = []
INTS_Err_QNAS_P = []

for n in range(24):
    INTS_Err_QNAS_N.append(ex_INTSQNAS[n]-QNAS_IntErr[n])
    INTS_Err_QNAS_P.append(ex_INTSQNAS[n]+QNAS_IntErr[n])   
    
Limb_velocity_pix_O4 = (Period_Uranus/1000)

Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []
for a in range(12):
    Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin((a)*np.pi/32.760514148112314))

for a in range(12):
    Limb_velocity_pix_O4_test2.append(Limb_velocity_pix_O4*np.sin((-1*a)*np.pi/32.760514148112314))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
    
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
SOUTH_SLIT = np.load('SOUTH_SLIT_CONFIG.npy', allow_pickle=True)
South_Pos = []
for n in range(23):
    South_Pos.append(SOUTH_SLIT[3][n] - 180)
Approx_Longs = South_Pos
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '12:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QNAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QNAS_N[0:23]), np.flip(Velocity_Err_QNAS_P[0:23]), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQNAS[0:23]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QNAS_N[0:23]), np.flip(INTS_Err_QNAS_P[0:23]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
#ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[0:49]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNASE.append(np.flip(Velocity_QNAS[0:23]))
Vel_UranusNASE_Err.append(np.flip(Velocity_QNAS_Err[0:23])) # 1. Clean the Data!!
Int_UranusNASE.append(np.flip(ex_INTSQNAS[0:23]))           # 2. Make averages of North, South and Equator!
Int_UranusNASE_Err.append(np.flip(QNAS_IntErr[0:23]))       # 3. Construct PRF plots

#%% Now we try to implement the velocities and Ints as an average together

#Tot_Weights = [8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8, 8]
WeightE = 4
WeightN = 2
WeightS = 2
WeightA = 2

Vels_NASE1 = np.sum([Vel_UranusNASE[0:2]], axis=1)/2
#Vels_NASE1_Err = np.sum([Vel_UranusNASE_Err[0:2]], axis=1)/2
Vels_NASE1_Err = np.sqrt((Vel_UranusNASE[0] - Vels_NASE1)**2 + (Vel_UranusNASE[1] - Vels_NASE1)**2)
Int_NASE1 = np.sum([Int_UranusNASE[0:2]], axis=1)/2
#Int_NASE1_Err = np.sum([Int_UranusNASE_Err[0:2]], axis=1)/2
Int_NASE1_Err = np.sqrt((Int_UranusNASE[0] - Int_NASE1)**2 + (Int_UranusNASE[1] - Int_NASE1)**2)

Vels_NASE2 = np.sum([Vel_UranusNASE[2:4]], axis=1)/2
Vels_NASE2_Err = np.sqrt((Vel_UranusNASE[2] - Vels_NASE2)**2 + (Vel_UranusNASE[3] - Vels_NASE2)**2)
#Vels_NASE2_Err = np.sum([Vel_UranusNASE_Err[2:4]], axis=1)/2
Int_NASE2 = np.sum([Int_UranusNASE[2:4]], axis=1)/2
#Int_NASE2_Err = np.sum([Int_UranusNASE_Err[2:4]], axis=1)/2
Int_NASE2_Err = np.sqrt((Int_UranusNASE[2] - Int_NASE2)**2 + (Int_UranusNASE[3] - Int_NASE2)**2)

Vels_NASN = np.sum([Vel_UranusNASN], axis=1)/WeightN
Vels_NASN_Err = np.sqrt((Vel_UranusNASN[0] - Vels_NASN)**2 + (Vel_UranusNASN[1] - Vels_NASN)**2)
#Vels_NASN_Err = np.sum([Vel_UranusNASN_Err], axis=1)/WeightN
Int_NASN = np.sum([Int_UranusNASN], axis=1)/WeightN
#Int_NASN_Err = np.sum([Int_UranusNASN_Err], axis=1)/WeightN
Int_NASN_Err = np.sqrt((Int_UranusNASN[0] - Int_NASN)**2 + (Int_UranusNASN[1] - Int_NASN)**2)

Vels_NASS = np.sum([Vel_UranusNASS], axis=1)/WeightS
Vels_NASS_Err = np.sqrt((Vel_UranusNASS[0] - Vels_NASS)**2 + (Vel_UranusNASS[1] - Vels_NASS)**2)
#Vels_NASS_Err = np.sum([Vel_UranusNASS_Err], axis=1)/WeightS
Int_NASS = np.sum([Int_UranusNASS], axis=1)/WeightS
Int_NASS_Err = np.sqrt((Int_UranusNASS[0] - Int_NASS)**2 + (Int_UranusNASS[1] - Int_NASS)**2)
#Int_NASS_Err = np.sum([Int_UranusNASS_Err], axis=1)/WeightS

Vels_NASA = np.sum([Vel_UranusNASA], axis=1)/WeightA
Vels_NASA_Err = np.sqrt((Vel_UranusNASA[0] - Vels_NASA)**2 + (Vel_UranusNASA[1] - Vels_NASA)**2)
#Vels_NASA_Err = np.sum([Vel_UranusNASA_Err], axis=1)/WeightA
Int_NASA = np.sum([Int_UranusNASA], axis=1)/WeightA
Int_NASA_Err = np.sqrt((Int_UranusNASA[0] - Int_NASA)**2 + (Int_UranusNASA[1] - Int_NASA)**2)
#Int_NASA_Err = np.sum([Int_UranusNASA_Err], axis=1)/WeightA

Perr_Int_AURORA = []
Nerr_Int_AURORA = []
Perr_Vel_AURORA = []
Nerr_Vel_AURORA= []

for n in range(25):
    Perr_Int_AURORA.append(Int_NASA[0][n] + Int_NASA_Err[0][n])
    Nerr_Int_AURORA.append(Int_NASA[0][n] - Int_NASA_Err[0][n])
    Perr_Vel_AURORA.append(Vels_NASA[0][n] + Vels_NASA_Err[0][n])
    Nerr_Vel_AURORA.append(Vels_NASA[0][n] - Vels_NASA_Err[0][n])
    
Perr_Int_AURORA1 = []
Nerr_Int_AURORA1 = []
Perr_Vel_AURORA1 = []
Nerr_Vel_AURORA1= []

for n in range(25):
    Perr_Int_AURORA1.append(Int_NASS[0][n] + Int_NASS_Err[0][n])
    Nerr_Int_AURORA1.append(Int_NASS[0][n] - Int_NASS_Err[0][n])
    Perr_Vel_AURORA1.append(Vels_NASS[0][n] + Vels_NASS_Err[0][n])
    Nerr_Vel_AURORA1.append(Vels_NASS[0][n] - Vels_NASS_Err[0][n])

Longs = np.load('LongitudesS2.npy')
Approx_Longs = []

for o in range(25):
    Approx_Longs.append((Longs[o]+Longs[o+1])/2)
    
Approx_LongsS = Approx_Longs
    
Planet_rotation_O4 = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(1000*13.098297869601417)
    
fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Int_NASS[0], color='b', label='8:59 UTC')
#ax.fill_between(Approx_Longs, Perr_Int_AURORA1, Nerr_Int_AURORA1, color='b', alpha=0.5)
#ax.plot(Approx_Longs, Int_NASA[0], color='r', label='11:57 UTC')
ax.fill_between(Approx_Longs, Perr_Int_AURORA, Nerr_Int_AURORA, color='r', alpha=0.5)
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(0, 0.75)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

fig, ax = plt.subplots(figsize=(10,8))
ax.plot(Approx_Longs, Vels_NASS[0], color='b', label='8:59 UTC')
#ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
#ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
ax.fill_between(Approx_Longs, Perr_Vel_AURORA1, Nerr_Vel_AURORA1, color='b', alpha=0.5)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation 2016', lw = 3)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '08:59 and 11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

plt.figure(figsize=(16,18))
plt.plot(Approx_Longs, Vels_NASS[0] - np.flip(Planet_rotation_O4), color='b', label='8:59 UTC')
plt.fill_between(Approx_Longs, Perr_Vel_AURORA - np.flip(Planet_rotation_O4), Nerr_Vel_AURORA - np.flip(Planet_rotation_O4), color='r', alpha=0.5)
#plt.plot(Approx_Longs, Vels_NASA[0] - np.flip(Planet_rotation_O4), color='r', label='11:57 UTC')
#plt.fill_between(Approx_Longs, Perr_Vel_AURORA1 - np.flip(Planet_rotation_O4), Nerr_Vel_AURORA1 - np.flip(Planet_rotation_O4), color='b', alpha=0.5)
plt.plot(Approx_Longs, np.zeros(25), color='k', ls = '--', label='Planetary Rotation 2016', lw = 5)
plt.xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
plt.xlim(-90, 90) #0, 3601
plt.ylim(-2.5, 2.5)
#ax2.set_ylim(0, 0.75)
plt.title('~1 hr exposure around ' + '08:59 and 11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
plt.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plt.grid(alpha = 0.7, linestyle = '--')

#Now we attempt to fit another rotation onto this
def poly_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

def lin_func(x, a , b):
    y = a* x + b
    return y

VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(1, 26, 1), np.flip(Vels_NASS[0]), p0 = [1, -2])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(1, 26), *VelS_popt)
Linefit2 = lin_func(np.arange(1, 26), *Pos_parms)
Linefit3 = lin_func(np.arange(1, 26), *Neg_parms)

Max_Super_rotation = (Pos_parms[0]/((2.325267499446119*2)/25)*100) -100 #(%)
Min_Super_rotation = (Neg_parms[0]/((2.325267499446119*2)/25)*100) -100 #(%)
print(Min_Super_rotation, Max_Super_rotation, '08:59')

Planet_rotations = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(1000*13.098297869601417)
      
fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(1, 26, 1), Vels_NASS[0], 'ko')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit), color = 'b', label = 'Best fit')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit2), color = 'r', label = 'Best fit + error')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit3), color = 'g', label = 'Best fit - error')
ax.plot(np.arange(1, 26, 1), np.flip(Planet_rotations), 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(0, 26) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '8:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

#%%
VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(1, 26, 1), np.flip(Vels_NASA[0]), p0 = [-1, 0])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(1, 26), *VelS_popt)
Linefit2 = lin_func(np.arange(1, 26), *Pos_parms)
Linefit3 = lin_func(np.arange(1, 26), *Neg_parms)

Max_Super_rotation = (Pos_parms[0]/((2.325267499446119*2)/25)*100) - 100#(%)
Min_Super_rotation = (Neg_parms[0]/((2.325267499446119*2)/25)*100) - 100 #(%)
print(Min_Super_rotation, Max_Super_rotation, '11:57')
      
fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(1, 26, 1), Vels_NASA[0], 'ko')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit), color = 'b', label = 'Best fit')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit2), color = 'r', label = 'Best fit + error')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit3), color = 'g', label = 'Best fit - error')
ax.plot(np.arange(1, 26, 1), np.flip(Planet_rotations), 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(0, 26) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Vels_NASA[0], color='b', label='11:57 UTC')
ax.fill_between(Approx_Longs, Vels_NASA[0]+Vels_NASA_Err[0], Vels_NASA[0]-Vels_NASA_Err[0], color='b', alpha=0.5)
# ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
# ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.set_title('~1 hr exposure around ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.grid(alpha = 0.7, linestyle = '--')

#%%
Approx_Longs = E_Longs

Planet_rotations = (np.arange(-12.5, 13.5, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(1000*13.098297869601417)
Planet_rotation_O5 = Planet_rotations

VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(1, 27, 1), np.flip(Vels_NASE1[0]), p0 = [-1, 0])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(1, 27), *VelS_popt)
Linefit2 = lin_func(np.arange(1, 27), *Pos_parms)
Linefit3 = lin_func(np.arange(1, 27), *Neg_parms)

Max_Super_rotation = (Pos_parms[0]/((2.325267499446119*2)/26)*100)-100 #(%) 
Min_Super_rotation = (Neg_parms[0]/((2.325267499446119*2)/26)*100)-100 #(%)
print(Min_Super_rotation, Max_Super_rotation, '08:00')
      
fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(1, 27, 1), Vels_NASE1[0], 'ko')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit), color = 'b', label = 'Best fit')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit2), color = 'r', label = 'Best fit + error')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit3), color = 'g', label = 'Best fit - error')
ax.plot(np.arange(1, 27, 1), np.flip(Planet_rotations), 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(0, 27) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '08:00' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:27], Vels_NASE1[0], color='b', label='8:00 UTC')
ax.fill_between(Approx_Longs[1:27], Vels_NASE1[0]+Vels_NASE1_Err[0], Vels_NASE1[0]-Vels_NASE1_Err[0], color='b', alpha=0.5)
# ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
# ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs[1:27], np.flip(Planet_rotations), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '08:00' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

#%%
VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(1, 27, 1), np.flip(Vels_NASE2[0]), p0 = [-1, 0])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(1, 27), *VelS_popt)
Linefit2 = lin_func(np.arange(1, 27), *Pos_parms)
Linefit3 = lin_func(np.arange(1, 27), *Neg_parms)

Max_Super_rotation = (Pos_parms[0]/((2.3479227476524116*2)/26)*100)-100 #(%)
Min_Super_rotation = (Neg_parms[0]/((2.3479227476524116*2)/26)*100)-100 #(%)
print(Min_Super_rotation, Max_Super_rotation, '10:58')

Planet_rotations = (np.arange(-13, 14, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(1000*13.098297869601417)
Planet_rotations_O6 = Planet_rotations      

fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(1, 27, 1), Vels_NASE2[0], 'ko')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit), color = 'b', label = 'Best fit')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit2), color = 'r', label = 'Best fit + error')
ax.plot(np.arange(1, 27, 1), np.flip(Linefit3), color = 'g', label = 'Best fit - error')
ax.plot(np.arange(1, 27, 1), np.flip(Planet_rotations[0:26]), 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(0, 27) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:27], Vels_NASE2[0], color='b', label='10:58 UTC')
ax.fill_between(Approx_Longs[1:27], Vels_NASE2[0]+Vels_NASE2_Err[0], Vels_NASE2[0]-Vels_NASE2_Err[0], color='b', alpha=0.5)
# ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
# ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs[1:27], np.flip(Planet_rotations[0:26]), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

plt.figure(figsize=(16,18))
plt.plot(Approx_Longs[1:27], Vels_NASE1[0] - np.flip(Planet_rotations[0:26]), color='b', label='8:00 UTC')
plt.fill_between(Approx_Longs[1:27], Vels_NASE1[0]+Vels_NASE1_Err[0] - np.flip(Planet_rotations[0:26]), Vels_NASE1[0]-Vels_NASE1_Err[0] - np.flip(Planet_rotations[0:26]), color='b', alpha=0.5)
plt.plot(Approx_Longs[1:27], Vels_NASE2[0] - np.flip(Planet_rotations[0:26]), color='r', label='10:58 UTC')
plt.fill_between(Approx_Longs[1:27], Vels_NASE2[0]+Vels_NASE2_Err[0] - np.flip(Planet_rotations[0:26]), Vels_NASE2[0]-Vels_NASE2_Err[0] - np.flip(Planet_rotations[0:26]), color='r', alpha=0.5)
plt.plot(Approx_Longs[1:27], np.zeros(26), color='k', ls = '--', label='Planetary Rotation 2016', lw = 5)
plt.xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
plt.xlim(-90, 90) #0, 3601
plt.ylim(-2.5, 2.5)
#ax2.set_ylim(0, 0.75)
plt.title('~1 hr exposure around ' + '08:00 and 10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
plt.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plt.grid(alpha = 0.7, linestyle = '--')

#%%
Approx_Longs = N_Longs
    
Approx_LongsN = Approx_Longs

VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(1, 26, 1), np.flip(Vels_NASN[0]), p0 = [-1, 0])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(1, 26), *VelS_popt)
Linefit2 = lin_func(np.arange(1, 26), *Pos_parms)
Linefit3 = lin_func(np.arange(1, 26), *Neg_parms)

Max_Super_rotation = (Pos_parms[0]/((2.3479227476524116*2)/25)*100)-100 #(%)
Min_Super_rotation = (Neg_parms[0]/((2.3479227476524116*2)/25)*100)-100 #(%)
print(Min_Super_rotation, Max_Super_rotation, '09:59')
      
fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(1, 26, 1), Vels_NASN[0], 'ko')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit), color = 'b', label = 'Best fit')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit2), color = 'r', label = 'Best fit + error')
ax.plot(np.arange(1, 26, 1), np.flip(Linefit3), color = 'g', label = 'Best fit - error')
ax.plot(np.arange(1, 26, 1), np.flip(Planet_rotations[1:26]), 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(0, 26) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Vels_NASN[0], color='b', label='9:59 UTC')
ax.fill_between(Approx_Longs, Vels_NASN[0]+Vels_NASN_Err[0], Vels_NASN[0]-Vels_NASN_Err[0], color='b', alpha=0.5)
# ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
# ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs, np.flip(Planet_rotations[1:26]), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '9:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

plt.figure(figsize=(16,18))
plt.plot(Approx_Longs, Vels_NASN[0] - np.flip(Planet_rotations[1:26]), color='b', label='09:59 UTC')
plt.fill_between(Approx_Longs, Vels_NASN[0]+Vels_NASN_Err[0] - np.flip(Planet_rotations[1:26]), Vels_NASN[0]-Vels_NASN_Err[0] - np.flip(Planet_rotations[1:26]), color='b', alpha=0.5)
#plt.plot(Approx_Longs, Vels_NASE2[0] - np.flip(Planet_rotations[0:26]), color='r', label='10:58 UTC')
#plt.fill_between(Approx_Longs, Vels_NASE2[0]+Vels_NASE2_Err[0] - np.flip(Planet_rotations[0:26]), Vels_NASE2[0]-Vels_NASE2_Err[0] - np.flip(Planet_rotations[0:26]), color='r', alpha=0.5)
plt.plot(Approx_Longs, np.zeros(25), color='k', ls = '--', label='Planetary Rotation 2016', lw = 5)
plt.xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
#plt.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
plt.xlim(-90, 90) #0, 3601
plt.ylim(-2.5, 2.5)
#ax2.set_ylim(0, 0.75)
plt.title('~1 hr exposure around ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
plt.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
plt.grid(alpha = 0.7, linestyle = '--')

#%%
fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Vels_NASS[0], color='b', label='8:59 UTC')
ax.fill_between(Approx_Longs, Perr_Vel_AURORA1, Nerr_Vel_AURORA1, color='b', alpha=0.5)
ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')


fig, ax = plt.subplots(figsize=(10,8))
#ax2 = ax.twinx()
# ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, Vels_NASS[0], color='b', label='8:59 UTC')
ax.fill_between(Approx_Longs, Perr_Vel_AURORA1, Nerr_Vel_AURORA1, color='b', alpha=0.5)
ax.plot(Approx_Longs, Vels_NASA[0], color='r', label='11:57 UTC')
ax.fill_between(Approx_Longs, Perr_Vel_AURORA, Nerr_Vel_AURORA, color='r', alpha=0.5)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', linewidth = 3)
#ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ax.grid(alpha = 0.7, linestyle = '--')

#Construct an average Int to minus from to make Mapping easier to see!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#%% Now lets do the planetary reference frame

