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

Wavelength_O4 = np.linspace(3.932951, 3.966128, 2048)
Wavelength_O5 = np.linspace(3.963083, 3.997173, 2048)
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

# plt.figure()
# plt.imshow(np.fliplr(np.flipud(image_data_alt)), cmap='gist_gray')
# plt.xlabel(r'Spectral pixel position (Pixel No.) ', fontsize=16)
# #label_x = 0, 100, 200, 300, 400, 500, 600, 700
# #plt.xticks(label_x, ('3.9506', '3.9564', '3.9622', '3.9680', '3.9738', '3.9796', '3.9854', '3.9919'), fontsize=15)
# plt.ylabel('Spatial position across Slit (Pixel No.)', fontsize=16)
# #plt.axvline(833, color='red', ls='-.', alpha=0.125)
# plt.annotate(r'Q(1,0${^-}$)', (862, 395), color='red', fontsize=12)
# #plt.axvline(688, color='blue', ls='-.', alpha=0.125)
# plt.annotate(r'Q(3,0${^-}$)', (717, 285), color='blue', fontsize=12)
# #plt.axhline(40, color='r', ls='--')
# plt.title(r'Single AB of Uranus with ISHELL on $11^{th}$ October 2016', fontsize=23)
# cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
# cbar.set_label(r'Intensity ($Wm^{-2}sr^{-1}$)', fontsize=20)

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
A_average = -4.1696909229209164e-05
B_average = -0.04916341231347669

No_list = 14, 17, 18, 21, 22, 25, 26, 29, 30, 33, 34, 37, 38, 41, 42, 45, 46, 49, 50, 53, 54, 57, 58, 61, 62, 65, 66, 69, 70, 73, 74, 77, 78, 81, 82, 85, 86, 89, 90, 93, 94, 97, 98, 101, 102
s = 0 #This will be used to create the ABBA pattern
All_IRTF_Data_Q1 = []
for n in range(45): #We use this list to create a list which holds all the data from Order19
    if n < 43:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.0' + str(No_list[n]) + '/uranus.0' + str(No_list[n]) + '.order.04.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        IRTF_Data = image_datai
        IRTF_Data[IRTF_Data == -inf] = 0
        IRTF_Data[IRTF_Data == inf] = 0
        IRTF_Data[np.isnan(IRTF_Data)] = 0
        All_IRTF_Data_Q1.append(IRTF_Data)
    else:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Uranus2/uranus.' + str(No_list[n]) + '/uranus.' + str(No_list[n]) + '.order.04.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        IRTF_Data = image_datai
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
A_average = 0.0001676459795899812
B_average = 0.008931745988305165

for n in range(45):
    if n == 0:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q1[0]
            Spec_row = Fig[i, :]
            ii = i - Approx_Mid
            corrective_shift_factor = ((ii**2)*A_average + ii*B_average)*-1 #Not sure why this isn't working
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
            ii = i - Approx_Mid
            corrective_shift_factor = ((ii**2)*A_average + ii*B_average)*-1 #Not sure why this isn't working
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
        
A_average = 0.0001676459795899812
B_average = 0.008931745988305165
IRTF_Data_Q3 = []
Approx_Mid = (90-6)/2 + 6
S_point = All_IRTF_Data_Q3[0]

for n in range(45):
    if n == 0:
        for i in range(len(S_point[:,0])):
            Fig = All_IRTF_Data_Q3[0]
            Spec_row = Fig[i, :]
            ii = i - Approx_Mid
            corrective_shift_factor = ((ii**2)*A_average + ii*B_average)*-1 #Not sure why this isn't working
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
            ii = i - Approx_Mid
            corrective_shift_factor = ((ii**2)*A_average + ii*B_average)*-1 #Not sure why this isn't working
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
        
# plt.figure()
# plt.imshow(Total_IRTF_DataQ3, cmap='gist_gray', vmax = 2*10**-1, vmin=-2*10**-1)
        
#%%
h3p = h3ppy.h3p()
wave2 = h3p.wavegen(np.nanmin(Wavelength_O4), np.nanmax(Wavelength_O4), 2048)
Total_IRTF_Data_Alt = Total_IRTF_Data/90
Total_IRTF_Data_Alt = acre(Total_IRTF_Data_Alt, width = 15, verbose = False)
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt > 0.025] = 0
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt < -0.025] = 0

plt.figure()
plt.imshow(Total_IRTF_Data, cmap='gist_gray', vmax = 0.1, vmin = -0.1)
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
plt.imshow(Total_IRTF_DataQ3, cmap='gist_gray', vmax = 0.1, vmin = -0.1)
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
    datai = dataI[15+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmin(datai[1347:1367]))
    a1_pointP = a1_pointP[0][0] - 1250
    ex_resultQ1P = gmodel.fit(datai[1250:1450], x=np.arange(200), a0=np.nanmin(datai[1347:1367]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    Pos_Q3Tn.append(ex_pQ1P.a1 + 1250)
    Err_Q3Tn.append(ex_eQ1P[1])
    datai = dataI[56+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1347:1367]))
    a1_pointP = a1_pointP[0][0] - 1250
    ex_resultQ1P = gmodel.fit(datai[1250:1450], x=np.arange(200), a0=np.nanmax(datai[1347:1367]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
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
    dataI = Total_IRTF_DataQ3
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
h3p = h3ppy.h3p()
wave2 = h3p.wavegen(np.nanmin(Wavelength_O4), np.nanmax(Wavelength_O4), 2048)

#%% Here use the advantage if the star seperation to help here
T_IRTFData = Total_IRTF_Data_Alt

Order_4_Offset = np.ones(2048)*41.2337 #These are in microns
Order_5_Offset = np.ones(2048)*41.6508

#Lets start shifting the positive line onto the negative line looping through all of x
Shape_Order4 = np.shape(IRTF_Data_Q1[:,:,0])
Combined_Q1 = []

for x in range(45):
    IRTFData = IRTF_Data_Q1[:,:,x]
    IRTFData = acre(IRTFData, width = 15, verbose = False)
    for n in range(Shape_Order4[1]):
        Spec_col = IRTFData[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_4_Offset[n]], mode='wrap')
        if n == 0:
            Order_4_shift = shifted_Spec_col
        else:
            Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))
    Order_4_shiftv1 = np.flipud(np.rot90(Order_4_shift))
    #plt.imshow(Order_4_shiftv1, cmap='gist_gray', vmax = 0.025, vmin = -0.025)
    for n in range(Shape_Order4[0]):
        Spec_col = Order_4_shiftv1[n,:]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*np.nanmean(Olap)], mode='wrap')
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
    IRTFData = acre(IRTFData, width = 15, verbose = False)
    for n in range(Shape_Order5[1]):
        Spec_col = IRTFData[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_5_Offset[n]], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shiftv1 = np.flipud(np.rot90(Order_5_shift))
    for n in range(Shape_Order4[0]):
        Spec_col = Order_5_shiftv1[n,:]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*np.nanmean(Olapv3)], mode='wrap')
        if n == 0:
            Order_5_shift = shifted_Spec_col
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))
    Order_5_shift = Order_5_shift
    Combined_Q3.append(Order_5_shift)

#%%
AB_Combined_Order4 = []
AB_Combined_Order5 = []

for x in range(45):
    IRTFData_Q1 = Combined_Q1[x] - IRTF_Data_Q1[:,:,x]
    AB_Combined_Order4.append(IRTFData_Q1)

for x in range(45):
    IRTFData_Q3 = Combined_Q3[x] - IRTF_Data_Q3[:,:,x]
    AB_Combined_Order5.append(IRTFData_Q3)
    
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
    dataI[dataI > 0.020] = 0.0025
    Q1_Data_Capture.append(dataI[0:50, :])
    o += 4
    # if x == 8:
    #     for m in range(25):
    #         datai = dataI[15+m, :]
    #         gmodel = Model(gauss_fit)
    #         a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    #         a1_pointP = a1_pointP[0][0] - 1110
    #         ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    #         ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    #         ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    #         A0 = ex_pQ1P.a0
    #         ex_A0_Q18.append(A0)
    #         A2 = ex_pQ1P.a2
    #         ex_A2_Q18.append((A2)*1.6199707031249897e-05)
    #         err_A0 = ex_eQ1P[0]
    #         ex_Err_A0_Q18.append(err_A0)
    #         err_A2 = ex_eQ1P[2]
    #         ex_Err_A2_Q18.append(err_A2*1.6199707031249897e-05)
    #         ex_A1_Q18.append(ex_pQ1P.a1)
    #         ex_Err_A1_Q18.append(ex_eQ1P[1]*1.6199707031249897e-05)
    #         ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
    #         ex_FWHMQ18.append(ex_FWHM)
    #         ex_INTSQ18.append(A0*ex_FWHM*(10**6))
    # if x == 3:
    #     for m in range(22):
    #         datai = dataI[16+m, :]
    #         gmodel = Model(gauss_fit)
    #         a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    #         a1_pointP = a1_pointP[0][0] - 1110
    #         ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    #         ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    #         ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    #         A0 = ex_pQ1P.a0
    #         ex_A0_Q18.append(A0)
    #         A2 = ex_pQ1P.a2
    #         ex_A2_Q18.append((A2)*1.6199707031249897e-05)
    #         err_A0 = ex_eQ1P[0]
    #         ex_Err_A0_Q18.append(err_A0)
    #         err_A2 = ex_eQ1P[2]
    #         ex_Err_A2_Q18.append(err_A2*1.6199707031249897e-05)
    #         ex_A1_Q18.append(ex_pQ1P.a1)
    #         ex_Err_A1_Q18.append(ex_eQ1P[1]*1.6199707031249897e-05)
    #         ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
    #         ex_FWHMQ18.append(ex_FWHM)
    #         ex_INTSQ18.append(A0*ex_FWHM*(10**6))
    #         # plt.figure()
    #         # plt.plot(XX, datai[1110:1310], 'o')
    #         # plt.plot(XX, ex_resultQ1P.best_fit, '-')
    # for m in range(22):
    #     datai = dataI[16+m, :]
    #     gmodel = Model(gauss_fit)
    #     a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    #     a1_pointP = a1_pointP[0][0] - 1110
    #     ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    #     ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    #     ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    #     A0 = ex_pQ1P.a0
    #     ex_A0_Q1.append(A0)
    #     A2 = ex_pQ1P.a2
    #     ex_A2_Q1.append((A2)*1.6199707031249897e-05)
    #     err_A0 = ex_eQ1P[0]
    #     ex_Err_A0_Q1.append(err_A0)
    #     err_A2 = ex_eQ1P[2]
    #     ex_Err_A2_Q1.append(err_A2*1.6199707031249897e-05)
    #     ex_A1_Q1.append(ex_pQ1P.a1)
    #     ex_Err_A1_Q1.append(ex_eQ1P[1]*1.6199707031249897e-05)
    #     ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
    #     ex_FWHMQ1.append(ex_FWHM)
    #     ex_INTSQ1.append(A0*ex_FWHM*(10**6))
    #     Pos_Q1.append(ex_pQ1P.a1 + 1110)
    #     # if x == 9:
    #     #     plt.figure()
    #     #     plt.plot(XX, datai[1110:1310], 'bo')   
    #     #     plt.plot(XX, ex_resultQ1P.best_fit, '--')
    #     #     print(ex_pQ1P.a1)
      
# Q1_IntErr = []
# #Slit is 15 so

# for o in range(242):
#     Q1IntErr = ex_INTSQ1[o]*np.sqrt((ex_Err_A0_Q1[o]/ex_A0_Q1[o])**2 + (ex_Err_A2_Q1[o]/ex_A2_Q1[o])**2)
#     Q1_IntErr.append(Q1IntErr)

o = 0
for x in range(11):
    dataI = (AB_Combined_Order5[o]+AB_Combined_Order5[o+1]+AB_Combined_Order5[o+2]+AB_Combined_Order5[o+3])/8
    dataI[dataI < -0.004] = 0
    Q3_Data_Capture.append(dataI[0:50, :])
    o += 4
#     for m in range(22):
#         datai = dataI[16+m, :]
#         gmodel = Model(gauss_fit)
#         a1_pointP = np.where(datai == np.nanmax(datai[1350:1370]))
#         a1_pointP = a1_pointP[0][0] - 1260
#         ex_resultQ3P = gmodel.fit(datai[1260:1460], x=XX, a0=np.nanmax(datai[1350:1370]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
#         ex_pQ3P = SimpleNamespace(**ex_resultQ3P.best_values)
#         ex_eQ3P = np.sqrt(np.diag(ex_resultQ3P.covar))
#         A0 = ex_pQ3P.a0
#         ex_A0_Q3.append(A0)
#         A2 = ex_pQ3P.a2
#         ex_A2_Q3.append((A2)*1.6199707031249897e-05)
#         err_A0 = ex_eQ1P[0]
#         ex_Err_A0_Q3.append(err_A0)
#         err_A2 = ex_eQ1P[2]
#         ex_Err_A2_Q3.append(err_A2*1.6199707031249897e-05)
#         ex_A1_Q3.append(ex_pQ3P.a1)
#         ex_Err_A1_Q3.append(ex_eQ3P[1]*1.6199707031249897e-05)
#         ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
#         ex_FWHMQ3.append(ex_FWHM)
#         ex_INTSQ3.append(A0*ex_FWHM*(10**6)) #Needs work with Q3 to skip difficult pixels?
#         Pos_Q3.append(ex_pQ3P.a1 + 1260)
#         # if x == 3 or x == 2:
#         #     if m > 5 and m < 15:
#         #         plt.figure()
#         #         plt.plot(XX, datai[1260:1460], 'bo')
#         #         z = (XX - ex_pQ3P.a1) / ex_pQ3P.a2
#         #         y = ex_pQ3P.a0 * np.exp(-z**2 / ex_pQ3P.a2) + ex_pQ3P.a3 + ex_pQ3P.a4 * XX + ex_pQ3P.a5 * XX**2
#         #         plt.plot(XX, y, 'g')
#         #         print(A0*ex_FWHM*(10**6))       
      
# Q3_IntErr = []
# #Slit is 15 so

# for o in range(242):
#     Q3IntErr = ex_INTSQ3[o]*np.sqrt((ex_Err_A0_Q3[o]/ex_A0_Q3[o])**2 + (ex_Err_A2_Q3[o]/ex_A2_Q3[o])**2)
#     Q3_IntErr.append(Q3IntErr)

#%% Now we try to pair Q3 and Q1 together 

for n in range(11):
    if n == 0:
        Total_IRTF_DataQ1 = Q1_Data_Capture[0]
    else:
        Total_IRTF_DataQ1 = np.add(Total_IRTF_DataQ1, Q1_Data_Capture[n])

for n in range(11):
    if n == 0:
        Total_IRTF_DataQ3 = Q3_Data_Capture[0]
        #Total_IRTF_DataQ3 = zoom(Q3_Data_Capture[0], (1.0, 0.9697589886378948), mode='constant')
    else:
        #Testing = zoom(Q3_Data_Capture[n], (1.0, 0.9697589886378948), mode='constant')
        #Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, Testing)
        Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, Q3_Data_Capture[n])

    #%%
Total_IRTF_Data_Alt = acre(Total_IRTF_DataQ3, width = 15, verbose = False)
Total_IRTF_Data_Alt2 = acre(Total_IRTF_DataQ1, width = 15, verbose = False)

ex_A1_Q1T = []
ex_Err_A1_Q1T = []
Pos_Q1T = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ1
    datai = dataI[15+m, :]
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
    datai = dataI[15+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[458:478]))
    a1_pointP = a1_pointP[0][0] - 440
    ex_resultQ2P = gmodel.fit(datai[440:490], x=np.arange(50), a0=np.nanmax(datai[458:478]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ2P = SimpleNamespace(**ex_resultQ2P.best_values)
    ex_eQ2P = np.sqrt(np.diag(ex_resultQ2P.covar))
    Pos_Q2.append(ex_pQ2P.a1 + 440)
    
ex_A1_Q3T = []
ex_Err_A1_Q3T = []
Pos_Q3T = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ3
    datai = dataI[15+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1350:1370]))
    a1_pointP = a1_pointP[0][0] - 1290
    ex_resultQ3P = gmodel.fit(datai[1290:1390], x=np.arange(100), a0=np.nanmax(datai[1350:1370]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ3P = SimpleNamespace(**ex_resultQ3P.best_values)
    ex_eQ3P = np.sqrt(np.diag(ex_resultQ3P.covar))
    Pos_Q3T.append(ex_pQ3P.a1 + 1290)

ex_A1_Q31 = []
ex_Err_A1_Q31 = []
Pos_Q31 = []

o = 0
for m in range(25):
    dataI = Total_IRTF_DataQ3
    datai = dataI[15+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1445:1465]))
    a1_pointP = a1_pointP[0][0] - 1420
    ex_resultQ31P = gmodel.fit(datai[1420:1480], x=np.arange(60), a0=np.nanmax(datai[1445:1465]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ31P = SimpleNamespace(**ex_resultQ31P.best_values)
    ex_eQ31P = np.sqrt(np.diag(ex_resultQ31P.covar))
    Pos_Q31.append(ex_pQ31P.a1 + 1420)
    
shift_list = []
shift_list2 = []
shift_list3 = []

for n in range(25):
    shift_list.append(Pos_Q3T[n] - Pos_Q1T[n])
    shift_list2.append(Pos_Q2[n] - Pos_Q1T[n])
    shift_list3.append(Pos_Q31[n] - Pos_Q1T[n])
    
shift_amount1 = np.nanmean(shift_list[0:23])
shift_amount2 = np.nanmean(shift_list2[0:23])
shift_amount3 = np.nanmean(shift_list3[0:23])

#%%
ii = 0
CombinedImage = []

for nn in range(11):
    Q1Wave = Q1_Data_Capture[nn]
    #Q3Wave = Total_IRTF_DataQ3 = zoom(Q3_Data_Capture[nn], (1.0, 0.9697589886378948), mode='constant')
    Q3Wave = Q3_Data_Capture[nn]
    for n in range(35):
        Shifted_Q3 = Q3Wave[n+10,:]
        shift_amount = shift_amount1
        shift_amount2 = shift_amount2
        shift_amount3 = shift_amount3
        shifted_Q3 = f.shift(Shifted_Q3, shift=[-shift_amount], mode='wrap') #41.6508 for Order 5
        shifted_Q2 = f.shift(Shifted_Q3, shift=[-shift_amount2], mode='wrap')
        shifted_Q31 = f.shift(Shifted_Q3, shift=[-shift_amount3], mode='wrap')
        if n == 0:
            Order_5_shift = np.vstack((Q3Wave[0:10, :], shifted_Q3))
            Order_5_shift1 = np.vstack((Q3Wave[0:10, :], shifted_Q2))
            Order_5_shift2 = np.vstack((Q3Wave[0:10, :], shifted_Q31))
        elif n == 34:
            Order_5_shift = np.vstack((Order_5_shift, Q3Wave[44:51, :]))
            Order_5_shift1 = np.vstack((Order_5_shift1, Q3Wave[44:51, :]))
            Order_5_shift2 = np.vstack((Order_5_shift2, Q3Wave[44:51, :]))
        else:
            Order_5_shift = np.vstack((Order_5_shift, shifted_Q3))
            Order_5_shift1 = np.vstack((Order_5_shift1, shifted_Q2))
            Order_5_shift2 = np.vstack((Order_5_shift2, shifted_Q31))
        ii += 1
    ComboFigure  = (Q1Wave + Order_5_shift + Order_5_shift1)/3
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
    datai = dataI[15+m, :]
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

for n in range(11):
    dataI = CombinedImage[n]
    dataI[dataI > 0.025] = 0.0025
    for m in range(25):
        datai = dataI[14+m, :]
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQP = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQP = SimpleNamespace(**ex_resultQP.best_values)
        ex_eQP = np.sqrt(np.diag(ex_resultQP.covar))
        A0 = ex_pQP.a0
        ex_A0_Q.append(A0)
        A2 = ex_pQP.a2
        ex_A2_Q.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQP[0]
        ex_Err_A0_Q.append(err_A0)
        err_A2 = ex_eQP[2]
        ex_Err_A2_Q.append(err_A2*1.4374243049645295e-05)
        ex_A1_Q.append(ex_pQP.a1+1110)
        ex_Err_A1_Q.append(ex_eQP[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQ.append(ex_FWHM)
        ex_INTSQ.append(A0*ex_FWHM*(10**6))
        # if n == 9:
        #     plt.figure()
        #     plt.plot(XX, datai[1110:1310], 'bo')
        #     plt.plot(XX, ex_resultQP.best_fit)

Q_IntErr = []
for o in range(275):
    QIntErr = ex_INTSQ[o]*np.sqrt((ex_Err_A0_Q[o]/ex_A0_Q[o])**2 + (ex_Err_A2_Q[o]/ex_A2_Q[o])**2)
    Q_IntErr.append(QIntErr)

#%% 
A1_Q = []

o = 0
for n in range(275):
    A1 = ((ex_A1_Q[n]-np.nanmean(ex_A1_Q[o:o+24]))*-1.4374243049645295e-05)
    A1_Q.append(A1)
    if (n + 1) % 25 == 0:
        o += 25

#%%
c = 299792458 
Velocity_Q = []
lamda_H3_O4 = 3.95295 
lamda_H3_O5 = 3.98558
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846

for n in range(275):
    V = ((A1_Q[n]/lamda_H3_O4)*c)
    Velocity_Q.append(V/1000)
    
Velocity_Q_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(275):
    Err = Velocity_Q[n]*np.sqrt((ex_Err_A1_Q[n]/A1_Q[n])**2)
    Velocity_Q_Err.append(np.sqrt(Err**2))
    
Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

R_Uranus = 25362*1000 #m
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

Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 25)
# Planet_rotation_O5 = np.linspace(Limb_velocity_pix_O5, -1*Limb_velocity_pix_O5, 22)

Velocity_Err_QP = []
Velocity_Err_QN = []

for n in range(275):
    Velocity_Err_QP.append(Velocity_Q[n] + Velocity_Q_Err[n])
    Velocity_Err_QN.append(Velocity_Q[n] - Velocity_Q_Err[n])
    
Ints_Err_QP = []
Ints_Err_QN = []
    
for n in range(275):
    Ints_Err_QP.append(ex_INTSQ[n] + Q_IntErr[n])
    Ints_Err_QN.append(ex_INTSQ[n] - Q_IntErr[n])    

TIME = ['07:45', '08:14', '08:44', '09:14', '09:44', '10:13', '10:43', '11:12', '11:42', '12:11', '12:42']

o = 0
for n in range(11):
    fig, ax = plt.subplots(figsize=(10,8))
    ax2 = ax.twinx()
    ax.set_title('~30 minute exposure at ' + str(TIME[n]) +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
    ax.plot(np.arange(25), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
    ax.plot(np.arange(25), np.flip(Velocity_Q[o:o+25]), color='b', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
    #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
    ax.fill_between(np.arange(25), np.flip(Velocity_Err_QN[o:o+25]), np.flip(Velocity_Err_QP[o:o+25]), color='b', alpha=0.5)
    #ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
    ax2.plot(np.arange(25), np.flip(ex_INTSQ[o:o+25]), color='r', label='IR Intensity', lw = 5)
    ax2.fill_between(np.arange(25), np.flip(Ints_Err_QP[o:o+25]), np.flip(Ints_Err_QN[o:o+25]), color='r', alpha=0.5)
    #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
    #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
    # ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
    #ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
    ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
    ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
    #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
    # ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    #ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    # ax.set_xlim(-90, 90) #0, 3601
    ax.set_ylim(-3, 3)
    ax2.set_ylim(0, 0.6)
    ax.legend(loc='upper left', fontsize=25)
    ax2.legend(loc='upper right', fontsize=25)
    o += 25

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

for n in range(0,2,1):
    dataI = CombinedImage[n+8]
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(25):
        if n == 0:
            datai = dataI[15+m, :]
        else:
            datai = dataI[15+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1200:1220]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1200:1220]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQAS.append(ex_FWHM)
        ex_INTSQAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(50):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(50):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[o+2:o+22]))*-1.4374243049645295e-05)
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
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
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
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 23)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(25), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(25), np.flip(Velocity_QAS[25:50]), color='b', ls= '--', label='Combined Q(1,0), Q(3,0) and (3,1) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(25), np.flip(Velocity_Err_QAS_N[25:50]), np.flip(Velocity_Err_QAS_P[25:50]), color='b', alpha=0.5)
ax2.plot(np.arange(25), np.flip(ex_INTSQAS[25:50]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(np.arange(25), np.flip(INTS_Err_QAS_N[25:50]), np.flip(INTS_Err_QAS_P[25:50]), color='r', alpha=0.5)
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
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(25), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(25), np.flip(Velocity_QAS[0:25]), color='b', ls= '--', label='Combined Q(1,0), Q(3,0) and (3,1) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(25), np.flip(Velocity_Err_QAS_N[0:25]), np.flip(Velocity_Err_QAS_P[0:25]), color='b', alpha=0.5)
ax2.plot(np.arange(25), np.flip(ex_INTSQAS[0:25]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(np.arange(25), np.flip(INTS_Err_QAS_N[0:25]), np.flip(INTS_Err_QAS_P[0:25]), color='r', alpha=0.5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

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

for n in range(0,2,1):
    dataI = CombinedImage[n]
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(26):
        if n == 0:
            datai = dataI[14+m, :]
        else:
            datai = dataI[14+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(52):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

A1_QNAS = []

o = 0
for n in range(52):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o+2:o+23]))*-1.4374243049645295e-05)
    A1_QNAS.append(A1)
    if n == 25:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(52):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(52):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(52):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 24)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '08:14' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[26:52]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[26:52]), np.flip(Velocity_Err_QNAS_P[26:52]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[26:52]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '07:45' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNAS = []
Vel_UranusNAS_Err = []
Int_UranusNAS = []
Int_UranusNAS_Err = []

Vel_UranusNAS.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNAS.append(np.flip(Velocity_QNAS[26:52]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[26:52]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[26:52]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[26:52]))

#%% Set 2 #2 #3 
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
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(26):
        if n == 0:
            datai = dataI[14+m, :]
        else:
            datai = dataI[14+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(52):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

A1_QNAS = []

o = 0
for n in range(52):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o+2:o+23]))*-1.4374243049645295e-05)
    A1_QNAS.append(A1)
    if n == 25:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(52):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(52):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(52):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 24)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '09:14' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[26:52]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[26:52]), np.flip(Velocity_Err_QNAS_P[26:52]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[26:52]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '08:44' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNAS.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNAS.append(np.flip(Velocity_QNAS[26:52]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[26:52]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[26:52]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[26:52]))

#%% Set 3 #4 #5
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
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(26):
        if n == 0:
            datai = dataI[14+m, :]
        else:
            datai = dataI[14+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(52):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

A1_QNAS = []

o = 0
for n in range(52):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o+2:o+23]))*-1.4374243049645295e-05)
    A1_QNAS.append(A1)
    if n == 25:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(52):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(52):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(52):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 24)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '10:13' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[26:52]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[26:52]), np.flip(Velocity_Err_QNAS_P[26:52]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[26:52]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '09:44' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNAS.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNAS.append(np.flip(Velocity_QNAS[26:52]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[26:52]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[26:52]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[26:52]))

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
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(26):
        if n == 0:
            datai = dataI[14+m, :]
        else:
            datai = dataI[14+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(52):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

A1_QNAS = []

o = 0
for n in range(52):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o+2:o+23]))*-1.4374243049645295e-05)
    A1_QNAS.append(A1)
    if n == 25:
        o = 25

Velocity_QNAS = []
c = 299792458 
for n in range(52):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(52):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(52):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 24)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:12' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[26:52]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[26:52]), np.flip(Velocity_Err_QNAS_P[26:52]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[26:52]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25) #Need to add labels for the limbs of the planet along with the errors on either side. Look into making an average to compare against

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '10:43' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNAS.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNAS.append(np.flip(Velocity_QNAS[26:52]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[26:52]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[26:52]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[0:26]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[26:52]))

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
    dataI[dataI < -0.004] = 0
    dataI[dataI > 0.01] = 0.001
    for m in range(26):
        if n == 0:
            datai = dataI[14+m, :]
        else:
            datai = dataI[14+m, :]            
        gmodel = Model(gauss_fit)
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0] - 1110
        ex_resultQPAS = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        ex_pQPAS = SimpleNamespace(**ex_resultQPAS.best_values)
        ex_eQPAS = np.sqrt(np.diag(ex_resultQPAS.covar))
        A0 = ex_pQPAS.a0
        ex_A0_QNAS.append(A0)
        A2 = ex_pQPAS.a2
        ex_A2_QNAS.append((A2)*1.4374243049645295e-05)
        err_A0 = ex_eQPAS[0]
        ex_Err_A0_QNAS.append(err_A0)
        err_A2 = ex_eQPAS[2]
        ex_Err_A2_QNAS.append(err_A2*1.4374243049645295e-05)
        ex_A1_QNAS.append(ex_pQPAS.a1 + 1110)
        ex_Err_A1_QNAS.append(ex_eQPAS[1]*1.4374243049645295e-05)
        ex_FWHM = A2*1.4374243049645295e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQNAS.append(ex_FWHM)
        ex_INTSQNAS.append(A0*ex_FWHM*(10**6))
        # if n == 0:
        #         plt.figure()
        #         plt.plot(XX, datai[1110:1310], 'bo')
        #         plt.plot(XX, ex_resultQPAS.best_fit, '-')            

for o in range(26):
    QIntErr = ex_INTSQNAS[o]*np.sqrt((ex_Err_A0_QNAS[o]/ex_A0_QNAS[o])**2 + (ex_Err_A2_QNAS[o]/ex_A2_QNAS[o])**2)
    QNAS_IntErr.append(QIntErr)

A1_QNAS = []

o = 0
for n in range(26):
    A1 = ((ex_A1_QNAS[n]-np.nanmean(ex_A1_QNAS[o+2:o+23]))*-1.4374243049645295e-05)
    A1_QNAS.append(A1)


Velocity_QNAS = []
c = 299792458 
for n in range(26):
    V = ((A1_QNAS[n]/lamda_H3_O4)*c)
    Velocity_QNAS.append(V/1000)
    
Velocity_QNAS_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(26):
    Err = Velocity_QNAS[n]*np.sqrt((ex_Err_A1_QNAS[n]/A1_QNAS[n])**2)
    Velocity_QNAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QNAS_N = []
Velocity_Err_QNAS_P = []

for n in range(26):
    Velocity_Err_QNAS_N.append(Velocity_QNAS[n]-Velocity_QNAS_Err[n])
    Velocity_Err_QNAS_P.append(Velocity_QNAS[n]+Velocity_QNAS_Err[n])

Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)
Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 24)
Planet_rotation_O4 = np.concatenate((np.zeros(1), Planet_rotation_O4))
Planet_rotation_O4 = np.concatenate((Planet_rotation_O4, np.zeros(1)))
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 25)
#Approx_Longitudes = [-84.35599219825002, -79.40499345802908, -72.07684175120923, -65.36564451357754, -58.84718985194722, -52.31004511289239, -45.62372872528929, -38.70178599365309, -31.49163010200687, -23.975094175512766, -16.17296685077457, -8.147891599318474, 0.0, 8.147891599318474, 16.17296685077457, 23.975094175512766, 31.49163010200687, 38.70178599365309, 45.62372872528929, 52.31004511289239, 58.84718985194722, 65.36564451357754, 72.07684175120923, 79.40499345802908, 84.35599219825002]

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '12:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), np.flip(Velocity_QNAS[0:26]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), np.flip(Velocity_Err_QNAS_N[0:26]), np.flip(Velocity_Err_QNAS_P[0:26]), color='b', alpha=0.5)
ax2.plot(np.arange(26), np.flip(ex_INTSQNAS[0:26]), color='r', label='IR Intensity', lw = 5)
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
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.6)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

Vel_UranusNAS.append(np.flip(Velocity_QNAS[0:26]))
Vel_UranusNAS_Err.append(np.flip(Velocity_QNAS_Err[0:26]))
Int_UranusNAS.append(np.flip(ex_INTSQNAS[0:26]))
Int_UranusNAS_Err.append(np.flip(QNAS_IntErr[0:26]))

#%% Now we try to implement the velocities and Ints as an average together

Tot_Weights = [0, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 8, 0]

Vels_NAS = np.sum([Vel_UranusNAS], axis=1)/Tot_Weights
Vels_NAS_Err = np.sum([Vel_UranusNAS_Err], axis=1)/Tot_Weights
Int_NAS = np.sum([Int_UranusNAS], axis=1)/Tot_Weights
Int_NAS_Err = np.sum([Int_UranusNAS_Err], axis=1)/Tot_Weights

Vels_NAS = np.array(Vels_NAS[0,:])
Vels_NAS_Err = np.array(Vels_NAS_Err[0,:])
Int_NAS = np.array(Int_NAS[0,:])
Int_NAS_Err = np.array(Int_NAS_Err[0,:])

Int_UranusNAS_P = []
Int_UranusNAS_N = []

for n in range(26):
    if n < 2 or n == 25:
        Int_UranusNAS_P.append(Int_NAS[n])
        Int_UranusNAS_N.append(Int_NAS[n])
    else:
        Int_UranusNAS_P.append(Int_NAS[n] + Int_NAS_Err[n])
        Int_UranusNAS_N.append(Int_NAS[n] - Int_NAS_Err[n])

Vels_NAS_P = []
Vels_NAS_N = []

for n in range(26):
    if n < 2 or n == 25:
        Vels_NAS_P.append(Vels_NAS[n])
        Vels_NAS_N.append(Vels_NAS[n])
    else:
        Vels_NAS_P.append(Vels_NAS[n] + Vels_NAS_Err[n])
        Vels_NAS_N.append(Vels_NAS[n] - Vels_NAS_Err[n])
        
Planet_rotation_O4b = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 23)
Planet_rotation_O4b = np.concatenate((np.zeros(2), Planet_rotation_O4b))
Planet_rotation_O4b = np.concatenate((Planet_rotation_O4b, np.zeros(1)))

Int_NAS_Errv2 = Int_NAS_Err
Int_NAS_Errv2[Int_NAS_Err == np.inf] = np.nan

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('Average 30 minute exposure on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(np.arange(26), Planet_rotation_O4b, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(np.arange(26), Vels_NAS, color='b', ls= '--', label='Combined Q(1,0), Q(3,0) and (3,1) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Vels_NAS_N, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(np.arange(26), Vels_NAS_N, Vels_NAS_P, color='b', alpha=0.5)
ax2.plot(np.arange(26), Int_NAS, color='r', label='IR Intensity', lw = 5)
ax2.fill_between(np.arange(26), Int_UranusNAS_N, Int_UranusNAS_P, color='r', alpha=0.5)
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
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Int_NAS_Errv2), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 0.5)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
