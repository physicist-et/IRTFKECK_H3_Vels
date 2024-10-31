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
from numpy import array

A_average = -0.00010685173336385038
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
        
A_average = -0.00010685173336385038
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
from matplotlib.patches import Circle

h3p = h3ppy.h3p()
wave2 = h3p.wavegen(np.nanmin(Wavelength_O4), np.nanmax(Wavelength_O4), 2048)
Total_IRTF_Data_Alt = Total_IRTF_Data/90
Total_IRTF_Data_Alt = acre(Total_IRTF_Data_Alt, width = 15, verbose = False)
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt > 0.025] = 0
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt < -0.025] = 0

plt.figure()
plt.imshow(Total_IRTF_Data_Alt, cmap='gist_gray')
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
plt.imshow(Total_IRTF_Data_AltQ3, cmap='gist_gray')
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

#%% Now plot a single line of data and make a line graph
#Total_IRTF_Data = acre(Total_IRTF_Data[330:551,:], width = 5, verbose = True)
T_IRTFData = Total_IRTF_Data_Alt
h3p = h3ppy.h3p()
wave2 = h3p.wavegen(np.nanmin(Wavelength_O4), np.nanmax(Wavelength_O4), 2048)
import more_itertools as mit

# Slit_NumbersP = []
# Slit_NumbersN = []
# Mid_PointsP = []
# Mid_PointsN = []
# SeperationN = []
# SeperationP = []

# #Condition = wave >= 3.96448
# #print(Condition)
# #print(np.extract(Condition, wave))

# #First lets find the exact positions of the slit ends

# DATA = T_IRTFData
# Noise = DATA[550:600, 500:600]
# Noise_Mean = np.mean(Noise)
# Noise_Std = np.std(Noise)
# Noise_Max = Noise_Mean
# Noise_Min = Noise_Mean
# for i in range(12):
#     Pixel_NumbersP = []
#     Pixel_NumbersN = []
#     Start = 1208
#     Line_No = Start + i
#     DATA_Line = DATA[570:645, Line_No]
#     for ii in range(75):
#         DATA_Pixel = DATA_Line[ii]
#         if DATA_Pixel > Noise_Max:
#             Pixel_NumbersP.append(ii)
#         elif DATA_Pixel < Noise_Min:
#             Pixel_NumbersN.append(ii)
#         else:
#             pass
#     GroupsP = [list(group) for group in mit.consecutive_groups(Pixel_NumbersP)]
#     Slit_NosP = max(GroupsP, key=len)
#     Start_PointP = min(Slit_NosP)
#     Finish_PointP = max(Slit_NosP)
#     Slit_NumbersP.append(Start_PointP+570)
#     Slit_NumbersP.append(Finish_PointP+570)
#     MidsP = int((Finish_PointP - Start_PointP)/2) + Start_PointP
#     Mid_PointsP.append(MidsP+570)
#     SepsP = int(Finish_PointP - Start_PointP)
#     SeperationP.append(SepsP)
#     GroupsN = [list(group) for group in mit.consecutive_groups(Pixel_NumbersN)]
#     Slit_NosN = max(GroupsN, key=len)
#     Start_PointN = min(Slit_NosN)
#     Finish_PointN = max(Slit_NosN)
#     Slit_NumbersN.append(Start_PointN+570)
#     Slit_NumbersN.append(Finish_PointN+570)
#     MidsN = int((Finish_PointN - Start_PointN)/2) + Start_PointN
#     SepsN = int(Finish_PointN - Start_PointN)
#     Mid_PointsN.append(MidsN+570)
#     SeperationN.append(SepsN)
    
# #To work out approximate pixel length of Uranus
# Start_length_UranusP = 579
# Finish_length_UranusP = 603.375

# Y_length_UranusP = Finish_length_UranusP - Start_length_UranusP
# X_length_UranusP = 5
# Slit_lengthP = np.sqrt((Y_length_UranusP)**2 + (X_length_UranusP)**2)

# Start_length_UranusN = 620
# Finish_length_UranusN = 644.375

# Y_length_UranusN = Finish_length_UranusN - Start_length_UranusN
# X_length_UranusN =  5
# Slit_lengthN = np.sqrt((Y_length_UranusN)**2 + (X_length_UranusN)**2)
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
    Order_4_shift = np.flipud(np.rot90(Order_4_shift))
    Combined_Q1.append(Order_4_shift)
        
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
    Order_5_shift = np.flipud(np.rot90(Order_5_shift))
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
    if x == 8:
        for m in range(25):
            datai = dataI[15+m, :]
            gmodel = Model(gauss_fit)
            a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
            a1_pointP = a1_pointP[0][0] - 1110
            ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
            ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
            ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
            A0 = ex_pQ1P.a0
            ex_A0_Q18.append(A0)
            A2 = ex_pQ1P.a2
            ex_A2_Q18.append((A2)*1.6199707031249897e-05)
            err_A0 = ex_eQ1P[0]
            ex_Err_A0_Q18.append(err_A0)
            err_A2 = ex_eQ1P[2]
            ex_Err_A2_Q18.append(err_A2*1.6199707031249897e-05)
            ex_A1_Q18.append(ex_pQ1P.a1)
            ex_Err_A1_Q18.append(ex_eQ1P[1]*1.6199707031249897e-05)
            ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
            ex_FWHMQ18.append(ex_FWHM)
            ex_INTSQ18.append(A0*ex_FWHM*(10**6))           
    if x == 3:
        for m in range(22):
            datai = dataI[16+m, :]
            gmodel = Model(gauss_fit)
            a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
            a1_pointP = a1_pointP[0][0] - 1110
            ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
            ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
            ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
            A0 = ex_pQ1P.a0
            ex_A0_Q18.append(A0)
            A2 = ex_pQ1P.a2
            ex_A2_Q18.append((A2)*1.6199707031249897e-05)
            err_A0 = ex_eQ1P[0]
            ex_Err_A0_Q18.append(err_A0)
            err_A2 = ex_eQ1P[2]
            ex_Err_A2_Q18.append(err_A2*1.6199707031249897e-05)
            ex_A1_Q18.append(ex_pQ1P.a1)
            ex_Err_A1_Q18.append(ex_eQ1P[1]*1.6199707031249897e-05)
            ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
            ex_FWHMQ18.append(ex_FWHM)
            ex_INTSQ18.append(A0*ex_FWHM*(10**6))
            # plt.figure()
            # plt.plot(XX, datai[1110:1310], 'o')
            # plt.plot(XX, ex_resultQ1P.best_fit, '-')
    for m in range(22):
        datai = dataI[16+m, :]
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
        # if x == 9:
        #     plt.figure()
        #     plt.plot(XX, datai[1110:1310], 'bo')   
        #     plt.plot(XX, ex_resultQ1P.best_fit, '--')
        #     print(ex_pQ1P.a1)
      
Q1_IntErr = []
#Slit is 15 so

for o in range(242):
    Q1IntErr = ex_INTSQ1[o]*np.sqrt((ex_Err_A0_Q1[o]/ex_A0_Q1[o])**2 + (ex_Err_A2_Q1[o]/ex_A2_Q1[o])**2)
    Q1_IntErr.append(Q1IntErr)

o = 0
for x in range(11):
    dataI = (AB_Combined_Order5[o]+AB_Combined_Order5[o+1]+AB_Combined_Order5[o+2]+AB_Combined_Order5[o+3])/8
    dataI[dataI < -0.004] = 0
    Q3_Data_Capture.append(dataI[0:50, :])
    o += 4
    for m in range(22):
        datai = dataI[16+m, :]
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
        err_A0 = ex_eQ1P[0]
        ex_Err_A0_Q3.append(err_A0)
        err_A2 = ex_eQ1P[2]
        ex_Err_A2_Q3.append(err_A2*1.6199707031249897e-05)
        ex_A1_Q3.append(ex_pQ3P.a1)
        ex_Err_A1_Q3.append(ex_eQ3P[1]*1.6199707031249897e-05)
        ex_FWHM = A2*1.6199707031249897e-05*2*np.sqrt(2*math.log(2))
        ex_FWHMQ3.append(ex_FWHM)
        ex_INTSQ3.append(A0*ex_FWHM*(10**6)) #Needs work with Q3 to skip difficult pixels?
        # if x == 3 or x == 2:
        #     if m > 5 and m < 15:
        #         plt.figure()
        #         plt.plot(XX, datai[1260:1460], 'bo')
        #         z = (XX - ex_pQ3P.a1) / ex_pQ3P.a2
        #         y = ex_pQ3P.a0 * np.exp(-z**2 / ex_pQ3P.a2) + ex_pQ3P.a3 + ex_pQ3P.a4 * XX + ex_pQ3P.a5 * XX**2
        #         plt.plot(XX, y, 'g')
        #         print(A0*ex_FWHM*(10**6))       
      
Q3_IntErr = []
#Slit is 15 so

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
Mapping_Q1B[Mapping_Q1B > 2.0] = 0
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
Q1Int = np.fliplr(Mapping_Q1B)
Q3Int = np.fliplr(Mapping_Q3B)
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
Temps = np.array(Temps)
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
plt.yticks(np.arange(10), TIMEv2[1:11])

#Temperature Errors
Temperature_Err = (np.nanmean(Mapping_Q1A)/np.nanmean(Mapping_Q3B))*np.sqrt((0.06/np.nanmean(Mapping_Q1A))**2+(0.05/np.nanmean(Mapping_Q3B))**2)
Temps_Err = Temperature_Err/(np.nanmean(Mapping_Q1A)/np.nanmean(Mapping_Q3B))
Temp_Err = np.nanmean(Temps)*(Temps_Err/math.log(np.nanmean(Mapping_Q1A)/np.nanmean(Mapping_Q3B)))

print(Temp_Err)

#%% Look to convert the temperature to CDs
# a0 = [-1.11391, 0.0581076, 0.000302967, -2.83724*10**-7, 2.31119*10**-10, -7.15895*10**-14, 1.00150*10**-17]
      
# Cols_Den = []

# for i in range(11): #etc
#     for m in range(22):
#         T = Temps[i]
#         Totalz = 0
#         for a in range(7):
#             Totalz += a0[a]*(np.log10(T[m])**a)
#         Imodel_Q1 = ((gQ1 * JQ1  *wQ1 * AQ1 * h * c * 100)/(4*np.pi*Totalz))*np.exp((h*c*(-2552.57)*100)/(kb*T[m]))
#         Cols_Den.append(Mapping_Q1A[i,m]/(Imodel_Q1*10**6))

#%% Now to find the A1 value and plot that

A1_Q1 = []
A1_Q3 = []
A1_Q18 = []
A1_Q38 = []
A1_Q12 = []

for n in range(242):
    A1 = ((ex_A1_Q1[n]-np.nanmean(ex_A1_Q1))*-1.6199707031249897e-05)
    A1_Q1.append(A1)

for n in range(242):
    A1 = ((ex_A1_Q3[n]-np.nanmean(ex_A1_Q3))*-1.6199707031249897e-05)
    A1_Q3.append(A1)

# for n in range(25):
#     A1 = ((ex_A1_Q18[n+22]-np.nanmean(ex_A1_Q18[25:44]))*1.5709811506855093e-05)
#     A1_Q18.append(A1)

for n in range(22):
    A1 = ((ex_A1_Q18[n]-np.nanmean(ex_A1_Q18[3:19]))*1.6199707031249897e-05)
    A1_Q12.append(A1)
    
# for n in range(25):
#     A1 = ((ex_A1_Q38[n]-np.nanmean(ex_A1_Q38))*1.6361537623956223e-05)
#     A1_Q13.append(A1)    

Q1_line_shift = np.reshape(A1_Q1, (11, 22))
Q3_line_shift = np.reshape(A1_Q3, (11, 22))

#Now lets calculate the planetary rotation
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
Limb_velocity_pix_O5 = (Period_Uranus/1000) #Need to see if there's slowing or speeding up of atmospheric/ionospheric winds at Uranus

c = 299792458 
lamda_H3_O4 = 3.95295 
lamda_H3_O5 = 3.98558
Velocity_O4 = []
Velocity_O5 = []
Velocity_O48 = []
Velocity_O58 = []
Velocity_O42 = []

for n in range(242):
    V = ((A1_Q1[n]/lamda_H3_O4)*c)
    Velocity_O4.append(V/1000)

for n in range(242):
    V = ((A1_Q3[n]/lamda_H3_O5)*c)
    Velocity_O5.append(V/1000)

# for n in range(25):
#     V = ((A1_Q18[n]/lamda_H3_O4)*c)
#     Velocity_O48.append(V/1000)

for n in range(22):
    V = ((A1_Q12[n]/lamda_H3_O4)*c)
    Velocity_O42.append(V/1000)   
    
# for n in range(25):
#     V = ((A1_Q38[n]/lamda_H3_O5)*c)
#     Velocity_O58.append(V/1000)
    
# for n in range(242):
#     V = ((A1_Q1[n]+lamda_H3_O4)*c)/(75000)
#     Velocity_O4.append(V/1000)

# for n in range(242):
#     V = ((A1_Q3[n]+lamda_H3_O5)*c)/(75000)
#     Velocity_O5.append(V/1000)
    
Velocity_O4_Err = []
Velocity_O48_Err = []
Velocity_O5_Err = []
Velocity_O58_Err = []
Velocity_O42_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(242):
    Err = Velocity_O4[n]*np.sqrt((ex_Err_A1_Q1[n]/A1_Q1[n])**2)
    Velocity_O4_Err.append(np.sqrt(Err**2))
    
# for n in range(25):
#     Err = Velocity_O48[n]*np.sqrt((ex_Err_A1_Q18[n+22]/A1_Q18[n])**2)
#     Velocity_O48_Err.append(np.sqrt(Err**2))

for n in range(242):
    Err = Velocity_O5[n]*np.sqrt((ex_Err_A1_Q3[n]/A1_Q3[n])**2)
    Velocity_O5_Err.append(np.sqrt(Err**2))
    
for n in range(22):
    Err = Velocity_O42[n]*np.sqrt((ex_Err_A1_Q18[n]/A1_Q12[n])**2)
    Velocity_O42_Err.append(np.sqrt(Err**2)) 

# for n in range(25):
#     Err = Velocity_O58[n]*np.sqrt((ex_Err_A1_Q38[n]/A1_Q38[n])**2)
#     Velocity_O58_Err.append(np.sqrt(Err**2))

Velocities_O4 = np.reshape(Velocity_O4, (11, 22))
Velocities_O4 = np.fliplr(np.flipud(Velocities_O4))

Velocity_Err_O4 = np.reshape(Velocity_O4_Err, (11, 22))
Velocity_Err_O4 = np.fliplr(np.flipud(Velocity_Err_O4))

Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 22)
Planet_rotation_O5 = np.linspace(Limb_velocity_pix_O5, -1*Limb_velocity_pix_O5, 22)
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 22)
Planet_rotation_O4_80 = np.linspace(Limb_velocity_pix_O4_80, -1*Limb_velocity_pix_O4_80, 22)

# Velocity_O48 = np.flip(Velocity_O48)
# Velocity_O42= np.flip(Velocity_O42)
# Velocity_O42_Err = np.flip(Velocity_O42_Err)
# Velocity_O48_Err= np.flip(Velocity_O48_Err)

TIME = ['07:45', '08:14', '08:44', '09:14', '09:44', '10:13', '10:43', '11:12', '11:42', '12:11', '12:42']

for a in range(11):
    k = a
    fig, ax = plt.subplots(figsize=(10,6))
    ax2 = ax.twinx()
    ax.set_title('Approximate ion flows from $H_{3}^{+}$ $Q(1,0^{-})$ over a ~30 minute exposure at ' +str(TIME[k]) +' UT', fontsize=30)
    ax.plot(np.arange(22), Planet_rotation_O4, color='k', label='Planetary Rotation')
    ax.plot(np.arange(22), Planet_rotation_O4_90, color='g', label='Reduced Planetary Rotation by ~(257m/s)')
    ax.plot(np.arange(22), Velocities_O4[k,:], color='b', label='IR Velocities', lw = 3)
    ax.fill_between(np.arange(22), Velocities_O4[k,:]-Velocity_Err_O4[k,:], Velocities_O4[k,:]+Velocity_Err_O4[k,:], color='b', alpha=0.5)
    ax2.plot(np.arange(22), np.flip(Mapping_Q1A[k,:]), color='r', label='IR Intensity', lw = 3)
    #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
    ax.set_xticks(np.linspace(0, 21, 7))
    ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'])
    # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
    ax.set_xlabel('Rough Arbitrary Uranian Longitude ($^\circ$)', fontsize=25)
    ax.set_ylabel('$H_{3}^{+}$ $Q(1,0^{-})$ ion flows across Uranus ($kms^{-1}$)', fontsize=25)
    ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
    ax.vlines(np.linspace(0, 21, 7), -200, 200, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    ax.hlines(np.linspace(-150, 150, 7), 0, 21, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    ax.tick_params(axis='both', labelsize=20)
    ax2.tick_params(axis='both', labelsize=20)
    ax.set_xlim(0, 21) #0, 3601
    ax.set_ylim(-3, 3)
    ax2.set_ylim(0, 1.0)
    ax.legend(loc = 'lower left', fontsize=30)
    ax2.legend(loc = 'lower right', fontsize=30)
#plt.title('Approximate ion velocities from $H_{3}^{+}$ $Q(1,0^{-})$ over a 30 to 40 minute exposure at ' +str(TIME[k]) +' UT', fontsize=30)
#%% Now to make an average profile of emission and velocities on the 15th #To do (Redo Lats and Longs to right values)

Velocity_Av = []
Velocity_Err = []
Ints_Av = []
Ints_Err = []

for a in range(22):
    Av_Vel = (Velocities_O4[0,a] + Velocities_O4[1,a] + Velocities_O4[4,a] + Velocities_O4[5,a] + Velocities_O4[8,a] + Velocities_O4[9,a])/6
    Velocity_Av.append(Av_Vel)
    Av_Vel_Err = np.sqrt(((Velocities_O4[0,a]-Av_Vel)**2+(Velocities_O4[1,a]-Av_Vel)**2+(Velocities_O4[4,a]-Av_Vel)**2+(Velocities_O4[5,a]-Av_Vel)**2+(Velocities_O4[8,a]-Av_Vel)**2+(Velocities_O4[9,a]-Av_Vel)**2)/5)
    Velocity_Err.append(Av_Vel_Err) 
    Av_Ints = (Mapping_Q1A[0,21-a] + Mapping_Q1A[1,21-a] + Mapping_Q1A[4,21-a] + Mapping_Q1A[5,21-a] + Mapping_Q1A[8,21-a] + Mapping_Q1A[9,21-a])/6
    Ints_Av.append(Av_Ints)
    Err_Ints = np.sqrt(((Mapping_Q1A[0,21-a]-Av_Ints)**2 + (Mapping_Q1A[1,21-a]-Av_Ints)**2 + (Mapping_Q1A[4,21-a]-Av_Ints)**2 + (Mapping_Q1A[5,21-a]-Av_Ints)**2 + (Mapping_Q1A[8,21-a]-Av_Ints)**2 + (Mapping_Q1A[9,21-a]-Av_Ints)**2)/5)
    Ints_Err.append(Err_Ints)
    
print(Planet_rotation_O4[0]-Planet_rotation_O4[-1] - (Velocity_Av[0]-Velocity_Av[-1]))

#%%Err Prop

TotErr = np.sqrt(2*(np.nanmean(Av_Vel_Err)**2))

#%%
PosErr = []
NegErr = []
PosErr2 = []
NegErr2 = []
Av_PosErr = []
Av_NegErr= []

#for n in range(25):
#     PosErr.append((Velocity_O48[n]+Velocity_O48_Err[n]))
#     NegErr.append((Velocity_O48[n]-Velocity_O48_Err[n]))
# #     PosErr.append(Velocities_O4[2,n]+Velocity_O4_Err[2,n])
# #     NegErr.append(Velocity_O48[n]-Velocity_O48_Err[n])

for n in range(22):
    Av_PosErr.append(Velocity_Av[n]+Velocity_Err[n])
    Av_NegErr.append(Velocity_Av[n]-Velocity_Err[n])

LongitudesS = np.load('LongitudesN.npy')
Mean_Longs = []

# for a in range(21):
#     if a == 0 or a == 20:
#         if a == 0:
#             Mean_Longs.append(int(-95))
#         else: 
#             Mean_Longs.append(int(95))
#     elif a == 1 or a == 19:
#         Mean_Longitude = (LongitudesS[a-1] + -1*LongitudesS[a+17])/2
#     else:
#         Mean_Longitude = (LongitudesS[a-1] + LongitudesS[a+17])/2
#     Mean_Longs.append(int(Mean_Longitude))

for a in range(22):
    Mean_Longitude = LongitudesS[a+22]
    if Mean_Longitude > 70:
        Mean_Longs.append(int(90))
    elif Mean_Longitude < -70:
        Mean_Longs.append(int(-90))
    else:
        Mean_Longs.append(int(Mean_Longitude))

# for n in range(22):
#     PosErr2.append(-1*(Velocity_O42[n]-Velocity_O42_Err[n]))
#     NegErr2.append(-1*(Velocity_O42[n]+Velocity_O42_Err[n])) #Return to figures and fix these
    
for n in range(22):
    PosErr2.append(Velocity_Av[n] + Velocity_Err[n])
    NegErr2.append(Velocity_Av[n] - Velocity_Err[n])

Vels_Set3 = Velocities_O4[3,:]
# IntsSet3 = ex_INTSQ18[22:47]
Errs_Set3 = Velocity_Err_O4[3,:]
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('Approximate ion velocities from $H_{3}^{+}$ $Q(1,0^{-})$ over a ~30 minute exposure at ' + '09:14' +' UT', fontsize=25)
ax.plot(Mean_Longs, np.zeros(22), color='k', ls = '--', label='Planetary Rotation')
ax.plot(Mean_Longs, Vels_Set3[0:22] - Planet_rotation_O4, color='b', label='Q(1,0) IR Velocities', lw = 5)
ax.plot(Mean_Longs, (Velocity_Av - Planet_rotation_O4), color='c', label='Mean IR Velocities', lw = 5, alpha = 0.5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs, (Vels_Set3[0:22] - Errs_Set3[0:22]) - Planet_rotation_O4, (Vels_Set3[0:22] + Errs_Set3[0:22]) - Planet_rotation_O4, color='b', alpha=0.5)
ax.fill_between(Mean_Longs, Av_NegErr - Planet_rotation_O4, Av_PosErr - Planet_rotation_O4, color='c', alpha = 0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs, np.flip(Mapping_Q1A[3,:]), color='r', label='IR Intensity', lw = 5)
ax2.plot(Mean_Longs, Ints_Av, color='m', label='Mean IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Ints_Err), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-1.5, 1.5, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-1.5, 1.5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='upper left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
#ax.annotate('Approximate Latitude = -23' + '$^\circ$' + u'\u00B1' + ' 9' +'$^\circ$', (12, 2.1), fontsize=25)

#Observer RF
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('Observations from 07:30 to 13:00 UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Mean_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
ax.plot(Mean_Longs, Planet_rotation_O4_90, color='darkgreen', ls = '--', label='Reduced Planetary Rotation by ~(257m/s)')
ax.plot(Mean_Longs, Planet_rotation_O4_80, color='lightgreen', ls = '--', label='Reduced Planetary Rotation by ~(514m/s)')
#ax.plot(Mean_Longs, Vels_Set3[0:22], color='b', label='Q(1,0) IR Velocities', lw = 5, ls='--')
ax.plot(Mean_Longs, Velocity_Av, color='c', label='Mean IR Velocities', lw = 5, alpha = 0.5, ls='--')
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(Mean_Longs, (Vels_Set3[0:22] - Errs_Set3[0:22]), (Vels_Set3[0:22] + Errs_Set3[0:22]), color='b', alpha=0.5)
ax.fill_between(Mean_Longs, Av_NegErr, Av_PosErr, color='c', alpha = 0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
#ax2.plot(Mean_Longs, np.flip(Mapping_Q1A[3,:]), color='r', label='IR Intensity', lw = 5)
ax2.plot(Mean_Longs, Ints_Av, color='m', label='Mean IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), fontsize=15)
ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=20)
ax.set_ylabel('LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=20)
ax2.set_ylabel(r'Intensity ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=20)
ax2.legend(loc='upper right', fontsize=20)
#ax.annotate('Approximate Latitude = -23' + '$^\circ$' + u'\u00B1' + ' 9' +'$^\circ$', (12, 2.1), fontsize=25)

#%% Lets tweak the velocities by the amount its off by to see what occur




#%%
Vels_Set8 = Velocities_O4[8,:]
IntsSet3 = ex_INTSQ18[22:47]
Errs_Set8 = Velocity_Err_O4[8,:]

Mean_Longs5 = []

for a in range(22):
    Mean_Longitude = LongitudesS[a+22]
    if Mean_Longitude > 70:
        Mean_Longs5.append(int(310))
    elif Mean_Longitude < -70:
        Mean_Longs5.append(int(140))
    else:
        Mean_Longs5.append(int(Mean_Longitude)+224)
        
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Mean_Longs5, np.zeros(22), color='k', ls = '--', label='Planetary Rotation')
ax.plot(Mean_Longs5, Vels_Set8[0:22] - Planet_rotation_O4, color='b', label='Q(1,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs5, (Vels_Set8[0:22] - Errs_Set8[0:22]) - Planet_rotation_O4, (Vels_Set8[0:22] + Errs_Set8[0:22]) - Planet_rotation_O4, color='b', alpha=0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs5, np.flip(Mapping_Q1A[8, :]), color='r', label='IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(150, 300, 6), fontsize=20)
ax.set_xticklabels(labels=['150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(150, 300, 6), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-1.5, 1.5, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(140, 300) #0, 3601
ax.set_ylim(-1.5, 1.5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='upper left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad = 25)
ax.plot(Mean_Longs5, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
ax.plot(Mean_Longs5, Planet_rotation_O4_90, color='darkgreen', ls = '--', label='Reduced Planetary Rotation by ~(257m/s)')
ax.plot(Mean_Longs5, Planet_rotation_O4_80, color='lightgreen', ls = '--', label='Reduced Planetary Rotation by ~(514m/s)')
ax.plot(Mean_Longs5, Vels_Set8[0:22], color='b', ls='--', label='Q(1,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs5, (Vels_Set8[0:22] - Errs_Set8[0:22]), (Vels_Set8[0:22] + Errs_Set8[0:22]), color='b', alpha=0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs5, np.flip(Mapping_Q1A[8, 0:22]), color='r', label='IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(150, 300, 6), fontsize=15)
ax.set_xticklabels(labels=['150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Horizon Uranian Longitude ($^\circ$)', fontsize=20)
ax.set_ylabel('LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=20)
ax2.set_ylabel(r'Intensity ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
ax.vlines(np.linspace(150, 300, 6), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-3, 3, 7), 140, 300, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(140, 300) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=20)
ax2.legend(loc='lower right', fontsize=20)
#ax.annotate('Approximate Latitude = -56' + '$^\circ$' + u'\u00B1' + ' 11' +'$^\circ$', (12, 2.1), fontsize=25)

#%% Add this event as it also has aurora
Vels_Set8 = Velocities_O4[9,:]
IntsSet3 = ex_INTSQ18[22:47]
Errs_Set8 = Velocity_Err_O4[9,:]
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Mean_Longs, np.zeros(22), color='k', ls = '--', label='Planetary Rotation')
ax.plot(Mean_Longs, Vels_Set8[0:22] - Planet_rotation_O4, color='b', label='Q(1,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs, (Vels_Set8[0:22] - Errs_Set8[0:22]) - Planet_rotation_O4, (Vels_Set8[0:22] + Errs_Set8[0:22]) - Planet_rotation_O4, color='b', alpha=0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs, np.flip(Mapping_Q1A[9, :]), color='r', label='IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), fontsize=20)
ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-1.5, 1.5, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-1.5, 1.5)
ax2.set_ylim(0, 1.0)
ax.legend(loc='upper left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~30 minute exposure at ' + '11:42' +' UTC on $11^{th}$ October 2016', fontsize=22, pad = 25)
ax.plot(Mean_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Mean_Longs, Planet_rotation_O4_90, color='darkgreen', ls = '--', label='Reduced Planetary Rotation by ~(257m/s)')
#ax.plot(Mean_Longs, Planet_rotation_O4_80, color='lightgreen', ls = '--', label='Reduced Planetary Rotation by ~(514m/s)')
ax.plot(Mean_Longs, Vels_Set8[0:22], color='b', ls='--', label='Q(1,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs, (Vels_Set8[0:22] - Errs_Set8[0:22]), (Vels_Set8[0:22] + Errs_Set8[0:22]), color='b', alpha=0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs, np.flip(Mapping_Q1A[9, 0:22]), color='r', label='IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), fontsize=15)
ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=20)
ax.set_ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=20)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=20)
ax2.legend(loc='lower right', fontsize=20)

#%%
# Velocities_O5 = np.reshape(Velocity_O5, (11, 22))
# Velocities_O5 = Velocities_O5.transpose()

# Velocity_Err_O5 = np.reshape(Velocity_O5_Err, (11, 22))
# Velocity_Err_O5 = Velocity_Err_O5.transpose()

# # for a in range(11):
# #     fig, ax = plt.subplots(figsize=(12,8))
# #     ax2 = ax.twinx()
# #     ax.set_title('Approximate ion velocities from $H_{3}^{+}$ $Q(1,0^{-})$ over a ~60 minute exposure at ' +str(TIME[a]) +' UT', fontsize=25)
# #     ax.plot(np.arange(22), np.flip(Planet_rotation_O5), color='k', label='Planetary Rotation')
# #     ax.plot(np.arange(22), np.flip(Velocities_O5[:,a]), color='b', label='IR Velocities', lw = 3)
# #     ax.fill_between(np.arange(22), np.flip(Velocities_O5[:,a])-np.flip(Velocity_Err_O5[:,a]), np.flip(Velocities_O5[:,a])+np.flip(Velocity_Err_O5[:,a]), color='b', alpha=0.5)
# #     ax2.plot(np.arange(22), np.flip(Mapping_Q3B[10-a, :]), color='r', label='IR Intensity', lw = 3)
# #     # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #     # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #     # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# #     ax.set_xticks(np.linspace(0, 21, 7))
# #     ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'])
# #     # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# #     ax.set_xlabel('Rough Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# #     ax.set_ylabel('Wind Speeds from $H_{3}^{+}$ $Q(3,0^{-})$ across Uranus (ULS) ($^\circ$)', fontsize=25, labelpad=15)
# #     ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) +- ' + str("{:.2f}".format(round(np.nanmean(Q3_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
# #     ax.vlines(np.linspace(0, 21, 7), -200, 200, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #     ax.hlines(np.linspace(-150, 150, 7), 0, 21, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #     ax.tick_params(axis='both', labelsize=15)
# #     ax2.tick_params(axis='both', labelsize=15)
# #     ax.set_xlim(0, 21) #0, 3601
# #     ax.set_ylim(-200, 200)
# #     ax2.set_ylim(0, 0.75)
# #     ax.legend(fontsize=15)
# #     ax2.legend(loc='upper left', fontsize=15)
# #     #plt.title('Approximate ion velocities from $H_{3}^{+}$ $Q(3,0^{-})$ over a 60 to 70 minute exposure at ' + str(TIME[k]) + ' UT', fontsize=30)

# #%%Add all lines #NEED TO CHECK FROM HERE ONWARDS!!!!!! Check the rotation rate and why rather 200m than 20000m 
# #For Q1

# Test_1 = (IRTF_Data_Q1[:,:,0] + IRTF_Data_Q1[:,:,1] + IRTF_Data_Q1[:,:,2] + IRTF_Data_Q1[:,:,3])/8

# plt.figure()
# plt.imshow(Test_1, cmap='gist_gray', vmax = np.nanmean(Test_1)+2*np.nanstd(Test_1), vmin = np.nanmean(Test_1)-2*np.nanstd(Test_1))

# Shape_Order4 = np.shape(Test_1)

# for n in range(Shape_Order4[1]):
#     Spec_col = Test_1[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
#     shifted_Spec_col = f.shift(Spec_col, shift=[-41.2337], mode='wrap') #41.6508 for Order 5
#     if n == 0:
#         Order_4_shift = shifted_Spec_col
#     else:
#         Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))

# Test_Total = np.flipud(np.rot90(Order_4_shift)) - Test_1

# Circle1 = plt.Circle((1214, 27), radius=12.049815498154981, color ='r', lw=2, ls='--', fill=False)
# fig, ax = plt.subplots()
# ax.imshow(Test_Total, cmap='gist_gray', vmax = np.nanmean(Test_1)+np.nanstd(Test_1), vmin = np.nanmean(Test_1)-np.nanstd(Test_1))
# ax.add_patch(Circle1)

# #For Q3

# Test_1 = (IRTF_Data_Q3[:,:,0] + IRTF_Data_Q3[:,:,1] + IRTF_Data_Q3[:,:,2] + IRTF_Data_Q3[:,:,3])/8

# plt.figure()
# plt.imshow(Test_1, cmap='gist_gray', vmax = np.nanmean(Test_1)+2*np.nanstd(Test_1), vmin = np.nanmean(Test_1)-2*np.nanstd(Test_1))

# Shape_Order5 = np.shape(Test_1)

# for n in range(Shape_Order5[1]):
#     Spec_col = Test_1[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
#     shifted_Spec_col = f.shift(Spec_col, shift=[-41.6508], mode='wrap') #41.6508 for Order 5
#     if n == 0:
#         Order_5_shift = shifted_Spec_col
#     else:
#         Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))

# Test_Total = np.flipud(np.rot90(Order_5_shift)) - Test_1

# Circle2 = plt.Circle((1358,26), radius=12.049815498154981, color ='r', lw=2, ls='--', fill=False)

# fig, ax = plt.subplots()
# ax.imshow(Test_Total, cmap='gist_gray', vmax = np.nanmean(Test_1)+np.nanstd(Test_1), vmin = np.nanmean(Test_1)-np.nanstd(Test_1))
# ax.add_patch(Circle2)

# #So lets complete this with all the data and see what we can get out in a loop

# Shape_Order4 = np.shape(IRTF_Data_Q1[:,:,0])
# Shape_Order5 = np.shape(IRTF_Data_Q3[:,:,0])

# k = 0

# ABBABins_Order4 = []
# ABBABins_Order5 = []

# for i in range(22):
#     ABBA_Bin_Q1 = (IRTF_Data_Q1[:,:,k] + IRTF_Data_Q1[:,:,k+1])/4
#     ABBA_Bin_Q3 = (IRTF_Data_Q3[:,:,k] + IRTF_Data_Q3[:,:,k+1])/4
#     for n in range(Shape_Order4[1]):
#         Spec_col = ABBA_Bin_Q1[:,n]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-41.2337], mode='wrap')
#         if n == 0:
#             Order_4_shift = shifted_Spec_col
#         else:
#             Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))
#     for n in range(Shape_Order5[1]):
#         Spec_col = ABBA_Bin_Q3[:,n]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-41.6508], mode='wrap')
#         if n == 0:
#             Order_5_shift = shifted_Spec_col
#         else:
#             Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))        
#     ABBA_Bin_Q1A = np.flipud(np.rot90(Order_4_shift[:,10:50])) - ABBA_Bin_Q1[10:50,:]
#     ABBA_Bin_Q3A = np.flipud(np.rot90(Order_5_shift[:,10:50])) - ABBA_Bin_Q3[10:50,:]
#     ABBABins_Order4.append(ABBA_Bin_Q1A)
#     ABBABins_Order5.append(ABBA_Bin_Q3A)
#     k += 2

# #%%Now we need to combine lines to find velocity!

# # for i in range(11):
# #     plt.figure()
# #     plt.imshow(ABBABins_Order4[i], cmap='gist_gray', vmax = np.nanmean(ABBABins_Order4[i]+np.std(ABBABins_Order4[i])), vmin = np.nanmean(ABBABins_Order4[i]-np.std(ABBABins_Order4[i])))
# #     plt.figure()
# #     plt.imshow(ABBABins_Order5[i], cmap='gist_gray', vmax = np.nanmean(ABBABins_Order5[i]+np.std(ABBABins_Order5[i])), vmin = np.nanmean(ABBABins_Order5[i]-np.std(ABBABins_Order5[i])))

# Total_Order4 = []
# Total_Order5 = []
# Total_Q1 = np.zeros((2048))
# Total_Q3 = np.zeros((2048))
# for i in range(22):
#     Q1_Strip = ABBABins_Order4[i]
#     Q3_Strip = ABBABins_Order5[i]
#     for n in range(35):
#         Total_Q1 = Q1_Strip[n,:] + Total_Q1
#         Total_Q3 = Q3_Strip[n,:] + Total_Q3
    
# Total_Order4 = Total_Q1
# Total_Order5 = Total_Q3

# # plt.figure()
# # plt.plot(np.arange(2048), Total_Order4, color = 'r')

# # plt.figure()
# # plt.plot(np.arange(2048), Total_Order5, color = 'b')

# #Now we find the positions in the 5th Order to be able to put them on top of one another

# Line_list = [468, 1359, 1455, 1955]
# A1_loc = []
# Err_A1_loc = []

# for i in range(4):
#     dataA1loc = np.where(Total_Order5 == np.max(Total_Order5[Line_list[i]-50:Line_list[i]+50]))
#     datailoc = dataA1loc[0][0] - Line_list[i] + 50
#     gmodel = Model(gauss_fit)
#     ex_resultQ3P = gmodel.fit(Total_Order5[Line_list[i]-50:Line_list[i]+50], x=np.arange(100), a0=np.nanmax(Total_Order5[Line_list[i]-50:Line_list[i]+50]), a1=datailoc, a2=1.75, a3=0, a4=0, a5=0)
#     p_Order5 = SimpleNamespace(**ex_resultQ3P.best_values)
#     ex_Order5 = np.sqrt(np.diag(ex_resultQ3P.covar))
#     A1_loc.append(p_Order5.a1 - 50 + Line_list[i])
#     Err_A1_loc.append(ex_Order5[1])

# #and for Q1
# Line_list = [1214]
# A1_loc_Q1 = []
# Err_A1_loc_Q1 = []

# for i in range(1):
#     dataA1loc = np.where(Total_Order4 == np.max(Total_Order4[Line_list[i]-50:Line_list[i]+50]))
#     datailoc = dataA1loc[0][0] - Line_list[i] + 50
#     gmodel = Model(gauss_fit)
#     ex_resultQ1P = gmodel.fit(Total_Order4[Line_list[i]-50:Line_list[i]+50], x=np.arange(100), a0=np.nanmax(Total_Order4[Line_list[i]-50:Line_list[i]+50]), a1=datailoc, a2=1.75, a3=0, a4=0, a5=0)
#     p_Order4 = SimpleNamespace(**ex_resultQ1P.best_values)
#     ex_Order4 = np.sqrt(np.diag(ex_resultQ1P.covar))
#     A1_loc_Q1.append(p_Order4.a1 - 50 + Line_list[i])
#     Err_A1_loc_Q1.append(ex_Order4[1])

# Total_speed = []
# #So to shift all onto the same bit if we set up 5 arrays shifted so they all appear at point 50 on the x and then overlay them
# for i in range(22):
#     Strip_Q1 = ABBABins_Order4[i]
#     Strip_Q3 = ABBABins_Order5[i]
#     for n in range(40):
#         Spec_col = Strip_Q1[n,:]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-A1_loc_Q1[0]+50], mode='wrap')
#         if n == 0:
#             Order_4_shift = shifted_Spec_col
#         else:
#             Order_4_shift = np.vstack((Order_4_shift, shifted_Spec_col))
#     ABBA_Bin_Q1A = Order_4_shift
#     for n in range(40):
#         Spec_col = Strip_Q3[n,:]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-A1_loc[0]+50], mode='wrap')
#         if n == 0:
#             Order_5_shift = shifted_Spec_col
#         else:
#             Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))                           
#     ABBA_Bin_Q3A = Order_5_shift
#     for n in range(40):
#         Spec_col = Strip_Q3[n,:]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-A1_loc[1]+50], mode='wrap')
#         if n == 0:
#             Order_5_shift = shifted_Spec_col
#         else:
#             Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))                           
#     ABBA_Bin_Q3B = Order_5_shift
#     for n in range(40):
#         Spec_col = Strip_Q3[n,:]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-A1_loc[2]+50], mode='wrap')
#         if n == 0:
#             Order_5_shift = shifted_Spec_col
#         else:
#             Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))                           
#     ABBA_Bin_Q3C = Order_5_shift
#     for n in range(40):
#         Spec_col = Strip_Q3[n,:]
#         shifted_Spec_col = f.shift(Spec_col, shift=[-A1_loc[3]+50], mode='wrap')
#         if n == 0:
#             Order_5_shift = shifted_Spec_col
#         else:
#             Order_5_shift = np.vstack((Order_5_shift, shifted_Spec_col))                           
#     ABBA_Bin_Q3D = Order_5_shift
#     #T_speed = (ABBA_Bin_Q1A + ABBA_Bin_Q3B)/2
#     T_speed = (ABBA_Bin_Q1A + ABBA_Bin_Q3A + ABBA_Bin_Q3B + ABBA_Bin_Q3C + ABBA_Bin_Q3D)/5
#     Total_speed.append(T_speed[:, 0:100])
    
# #%% Now to calculate the velocities approximately planet is 6 - 28 pixels so
# A1_loc_T = []
# Err_A1_loc_T = []
# Pix_A1_diff = []
# for i in range(22):
#     Data = Total_speed[i]
#     for n in range(22):
#         if n == 0:
#             Datai = Data[n+6,:]
#             dataA1loc = np.where(Datai == np.max(Datai[40:60]))
#             datailoc = dataA1loc[0][0]
#             gmodel = Model(gauss_fit)
#             result = gmodel.fit(Datai, x=np.arange(100), a0=np.nanmax(Datai[40:60]), a1=datailoc, a2=1.75, a3=0, a4=0, a5=0)
#             p_Order45 = SimpleNamespace(**result.best_values)
#             ex_Order45 = np.sqrt(np.diag(result.covar))
#             A1_loc_T.append(p_Order45.a1)
#             Err_A1_loc_T.append(ex_Order45[1])             
#         elif i == 8 and n == 4 or n == 5:
#             A1_loc_T.append(p_Order45.a1)
#             Err_A1_loc_T.append(ex_Order45[1])
#         elif i == 15 and n == 18 or n == 19:
#             A1_loc_T.append(p_Order45.a1)
#             Err_A1_loc_T.append(ex_Order45[1])            
#         else:
#             Datai = Data[n+6,:]
#             dataA1loc = np.where(Datai == np.max(Datai[40:60]))
#             datailoc = dataA1loc[0][0]
#             gmodel = Model(gauss_fit)
#             result = gmodel.fit(Datai, x=np.arange(100), a0=np.nanmax(Datai[40:60]), a1=datailoc, a2=1.75, a3=0, a4=0, a5=0)
#             p_Order45 = SimpleNamespace(**result.best_values)
#             ex_Order45 = np.sqrt(np.diag(result.covar))
#             A1_loc_T.append(p_Order45.a1)
#             Err_A1_loc_T.append(ex_Order45[1])    

# #So now lets calculate the velocities
# A1_T = []

# for n in range(484):
#     A1 = ((A1_loc_T[n]-np.nanmean(A1_loc_T))*1.6273849607182894e-05)
#     A1_T.append(A1)

# Line_shift = np.reshape(A1_T, (22, 22))
# lamda_H3 = (3.95295 + 3.98558 + 3.994 + 3.971)/4

# Velocity_T = []

# for n in range(484):
#     V = ((A1_T[n]/(lamda_H3_O4+lamda_H3_O5))*c)/2
#     Velocity_T.append(V/1000)

# Velocities_T = np.reshape(Velocity_T, (22, 22))
# Velocities_T = Velocities_T.transpose()

# Velocity_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

# Err = (Err1 + Err3)/2
# Uranus_width_pix = (Uranus_width_pix_O4 + Uranus_width_pix_O5)/2

# ErrorA = (Err/Uranus_width_pix)

# for n in range(484):
#     Err = Velocity_T[n]*np.sqrt((Err_A1_loc_T[n]/A1_T[n])**2)
#     Velocity_Err.append(np.sqrt(Err**2)/1000)
       
# Velocity_Err = np.reshape(Velocity_Err, (22, 22))
# Velocity_Err = Velocity_Err.transpose()

# #%%
# k = 2

# TIME = ['07:30', '07:45', '08:00', '08:14', '08:29', '08:44', '08:59', '09:14', '09:29', '09:44', '09:59', '10:13', '10:28', '10:43', '10:58', '11:12', '11:27', '11:42',  '11:57', '12:11', '12:27', '12:42']
# for i in range(22):
#     k = i
#     fig, ax = plt.subplots(figsize=(12,8))
#     ax2 = ax.twinx()
#     ax.set_title('Approximate ion velocities from all $H_{3}^{+}$ lines over a ~15 minute exposure at ~' +str(TIME[k]) +' UT', fontsize=25)
#     ax.plot(np.arange(22), np.flip((Planet_rotation_O4+Planet_rotation_O5)/2), color='k', label='Planetary Rotation')
#     ax.plot(np.arange(22), np.flip(Velocities_T[:,k]), color='b', label='IR Velocities', lw = 3)
#     ax.fill_between(np.arange(22), np.flip(Velocities_T[:,k])-np.flip(Velocity_Err[:,k]), np.flip(Velocities_T[:,k])+np.flip(Velocity_Err[:,k]), color='b', alpha=0.5)
#     ax2.plot(np.arange(22), np.flip(Mapping_Q1A[10-k,:]), color='r', label='IR Intensity', lw = 3)
#     #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
#     # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
#     # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#     ax.set_xticks(np.linspace(0, 21, 7))
#     ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'])
#     # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
#     ax.set_xlabel('Rough Arbitrary Uranian Longitude ($^\circ$)', fontsize=20, labelpad=15)
#     ax.set_ylabel('Wind Speeds from all $H_{3}^{+}$ lines across Uranus (ULS) ($^\circ$)', fontsize=20, labelpad=15)
#     #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) +- ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr[50:100]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
#     ax.vlines(np.linspace(0, 21, 7), -200, 200, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#     ax.hlines(np.linspace(-150, 150, 7), 0, 21, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#     ax.tick_params(axis='both', labelsize=15)
#     ax2.tick_params(axis='both', labelsize=15)
#     ax.set_xlim(0, 21) #0, 3601
#     ax.set_ylim(-3, 3)
#     ax2.set_ylim(0, 1.0)
#     ax.legend(fontsize=15)
    #ax2.legend(loc='upper left', fontsize=15)
    #plt.title('Approximate ion velocities from $H_{3}^{+}$ $Q(1,0^{-})$ over a 30 to 40 minute exposure at ' +str(TIME[k]) +' UT', fontsize=30)
    
#%%
dataI = (AB_Combined_Order4[32]+AB_Combined_Order4[33]+AB_Combined_Order4[34]+AB_Combined_Order4[35]+AB_Combined_Order4[36]+AB_Combined_Order4[37]+AB_Combined_Order4[38]+AB_Combined_Order4[39])/16
dataI[dataI > 0.020] = 0.0025

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

for m in range(22):
    datai = dataI[16+m, :]
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
    a1_pointP = a1_pointP[0][0] - 1110
    ex_resultQ1P = gmodel.fit(datai[1110:1310], x=XX, a0=np.nanmax(datai[1205:1225]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    ex_pQ1P = SimpleNamespace(**ex_resultQ1P.best_values)
    ex_eQ1P = np.sqrt(np.diag(ex_resultQ1P.covar))
    A0 = ex_pQ1P.a0
    ex_A0_Q1.append(A0)
    A2 = ex_pQ1P.a2
    ex_A2_Q18.append((A2)*1.5709811506855093e-05)
    err_A0 = ex_eQ1P[0]
    ex_Err_A0_Q18.append(err_A0)
    err_A2 = ex_eQ1P[2]
    ex_Err_A2_Q18.append(err_A2*1.5709811506855093e-05)
    ex_A1_Q18.append(ex_pQ1P.a1)
    ex_Err_A1_Q18.append(ex_eQ1P[1]*1.5709811506855093e-05)
    ex_FWHM = A2*1.5709811506855093e-05*2*np.sqrt(2*math.log(2))
    ex_FWHMQ18.append(ex_FWHM)
    ex_INTSQ18.append(A0*ex_FWHM*(10**6))

A1_Q1 = []
A1_Q3 = []
A1_Q18 = []
A1_Q38 = []
A1_Q12 = []

for n in range(22):
    A1 = ((ex_A1_Q18[n]-np.nanmean(ex_A1_Q18))*1.5709811506855093e-05)
    A1_Q18.append(A1)
# for n in range(25):
#     A1 = ((ex_A1_Q38[n]-np.nanmean(ex_A1_Q38))*1.6361537623956223e-05)
#     A1_Q13.append(A1)    

#Now lets calculate the planetary rotation
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
Limb_velocity_pix_O5 = (Period_Uranus/1000) #Need to see if there's slowing or speeding up of atmospheric/ionospheric winds at Uranus

c = 299792458 
lamda_H3_O4 = 3.95295 
lamda_H3_O5 = 3.98558
Velocity_O4 = []
Velocity_O5 = []
Velocity_O48 = []
Velocity_O58 = []
Velocity_O42 = []

for n in range(22):
    V = ((A1_Q18[n]/lamda_H3_O4)*c)
    Velocity_O48.append(V/1000)
    
Velocity_O4_Err = []
Velocity_O48_Err = []
Velocity_O5_Err = []
Velocity_O58_Err = []
Velocity_O42_Err = []

Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)
   
for n in range(22):
    Err = Velocity_O48[n]*np.sqrt((ex_Err_A1_Q18[n]/A1_Q18[n])**2)
    Velocity_O48_Err.append(np.sqrt(Err**2))

Planet_rotation_O4 = np.linspace(Limb_velocity_pix_O4, -1*Limb_velocity_pix_O4, 22)
Planet_rotation_O5 = np.linspace(Limb_velocity_pix_O5, -1*Limb_velocity_pix_O5, 22)
Planet_rotation_O4_90 = np.linspace(Limb_velocity_pix_O4_90, -1*Limb_velocity_pix_O4_90, 22)
Planet_rotation_O4_80 = np.linspace(Limb_velocity_pix_O4_80, -1*Limb_velocity_pix_O4_80, 22)

PosErr4 = []
NegErr4 = []

for n in range(22):
    PosErr4.append(Vels_Set8[n]+Velocity_O48_Err[n])
    NegErr4.append(Vels_Set8[n]-Velocity_O48_Err[n])
    
Vels_Set8 = Velocity_O48
IntsSet8 = ex_INTSQ18
Errs_Set8 = Velocity_O48_Err

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~60 minute exposure between 11:00 and 12:00 UTC on $11^{th}$ October 2016', fontsize=22, pad = 25)
ax.plot(Mean_Longs, Planet_rotation_O4, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Mean_Longs, Planet_rotation_O4_90, color='darkgreen', ls = '--', label='Reduced Planetary Rotation by ~(257m/s)')
#ax.plot(Mean_Longs, Planet_rotation_O4_80, color='lightgreen', ls = '--', label='Reduced Planetary Rotation by ~(514m/s)')
ax.plot(Mean_Longs, Vels_Set8, color='b', ls='--', label='Q(1,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
ax.fill_between(Mean_Longs, NegErr4, PosErr4, color='b', alpha=0.5)
#ax.fill_between(np.arange(25), NegErrQ3 - Planet_rotation_O58, PosErrQ3 - Planet_rotation_O58, color='b', alpha=0.5)
ax2.plot(Mean_Longs, np.flip(ex_INTSQ18), color='r', label='IR Intensity', lw = 5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), fontsize=15)
ax.set_xticklabels(labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=15)
ax.tick_params(axis='both', which='major', labelsize=15)
ax2.tick_params(axis='both', which='major', labelsize=15)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=20)
ax.set_ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ $Q(1,0^{-})$ (kms$^{-1}$)', fontsize=20)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(Q1_IntErr), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=20)
ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-3, 3)
ax2.set_ylim(0, 1.0)
ax.legend(loc='lower left', fontsize=20)
ax2.legend(loc='lower right', fontsize=20)

#%% Need to add Q1 and Q3 together to do that translate one to the other and see if other lines are possible to enhance signal!