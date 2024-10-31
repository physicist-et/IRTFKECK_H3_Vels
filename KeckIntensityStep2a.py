# -*- coding: utf-8 -*-
"""
Created on Wed Jan 26 15:55:53 2022

@author: snowy
"""

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits #import the relevant directories to read a fits file from your directory and plot it
from scipy.optimize import curve_fit
from KeckDataReductionStep1 import Cal_Flats, ABBASep
import scipy.ndimage as f 
import warnings
from matplotlib.patches import Ellipse
import h3ppy
import math
from lmfit import Model
from types import SimpleNamespace

#%%
image_file1 = 'C:/Users/snowy/OneDrive/Documents/Python work/Keck 13OCT14/order1/order1s0041.fits.gz' #imports a single fits file from the cwd + the file path in the commas
fits.open(image_file1)

image_data1 = fits.getdata(image_file1, ext=0) #Focuses specifically on the array data contained within the fits file

Keck_Data = [] #First create lists in which the arrays of Keck data will go into
Keck_DataABBA = []
Keck_Data_Total = np.zeros((image_data1.shape[0], image_data1.shape[1]))
s = 0 #This will be used to create the ABBA pattern

for n in range(72): #We use this list to create a list which holds all the data from Order19
    num = n + 45
    if num < 100:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/order1/order1s00' + str(num) + '.fits.gz'
        fits.open(image_filei)
        image_datai = fits.getdata(image_filei, ext=0)
        Keck_Data.append(image_datai)
    elif num >= 100:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/order1/order1s0' + str(num) + '.fits.gz'
        fits.open(image_filei)
        image_datai = fits.getdata(image_filei, ext=0)
        Keck_Data.append(image_datai)
    else:
        pass
    
ABsubs = len(Keck_Data)//4 #We create this variable to create the ABBA pattern

for n in range(ABsubs): #This for loop uses the A - B - B + A and adds them into a list Keck_DataABBA before creating a total spectra 
    image_dataSUB = ((Keck_Data[s]/30 - Cal_Flats) - (Keck_Data[(s+1)]/30 - Cal_Flats) - (Keck_Data[(s+2)]/30 - Cal_Flats) + (Keck_Data[(s+3)]/30 - Cal_Flats))/4
    #image_dataSUB = ((Keck_Data[s]/30) - (Keck_Data[(s+1)]/30) - (Keck_Data[(s+2)]/30) + (Keck_Data[(s+3)]/30))/4
    Keck_DataABBA.append(image_dataSUB)
    s += 4
    Keck_Data_Total += Keck_DataABBA[n]
    #filename = 'KeckDataABBASet' + str(n) + '.txt'
    #np.savetxt(filename, Keck_DataABBA[n])
    
#%% Add the corrction here

corrective_x = np.load('Corrective_Array.npy')
IRTF_Data_Q1 = []
IRTF_Data_Q3 = []

for n in range(ABsubs):
    if n == 0:
        for i in range(183):
            Fig = Keck_DataABBA[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][144:151]) #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q1 = shift_Spec_row
    else:
        for i in range(183):
            Fig = Keck_DataABBA[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][144:151])  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q1 = np.dstack((IRTF_Data_Q1, shift_Spec_row))
        
for n in range(ABsubs):
    if n == 0:
        Total_IRTF_DataQ1 = IRTF_Data_Q1[:,:,0]
    else:
        Total_IRTF_DataQ1 = np.add(Total_IRTF_DataQ1, IRTF_Data_Q1[:,:,n])

for n in range(ABsubs):
    if n == 0:
        for i in range(183):
            Fig = Keck_DataABBA[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][694:700]) #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q3 = shift_Spec_row
    else:
        for i in range(183):
            Fig = Keck_DataABBA[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][694:700])  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q3 = np.dstack((IRTF_Data_Q3, shift_Spec_row))
        
for n in range(ABsubs):
    if n == 0:
        Total_IRTF_DataQ3 = IRTF_Data_Q3[:,:,0]
    else:
        Total_IRTF_DataQ3 = np.add(Total_IRTF_DataQ3, IRTF_Data_Q3[:,:,n])
        
#%% Lets do Q2 and Q32 if possible!!
#Q2
IRTF_Data_Q2 = []
IRTF_Data_Q32 = []
for n in range(ABsubs):
    if n == 0:
        for i in range(183):
            Fig = Keck_DataABBA[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][445:451]) #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q2 = shift_Spec_row
    else:
        for i in range(183):
            Fig = Keck_DataABBA[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][445:451])  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q2 = np.dstack((IRTF_Data_Q2, shift_Spec_row))
        
for n in range(ABsubs):
    if n == 0:
        Total_IRTF_DataQ2 = IRTF_Data_Q2[:,:,0]
    else:
        Total_IRTF_DataQ2 = np.add(Total_IRTF_DataQ2, IRTF_Data_Q2[:,:,n])

#Q32
for n in range(ABsubs):
    if n == 0:
        for i in range(183):
            Fig = Keck_DataABBA[0]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][853:860]) #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row)) 
        IRTF_Data_Q32 = shift_Spec_row
    else:
        for i in range(183):
            Fig = Keck_DataABBA[n]
            Spec_row = Fig[i, :]
            corrective_shift_factor = np.nanmean(corrective_x[i][853:860])  #Not sure why this isn't working
            shifted_Spec_row = f.shift(Spec_row, shift=[corrective_shift_factor], mode='wrap')
            if i == 0:
                shift_Spec_row = shifted_Spec_row
            else:
                shift_Spec_row = np.vstack((shift_Spec_row, shifted_Spec_row))
        IRTF_Data_Q32 = np.dstack((IRTF_Data_Q32, shift_Spec_row))
        
for n in range(ABsubs):
    if n == 0:
        Total_IRTF_DataQ32 = IRTF_Data_Q32[:,:,0]
    else:
        Total_IRTF_DataQ32 = np.add(Total_IRTF_DataQ32, IRTF_Data_Q32[:,:,n])
        
#%% Now lets calibrate the data #HERE TO COMBO DATA AND SEE WHAT THE VELOCITIES LOOK LIKE! DO DARKS AND POS'

Total_IRTF_DataQ1A = np.nanmean((IRTF_Data_Q1[:,:,0:12]), axis=(2))

A1_P = []
E1_P = []
A1_N = []
E1_N = []

def gaussian_fit(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

gmodel = Model(gaussian_fit)

for a in range(34):
    EM_1P = Total_IRTF_DataQ1A[28+a, 125:175]
    A01P = np.nanmax(EM_1P[20:30])
    a1_pointO1 = np.where(EM_1P == np.nanmax(EM_1P[20:30]))
    a1_pointO1 = a1_pointO1[0][0]
    resultQ1P = gmodel.fit(EM_1P, x=np.arange(50), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
    pQ1P = SimpleNamespace(**resultQ1P.best_values)
    eQ1P = np.sqrt(np.diag(resultQ1P.covar))
    A1_P.append(pQ1P.a1)
    E1_P.append(eQ1P[1])

wave = np.load('Wavelength_Order.npy')
Shiftings  = []
Vels_roughP = []
Errs_P = []

for a in range(34):
    Total_shifts = (A1_P[33-a]-np.nanmean(A1_P[2:-2]))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.5295))*299792.458
    Vels_roughP.append(Vels)
    Errs_P.append(np.sqrt((Vels*(np.sqrt((E1_P[a]*(np.nanmean(np.gradient(wave[145:150])))/Total_shifts)**2+(np.std(A1_P)/np.nanmean(A1_P))**2)))**2))

for a in range(33):
    EM_1P = Total_IRTF_DataQ1A[111+a, 125:175]
    A01P = np.nanmin(EM_1P[20:30])
    a1_pointO1 = np.where(EM_1P == np.nanmin(EM_1P[20:30]))
    a1_pointO1 = a1_pointO1[0][0]
    resultQ1P = gmodel.fit(EM_1P, x=np.arange(50), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
    pQ1P = SimpleNamespace(**resultQ1P.best_values)
    eQ1P = np.sqrt(np.diag(resultQ1P.covar))
    A1_N.append(pQ1P.a1)
    E1_N.append(eQ1P[1])
    
Shiftings  = []
Vels_roughN = []
Errs_N = []

for a in range(33):
    Total_shifts = (A1_N[32-a]-np.nanmean(A1_N[2:-2]))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.5295))*299792.458
    Vels_roughN.append(Vels)
    Errs_N.append(np.sqrt((Vels*(np.sqrt((E1_N[a]*(np.nanmean(np.gradient(wave[145:150])))/Total_shifts)**2+(np.std(A1_N)/np.nanmean(A1_N))**2)))**2))

#%%
# Line_Pos = Total_IRTF_DataQ1A[25:65,:]
# Line_Neg = Total_IRTF_DataQ1A[107:147,:]
import matplotlib.ticker as tkr

#Shift the positive line to the negative line
Pos_Moved_Spectra = f.shift(Total_IRTF_DataQ1A, shift=(np.nanmean(ABBASep[145:150]), 0), mode='wrap')

Fin_IRTF_Q1A = (-1*Total_IRTF_DataQ1A + Pos_Moved_Spectra)/2
#Fin_IRTF_Q1A = Fin_IRTF_Q1A[84:129]

plt.figure()
plt.imshow(np.flipud(Fin_IRTF_Q1A[85:-5,100:900]), cmap='gist_heat', vmax = 1.5, vmin = 0)
plt.title('Averaged $H_{3}^{+}$ emission spectrum for Keck II NIRSPEC in October 2014', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.2f'))
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('CCD counts (counts per second)', fontsize = 25)
plt.xlabel('Wavelength ($\mu$m)', fontsize = 25)
plt.ylabel('Spatial axis row position (pixels)', fontsize = 25)
plt.yticks(fontsize = 20)
#plt.xticks(ticks = (0, 799), labels=("{:.3f}".format(round(wave[100], 2)), "{:.3f}".format(round(wave[899], 2))), fontsize = 20)
plt.vlines((40, 55, 590, 605), ymin = 0, ymax = 183, color='w', lw = 3, ls = 'dashed')
plt.ylim(0, 85)
plt.text(22, 28, '$Q(1,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(329, 28, '$Q(2,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(574, 28, '$Q(3,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(739, 28, '$Q(3,2^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(630, 28, '$Q(3,1^{-})$', fontsize = 15, color = 'w', rotation = 90)

#%% Extract out the line emission
Fc_HR718 = np.load('HR_718Flux.npy')
wave = np.load('Wavelength_Order.npy')
image = np.rot90(np.flipud(Fin_IRTF_Q1A[85:-5,135:161]), k = 1)*Fc_HR718[135:161, None]*((4.2545*(10**10))/(0.288*0.1445))
image2 = np.rot90(np.flipud(Fin_IRTF_Q1A[84:-6,684:710]), k = 1)*Fc_HR718[684:710, None]*((4.2545*(10**10))/(0.288*0.1445))

plt.figure(figsize = (16,7))
plt.imshow(image*10**3, cmap='gist_heat', vmax = 2, vmin = 0.00)
plt.title('a) Averaged $H_{3}^{+}$ emission line', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.2f'), shrink = 0.5)
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Intensity (m$Wm^{-2}sr^{-1}$$\mu$m$^{-1}$)', fontsize = 25)
labels = str("{:.5f}".format(wave[0])), str("{:.5f}".format(wave[5])), str("{:.5f}".format(wave[10])), str("{:.5f}".format(wave[15])), str("{:.5f}".format(wave[20])), str("{:.5f}".format(wave[25]))
plt.yticks([0, 5, 10, 15, 20, 25], labels = labels, fontsize=20)
plt.xticks(fontsize = 20)
plt.ylabel('Wavelength ($\mu$m)', fontsize = 25, labelpad = 10)
plt.xlabel('Spatial axis row position (pixels)', fontsize = 25)

#%%
def gaussian_fit(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y
import random 
# gmodel = Model(gaussian_fit)
# Int_Ex = []

# for a in range(94):
#     if a < 33:
#         Int_Ex.append((random.randrange(0, 10)*0.000000004))
#     elif a > 67:
#         Int_Ex.append((random.randrange(0, 10)*0.000000004))
#     else:
#         EM_1P = image[:,a]
#         A01P = np.nanmax(EM_1P[10:16])
#         a1_pointO1 = np.where(EM_1P == np.nanmax(EM_1P[10:16]))
#         a1_pointO1 = a1_pointO1[0][0]
#         resultQ1P = gmodel.fit(EM_1P, x=np.arange(26), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
#         pQ1P = SimpleNamespace(**resultQ1P.best_values)
#         eQ1P = np.sqrt(np.diag(resultQ1P.covar))
#         FWHM = pQ1P.a2*5.705793605517151e-05*2*np.sqrt(2*math.log(2))
#         Int_Ex.append(FWHM*pQ1P.a0)

# Int_Ex3 = []
# for a in range(94):
#     if a < 33:
#         Int_Ex3.append((random.randrange(0, 10)*0.000000004))
#     elif a > 67:
#         Int_Ex3.append((random.randrange(0, 10)*0.000000004))
#     else:
#         EM_1P = image2[:,a]
#         A01P = np.nanmax(EM_1P[10:16])
#         a1_pointO1 = np.where(EM_1P == np.nanmax(EM_1P[10:16]))
#         a1_pointO1 = a1_pointO1[0][0]
#         resultQ1P = gmodel.fit(EM_1P, x=np.arange(26), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
#         pQ1P = SimpleNamespace(**resultQ1P.best_values)
#         eQ1P = np.sqrt(np.diag(resultQ1P.covar))
#         FWHM = pQ1P.a2*5.705793605517151e-05*2*np.sqrt(2*math.log(2))
#         Int_Ex3.append(FWHM*pQ1P.a0)

# Ratio = np.array(Int_Ex[33:67])/np.array(Int_Ex3[33:67])
# P1 = np.zeros(33)
# P2= np.zeros(27)
# Ratio = np.concatenate((P1, Ratio))
# Ratio = np.concatenate((Ratio, P2))

# fig, ax = plt.subplots(figsize=(14,10))
# ax2 = ax.twinx()
# ax.plot(np.arange(94), np.array(Int_Ex)/(10**-6), color = 'gray', linewidth = 3, linestyle = '--', label = 'Q1')
# ax.plot(np.arange(94), np.array(Int_Ex3)/(10**-6), color = 'k', linewidth = 3, label = 'Q3')
# ax2.plot(np.arange(94), Ratio, color = 'c', linewidth = 2, label = 'Q1 to Q3 ratio')
# labels = '85', '105', '125', '145', '165'
# ax.set_xticks([0, 20, 40, 60, 80], labels = labels)
# ax.tick_params(axis='both', labelsize = 20)
# ax2.tick_params(axis='both', labelsize = 20)
# ax.set_ylabel('Intensity (m$Wm^{-2}sr^{-1}$$\mu$m$^{-1}$)', fontsize = 25)
# ax2.set_ylabel('Q(1,0$^{-}$) to Q(3,0$^{-}$) emission line ratio', fontsize = 25, labelpad = 10)
# ax.set_xlabel('Spatial Axis row position (0.144" per pixel)', fontsize = 25)
# ax.grid(alpha=0.5, linestyle='--')
# ax.set_xlim(0, 93)
# ax2.set_ylim(0, 4)
# ax.legend(fontsize= 20, loc = 'upper left')
# ax2.legend(fontsize= 20, loc = 'upper right')


# #%% Now lets see what LOS does to it
# r_planet = 32.560553633218/2
# Pixels = np.arange(-17, 17.5, 1)
# LOS = []
# for i in range(len(Pixels)):
#     if Pixels[i] > r_planet:
#         LOS.append(np.cos(1))
#     elif Pixels[i] < -16:
#         LOS.append(np.cos(1))
#     else:
#         LOS.append(np.cos(Pixels[i]/r_planet))

# fig, ax = plt.subplots(figsize=(14,10))
# ax2 = ax.twinx()
# ax.plot(np.arange(28,67), np.array(Int_Ex[30:69])/(10**-6), color = 'gray', linewidth = 3, linestyle = '--', label = 'Intensity')
# ax2.plot(np.arange(30, 65), np.array(LOS), color = 'r', linewidth = 2, label = 'LOS correction')
# Ratio = np.array(LOS)*np.array(Int_Ex[32:67])/(10**-6)
# ax.plot(np.arange(30, 65), Ratio, color = 'k', linewidth = 2, label = 'LOS corrected Intensity')
# ax.grid(alpha=0.5, linestyle='--')
# ax.set_ylim(0, 0.6)
# ax2.set_ylim(0, 1.1)
# ax.set_xlim(29, 65)
# labels = ['-1.0', '-0.5', '0.0', '0.5', '1.0']
# ax.set_xticks([30, 38.5, 47, 55.5, 64], labels = labels)
# ax.tick_params(axis='both', labelsize = 20)
# ax2.tick_params(axis='both', labelsize = 20)
# ax.set_ylabel('Intensity (m$Wm^{-2}sr^{-1}$$\mu$m$^{-1}$)', fontsize = 25)
# ax2.set_ylabel('LOS correction factor', fontsize = 25, labelpad = 10)
# ax.set_xlabel('Distance from central point of Uranus (R$_{U}$)', fontsize = 25)

# #%% Now Temps
# h = 6.63*(10**-34)
# c = 2.99*(10**8)
# kb = 1.38*(10**-23)

# y1 = 4*7*h*c*2529.73*128.7
# y2 = 4*19*h*c*2509.08*123.2
# y = y1/y2

# P1 = np.log(y) - np.log(Ratio)
# P2 = ((2552.57-2961.84)*100*h*c)/kb

# Temperatures = P2/P1
# Temperatures[Temperatures > 600] = 600

# fig, ax = plt.subplots(figsize=(14,10))
# ax2 = ax.twinx()
# ax.plot(np.arange(94), np.array(Int_Ex)/(10**-6), color = 'gray', linewidth = 3, linestyle = '--', label = 'Q1')
# ax.plot(np.arange(94), np.array(Int_Ex3)/(10**-6), color = 'k', linewidth = 3, label = 'Q3')
# ax2.plot(np.arange(94), Temperatures, color = 'r', linewidth = 2, label = 'Temperature')
# labels = '85', '105', '125', '145', '165'
# ax.set_xticks([0, 20, 40, 60, 80], labels = labels)
# ax.tick_params(axis='both', labelsize = 20)
# ax2.tick_params(axis='both', labelsize = 20)
# ax.set_ylabel('Intensity (m$Wm^{-2}sr^{-1}$$\mu$m$^{-1}$)', fontsize = 25)
# ax2.set_ylabel('Temperature (K)', fontsize = 25, labelpad = 10)
# ax.set_xlabel('Spatial Axis row position (0.144" per pixel)', fontsize = 25)
# ax.grid(alpha=0.5, linestyle='--')
# ax.set_xlim(0, 93)
# ax2.set_ylim(350, 800)
# ax.legend(fontsize= 20, loc = 'upper left')
# ax2.legend(fontsize= 20, loc = 'upper right')
#%%
def gaussian_fit(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

gmodel = Model(gaussian_fit)

A1_P = []
E1_P = []
A1_N = []
E1_N = []

for a in range(33):
    EM_1P = Fin_IRTF_Q1A[111+a, 125:175]
    A01P = np.nanmax(EM_1P[20:30])
    a1_pointO1 = np.where(EM_1P == np.nanmax(EM_1P[20:30]))
    a1_pointO1 = a1_pointO1[0][0]
    resultQ1P = gmodel.fit(EM_1P, x=np.arange(50), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
    pQ1P = SimpleNamespace(**resultQ1P.best_values)
    eQ1P = np.sqrt(np.diag(resultQ1P.covar))
    A1_P.append(pQ1P.a1)
    E1_P.append(eQ1P[1])

wave = np.load('Wavelength_Order.npy')
Shiftings  = []
Vels_rough = []
Q1_Velocities = []

for a in range(33):
    Total_shifts = (A1_P[32-a]-np.nanmean(A1_P))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.53006))*299792.458
    Vels_rough.append(Vels)

Q1_Velocities.append(Vels_rough)
plt.figure()
plt.plot(np.arange(33), np.flip(Vels_rough[:]), 'r')
#plt.fill_between(np.arange(30), Vels_rough[:]+np.array(E1_P), Vels_rough[:]-np.array(E1_P), 'r', alpha = 0.5)

# Shiftings  = []
# Vels_rough = []

# for a in range(31):
#     Total_shifts = -1*(A1_N[a]-np.nanmean(A1_N))*np.nanmean(np.gradient(wave[145:150]))
#     Shiftings.append(Total_shifts)
#     Vels = (Total_shifts/(3.53006))*299792.458
#     Vels_rough.append(Vels)

# plt.plot(np.arange(31), np.flip(Vels_rough[:]), 'b')
# Q1_Velocities.append(Vels_rough)
#plt.fill_between(np.arange(30), Vels_rough[:]+np.array(E1_N), Vels_rough[:]-np.array(E1_N), 'b', alpha = 0.5)
#%%
import matplotlib.ticker as tkr
Total_IRTF_DataQ3A = np.nanmean((IRTF_Data_Q3[:,:,0:12]), axis=(2))

# Line_Pos = Total_IRTF_DataQ1A[25:65,:]
# Line_Neg = Total_IRTF_DataQ1A[107:147,:]

#Shift the positive line to the negative line
Pos_Moved_Spectra = f.shift(Total_IRTF_DataQ3A, shift=(np.nanmean(ABBASep[695:700])+0.5, 0), mode='wrap')

Fin_IRTF_Q3A = (-1*Total_IRTF_DataQ3A + Pos_Moved_Spectra)/2
#Fin_IRTF_Q1A = Fin_IRTF_Q1A[84:129]

plt.figure()
plt.imshow(np.flipud(Fin_IRTF_Q3A[85:,100:900])/72, cmap='gist_heat', vmax = 0.02, vmin = 0)
plt.title('c) Averaged $H_{3}^{+}$ emission spectrum for Keck II NIRSPEC in October 2014', fontsize = 30)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.2f'))
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Intensity ${\mu}Wm^{-2}sr^{-1}$', fontsize = 25)
plt.xlabel('Wavelength ($\mu$m)', fontsize = 25)
plt.ylabel('Spatial axis row position (pixels)', fontsize = 25)
plt.yticks(fontsize = 20)
plt.xticks(ticks = (0, 799), labels=("{:.3f}".format(round(wave[100], 2)), "{:.3f}".format(round(wave[899], 2))), fontsize = 20)
plt.vlines((40, 55, 590, 605), ymin = 0, ymax = 183, color='w', lw = 3, ls = 'dashed')
plt.ylim(0, 85)
plt.text(22, 28, '$Q(1,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(329, 28, '$Q(2,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(574, 28, '$Q(3,0^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(739, 28, '$Q(3,2^{-})$', fontsize = 15, color = 'w', rotation = 90)
plt.text(630, 28, '$Q(3,1^{-})$', fontsize = 15, color = 'w', rotation = 90)

A1_P = []
E1_P = []
A1_N = []
E1_N = []

for a in range(33):
    EM_1P = Fin_IRTF_Q3A[111+a, 125:175]
    A01P = np.nanmax(EM_1P[20:30])
    a1_pointO1 = np.where(EM_1P == np.nanmax(EM_1P[20:30]))
    a1_pointO1 = a1_pointO1[0][0]
    resultQ1P = gmodel.fit(EM_1P, x=np.arange(50), a0=A01P, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
    pQ1P = SimpleNamespace(**resultQ1P.best_values)
    eQ1P = np.sqrt(np.diag(resultQ1P.covar))
    A1_P.append(pQ1P.a1)
    E1_P.append(eQ1P[1])
    
Shiftings  = []
Vels_rough = []
Q3_Velocities = []

for a in range(33):
    Total_shifts = (A1_P[32-a]-np.nanmean(A1_P))*np.nanmean(np.gradient(wave[695:700]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.9860))*299792.458
    Vels_rough.append(Vels)

Q3_Velocities.append(Vels_rough)
plt.figure()
plt.plot(np.arange(33), np.flip(Vels_rough[:]), 'r')

# Shiftings  = []
# Vels_rough = []

# for a in range(33):
#     Total_shifts = -1*(A1_N[a]-np.nanmean(A1_N[5:-5]))*np.nanmean(np.gradient(wave[145:150]))
#     Shiftings.append(Total_shifts)
#     Vels = (Total_shifts/(3.95295))*299792.458
#     Vels_rough.append(Vels)

# plt.plot(np.arange(1, 33), np.flip(Vels_rough[:]), 'g')

#%%
Shape_Order1 = np.shape(image_data1)
Order_1_Offset = ABBASep

COADD_ABBA = []
COADD_ABBA_Total = np.zeros((75, Shape_Order1[1]))

Fc_HR718 = np.load('HR_718Flux.npy')

for xx in range(ABsubs):
    ABBA = IRTF_Data_Q1[:,:,xx]*Fc_HR718*((4.2545*(10**10))/(0.288*0.1445))
    for n in range(Shape_Order1[1]):
        Spec_col = ABBA[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_1_Offset[n]], mode='wrap')
        if n == 0:
            Order_1_shift = shifted_Spec_col
        else:
            Order_1_shift = np.vstack((Order_1_shift, shifted_Spec_col))

    Order_1_shift = np.flipud(np.rot90(Order_1_shift))
    FIN_ABBA = (ABBA - Order_1_shift)/2
    COADD_ABBA.append(FIN_ABBA[0:75,:])
    COADD_ABBA_Total += COADD_ABBA[xx]
    
COADD_ABBAQ3 = []
COADD_ABBA_TotalQ3 = np.zeros((75, Shape_Order1[1]))

for xx in range(ABsubs):
    ABBA = IRTF_Data_Q3[:,:,xx]*Fc_HR718*((4.2545*(10**10))/(0.288*0.1445))
    for n in range(Shape_Order1[1]):
        Spec_col = ABBA[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_1_Offset[n]], mode='wrap')
        if n == 0:
            Order_1_shift = shifted_Spec_col
        else:
            Order_1_shift = np.vstack((Order_1_shift, shifted_Spec_col))

    Order_1_shift = np.flipud(np.rot90(Order_1_shift))
    FIN_ABBA = (ABBA - Order_1_shift)/2
    COADD_ABBAQ3.append(FIN_ABBA[0:75,:])
    COADD_ABBA_TotalQ3 += COADD_ABBAQ3[xx]

COADD_ABBAQ2 = []
COADD_ABBA_TotalQ2 = np.zeros((75, Shape_Order1[1]))

for xx in range(ABsubs):
    ABBA = IRTF_Data_Q2[:,:,xx]*Fc_HR718*((4.2545*(10**10))/(0.288*0.1445))
    for n in range(Shape_Order1[1]):
        Spec_col = ABBA[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_1_Offset[n]], mode='wrap')
        if n == 0:
            Order_1_shift = shifted_Spec_col
        else:
            Order_1_shift = np.vstack((Order_1_shift, shifted_Spec_col))

    Order_1_shift = np.flipud(np.rot90(Order_1_shift))
    FIN_ABBA = (ABBA - Order_1_shift)/2
    COADD_ABBAQ2.append(FIN_ABBA[0:75,:])
    COADD_ABBA_TotalQ2 += COADD_ABBA[xx]

COADD_ABBAQ32 = []
COADD_ABBA_TotalQ32 = np.zeros((75, Shape_Order1[1]))

for xx in range(ABsubs):
    ABBA = IRTF_Data_Q32[:,:,xx]*Fc_HR718*((4.2545*(10**10))/(0.288*0.1445))
    for n in range(Shape_Order1[1]):
        Spec_col = ABBA[:,n]     #This selects a column at a time to shift up in the y direction of Order 4
        shifted_Spec_col = f.shift(Spec_col, shift=[-1*Order_1_Offset[n]], mode='wrap')
        if n == 0:
            Order_1_shift = shifted_Spec_col
        else:
            Order_1_shift = np.vstack((Order_1_shift, shifted_Spec_col))

    Order_1_shift = np.flipud(np.rot90(Order_1_shift))
    FIN_ABBA = (ABBA - Order_1_shift)/2
    COADD_ABBAQ32.append(FIN_ABBA[0:75,:])
    COADD_ABBA_TotalQ32 += COADD_ABBA[xx]
# for a in range(ABsubs):
#     plt.figure()
#     plt.imshow(COADD_ABBA[a], cmap='gist_gray', vmax = 0.025, vmin = -0.025)
#     plt.title(r'First ABBA observation set of Uranus with NIRSPEC on $12^{th}$ October 2014', fontsize=23)
#     plt.xlabel(r'Wavelength ($\mu$m)', fontsize=20)
#     #label_x = 0, 200, 400, 600, 800, 1000
#     #plt.xticks(label_x, ('3.9445', '3.9565', '3.9685', '3.9804', '3.9924', '4.0044'), fontsize=15)
#     plt.ylabel('Spatial position across Slit (Pixel No.)', fontsize=16) 
#     cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
#     #cbar.set_label('Uncalibrated Intensity (CCD Counts)', fontsize=20)
    
#%% Next is to work out the intensities and produce a simple colour graph of the slit changing with time (and confirm the orientation of the slit)

SLIT_ANG = []
warnings.filterwarnings('ignore', 'The following header keyword is invalid or follows an unrecognized non-standard convention',)

#Find the orientation of the slit
for n in range(72): #We use this list to create a list which holds all the data from Order19
    num = n + 45
    if num < 100:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s00' + str(num) + '.fits'
        hdul = fits.open(image_filei, ignore_missing_end=True)
        hdr = hdul[0].header
        SLIT_ANG.append(hdr['SLITANG'])
        #print(hdr['UTC'])
    elif num >= 100:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s0' + str(num) + '.fits'
        hdul = fits.open(image_filei, ignore_missing_end=True)
        hdr = hdul[0].header
        SLIT_ANG.append(hdr['SLITANG'])
        #print(hdr['UTC'])
    else:
        pass

uranus_seangleK = (27.352772 + 27.338911)/2
#stretch yy to become a sphere
flattening = 0.0229
uranus_seangleK_Rads = (uranus_seangleK*np.pi)/180
losflattening=flattening*(1-np.sin(uranus_seangleK_Rads))
eq_po_ratioK=1-losflattening

#Slit is 8.1 degrees so slightly off from E-W alignment
fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
ax.imshow(COADD_ABBA_Total/ABsubs, cmap='gist_gray', vmax = 0.005, vmin = -0.005)
EllsQ1 = Ellipse((147.5, 44.5), 31, 33, angle=0, linewidth=1, color = 'b', fill=False)
EllsQ3 = Ellipse((697, 44.5), 31, 33, angle=0, linewidth=1, color ='g', fill=False)
ax.add_artist(EllsQ1)
ax.add_artist(EllsQ3)

fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
ax.imshow(COADD_ABBA_TotalQ3/ABsubs, cmap='gist_gray', vmax = 0.005, vmin = -0.005)
EllsQ1 = Ellipse((147.5, 44.5), 31, 33, angle=0, linewidth=1, color = 'b', fill=False)
EllsQ3 = Ellipse((697, 44.5), 31, 33, angle=0, linewidth=1, color ='g', fill=False)
ax.add_artist(EllsQ1)
ax.add_artist(EllsQ3)

#%% Now lets make a picture placing this onto a fake Uranus
image = 1000*np.array(COADD_ABBA_Total[14:,130:166])/ABsubs
image[image < 0] = 0

fig, ax = plt.subplots(figsize=(12,9))
im = ax.imshow(np.rot90(image, k=3), cmap='gist_heat', vmax = np.nanmean(image)+(10*np.nanstd(image)), vmin = 0)
EllsQ1 = Ellipse((30, 17.28), 32, 31, angle=0, linewidth=4, color = 'c', linestyle = '--', fill=False)
EllsQ2 = Ellipse((30, 17.28), 28.65, 28, angle=0, linewidth=4, color = 'm', linestyle = '--', fill=False)
EllsQ3 = Ellipse((30, 17.28), 25.6, 25, angle=0, linewidth=4, color ='g', linestyle = '--', fill=False)
ax.add_artist(EllsQ1)
ax.add_artist(EllsQ2)
ax.add_artist(EllsQ3)
ax.set_ylabel('Spatial Axis (arc seconds)', fontsize = 25)
ax.set_xlabel('Spatial Axis (arc seconds)', fontsize = 25)
cbar = fig.colorbar(im, format=tkr.FormatStrFormatter('%.2f'))
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Spectral Radiance (mWm$^{-2}$sr$^{-1}$${\mu}$m$^{-1}$)', fontsize = 25)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.grid(color = 'w', linestyle = '--', alpha = 0.5)
ax.set_xlim(0,60)
ax.set_ylim(0, 35)

#%%
#Rough size for Uranus is 25.636988156280573 pixels (probably more need to check with script)
#Lets use h3ppy with wave to find out the values for the first scan

XX = np.arange(300)
X = np.arange(100)
AVALS = []
COVALS = []
SDs = []

ABBAi = []
TLDR = []
i = 0

def gaussian_fit(x, a0, a1, a2, a3, a4, a5):
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

wave = np.load('Wavelength_Order.npy')
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)
for xx in range(6):
    ABBAi = (COADD_ABBA[xx*2] + COADD_ABBA[(xx*2)+1])/2
    plt.figure()
    plt.imshow(ABBAi, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    # if xx == 1:
    #     ABBAi[42,146] = np.nanmean(ABBAi) - 1e-05
    #     # ABBAi[45:52,149:152] = np.nanmean(ABBAi) - 1e-05
    #     ABBAi[29:32,150:152] = np.nanmean(ABBAi) - 1e-05
    #     # ABBAi[32:36,146] = np.nanmean(ABBAi) - 1e-05
    #     plt.figure()
    #     plt.imshow(ABBAi, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    for i in range(31): #28
        if xx == 0:
            ii = i + 27
        else:
            ii = i + 28
        ABBAt = (ABBAi[ii,:])
        A01 = np.nanmax(ABBAt[135:155])
        a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[135:155]))
        a1_pointO1 = a1_pointO1[0][0] - 100
        resultQ1 = gmodel.fit(ABBAt[100:200], x=np.arange(100), a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
        pQ1 = SimpleNamespace(**resultQ1.best_values)
        eQ1 = np.sqrt(np.diag(resultQ1.covar))
        A0_Q1.append(pQ1.a0)
        A1_Q1.append(pQ1.a1)
        A2_Q1.append(pQ1.a2*6.119575999463667e-05)
        Err_A0_Q1.append(eQ1[0])
        Err_A1_Q1.append(eQ1[1]*6.119575999463667e-05)
        Err_A2_Q1.append(eQ1[2]*6.119575999463667e-05)
        FWHM = pQ1.a2*6.119575999463667e-05*2*np.sqrt(2*math.log(2))
        FWHMQ1.append(FWHM)
        INTSQ1.append(pQ1.a0*FWHM)
    ABBAi = (COADD_ABBAQ3[xx*2] + COADD_ABBAQ3[(xx*2)+1])/2    
    for i in range(31): # 28
        if xx == 5:
            ii = i + 28
        else:
            ii = i + 28 #30
        ABBAt = (ABBAi[ii,:])
        A03 = np.nanmax(ABBAt[690:705])
        a1_pointO3 = np.where(ABBAt == np.nanmax(ABBAt[690:705]))
        a1_pointO3 = a1_pointO3[0][0] - 680
        resultQ3 = gmodel.fit(ABBAt[680:720], x=np.arange(40), a0=A03, a1=a1_pointO3, a2=1.3, a3=0, a4=0, a5=0)
        pQ3 = SimpleNamespace(**resultQ3.best_values)
        eQ3 = np.sqrt(np.diag(resultQ3.covar))
        A0_Q3.append(pQ3.a0)
        A1_Q3.append(pQ3.a1)
        A2_Q3.append(pQ3.a2*5.705793605517151e-05)
        Err_A0_Q3.append(eQ3[0])
        Err_A1_Q3.append(eQ3[1]*5.705793605517151e-05) #If you look for the gradient on wave at Q1 and Q3 it differs!
        Err_A2_Q3.append(eQ3[2]*5.705793605517151e-05)
        FWHM = pQ3.a2*5.705793605517151e-05*2*np.sqrt(2*math.log(2))
        FWHMQ3.append(FWHM)
        INTSQ3.append(pQ3.a0*FWHM)
        # if xx == 0:
        # # #     pass
        # #     # plt.plot(np.arange(1024), ABBAt)
        #     plt.figure()
        #     plt.plot(XX, ABBAt[0:300], 'bo')
        #     plt.plot(XX, resultQ1.best_fit)
                
    print('ABBA Set ' + str(xx) + ' completed!')

Q1_IntErr = []
Q3_IntErr = []

for o in range(len(INTSQ1)):
    Q1IntErr = INTSQ1[o]*np.sqrt((Err_A0_Q1[o]/A0_Q1[o])**2 + (Err_A2_Q1[o]/A2_Q1[o])**2)
    Q1_IntErr.append(Q1IntErr)
    Q3IntErr = INTSQ3[o]*np.sqrt((Err_A0_Q3[o]/A0_Q3[o])**2 + (Err_A2_Q3[o]/A2_Q3[o])**2)
    Q3_IntErr.append(Q3IntErr)
    
T_Int_Err = np.sqrt((np.nanmean(Q1_IntErr)**2) + (np.nanmean(Q3_IntErr)**2))
print('Total Q1 + Q3 Err is =' + str(T_Int_Err*10**6))   

INTSQ1_14_Set1 = np.transpose(np.reshape(INTSQ1, (6, 31)))
Errs_Int_Q1_Set1 = np.transpose(np.reshape(Q1_IntErr, (6, 31)))
INTSQ3_14_Set1 = np.transpose(np.reshape(INTSQ3, (6, 31)))
Errs_Int_Q3_Set1 = np.transpose(np.reshape(Q3_IntErr, (6, 31)))
Total_shift_Set1 = np.fliplr(np.reshape(A1_Q1, (6, 31)))
Err_shift_Set1 = np.fliplr(np.reshape(Err_A1_Q1, (6, 31)))

#%% Find out velocities here
Shiftings = []
Vels_rough = []

for xx in range(6):
    Total_shifts = -1*(Total_shift_Set1[xx]-np.nanmean(Total_shift_Set1[xx, 2:-2]))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.53006))*299792.458
    Vels_rough.append(Vels)

#%%
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)
for xx in range(6,7):
    ABBAi = (COADD_ABBA[xx*2] + COADD_ABBA[(xx*2)+1])/2
    plt.figure()
    plt.imshow(ABBAi, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    for i in range(28): #28
        ii = i + 30 #30
        ABBAt = (ABBAi[ii,:])
        A01 = np.nanmax(ABBAt[135:155])
        a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[135:155]))
        a1_pointO1 = a1_pointO1[0][0]
        resultQ1 = gmodel.fit(ABBAt[0:300], x=XX, a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
        pQ1 = SimpleNamespace(**resultQ1.best_values)
        eQ1 = np.sqrt(np.diag(resultQ1.covar))
        A0_Q1.append(pQ1.a0)
        A1_Q1.append(pQ1.a1)
        A2_Q1.append(pQ1.a2*6.119575999463667e-05)
        Err_A0_Q1.append(eQ1[0])
        Err_A1_Q1.append(eQ1[1]*6.119575999463667e-05)
        Err_A2_Q1.append(eQ1[2]*6.119575999463667e-05)
        FWHM = pQ1.a2*6.119575999463667e-05*2*np.sqrt(2*math.log(2))
        FWHMQ1.append(FWHM)
        INTSQ1.append(pQ1.a0*FWHM)
    ABBAi = (COADD_ABBAQ3[xx*2] + COADD_ABBAQ3[(xx*2)+1])/2    
    for i in range(28): # 28
        ii = i + 30 #30
        ABBAt = (ABBAi[ii,:])
        A03 = np.nanmax(ABBAt[690:705])
        a1_pointO3 = np.where(ABBAt == np.nanmax(ABBAt[690:705]))
        a1_pointO3 = a1_pointO3[0][0] - 680
        resultQ3 = gmodel.fit(ABBAt[680:720], x=np.arange(40), a0=A03, a1=a1_pointO3, a2=1.3, a3=0, a4=0, a5=0)
        pQ3 = SimpleNamespace(**resultQ3.best_values)
        eQ3 = np.sqrt(np.diag(resultQ3.covar))
        A0_Q3.append(pQ3.a0)
        A1_Q3.append(pQ3.a1)
        A2_Q3.append(pQ3.a2*5.705793605517151e-05)
        Err_A0_Q3.append(eQ3[0])
        Err_A1_Q3.append(eQ3[1]*5.705793605517151e-05) #If you look for the gradient on wave at Q1 and Q3 it differs!
        Err_A2_Q3.append(eQ3[2]*5.705793605517151e-05)
        FWHM = pQ3.a2*5.705793605517151e-05*2*np.sqrt(2*math.log(2))
        FWHMQ3.append(FWHM)
        INTSQ3.append(pQ3.a0*FWHM)
        # if xx == 0:
        # # #     pass
        # #     # plt.plot(np.arange(1024), ABBAt)
        #     plt.figure()
        #     plt.plot(XX, ABBAt[0:300], 'bo')
        #     plt.plot(XX, resultQ1.best_fit)
                
    print('ABBA Set ' + str(xx) + ' completed!')

Q1_IntErr = []
Q3_IntErr = []

for o in range(len(INTSQ1)):
    Q1IntErr = INTSQ1[o]*np.sqrt((Err_A0_Q1[o]/A0_Q1[o])**2 + (Err_A2_Q1[o]/A2_Q1[o])**2)
    Q1_IntErr.append(Q1IntErr)
    Q3IntErr = INTSQ3[o]*np.sqrt((Err_A0_Q3[o]/A0_Q3[o])**2 + (Err_A2_Q3[o]/A2_Q3[o])**2)
    Q3_IntErr.append(Q3IntErr)
    
T_Int_Err = np.sqrt((np.nanmean(Q1_IntErr)**2) + (np.nanmean(Q3_IntErr)**2))
print('Total Q1 + Q3 Err is =' + str(T_Int_Err*10**6))   

INTSQ1_14_Set2 = np.transpose(np.reshape(INTSQ1, (1, 28)))
Errs_Int_Q1_Set2 = np.transpose(np.reshape(Q1_IntErr, (1, 28)))
INTSQ3_14_Set2 = np.transpose(np.reshape(INTSQ3, (1, 28)))
Errs_Int_Q3_Set2 = np.transpose(np.reshape(Q3_IntErr, (1, 28)))
Total_shift_Set2 = np.fliplr(np.reshape(A1_Q1, (1, 28)))
Err_shift_Set2 = np.fliplr(np.reshape(Err_A1_Q1, (1, 28)))

#%%
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)
for xx in range(7,9):
    ABBAi = (COADD_ABBA[xx*2] + COADD_ABBA[(xx*2)+1])/2
    # plt.figure()
    # plt.imshow(ABBAi, cmap='gist_gray')
    # if xx == 1:
    #     ABBAi[42,146] = np.nanmean(ABBAi) - 1e-05
    #     # ABBAi[45:52,149:152] = np.nanmean(ABBAi) - 1e-05
    #     ABBAi[29:32,150:152] = np.nanmean(ABBAi) - 1e-05
    #     # ABBAi[32:36,146] = np.nanmean(ABBAi) - 1e-05
    #     plt.figure()
    #     plt.imshow(ABBAi, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
    for i in range(31): #28
        ii = i + 28 #30
        ABBAt = (ABBAi[ii,:])
        A01 = np.nanmax(ABBAt[135:155])
        a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[135:155]))
        a1_pointO1 = a1_pointO1[0][0]
        resultQ1 = gmodel.fit(ABBAt[0:300], x=XX, a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
        pQ1 = SimpleNamespace(**resultQ1.best_values)
        eQ1 = np.sqrt(np.diag(resultQ1.covar))
        A0_Q1.append(pQ1.a0)
        A1_Q1.append(pQ1.a1)
        A2_Q1.append(pQ1.a2*6.119575999463667e-05)
        Err_A0_Q1.append(eQ1[0])
        Err_A1_Q1.append(eQ1[1]*6.119575999463667e-05)
        Err_A2_Q1.append(eQ1[2]*6.119575999463667e-05)
        FWHM = pQ1.a2*6.119575999463667e-05*2*np.sqrt(2*math.log(2))
        FWHMQ1.append(FWHM)
        INTSQ1.append(pQ1.a0*FWHM)
    ABBAi = (COADD_ABBAQ3[xx*2] + COADD_ABBAQ3[(xx*2)+1])/2    
    for i in range(28): # 28
        ii = i + 29 #30
        ABBAt = (ABBAi[ii,:])
        A03 = np.nanmax(ABBAt[690:705])
        a1_pointO3 = np.where(ABBAt == np.nanmax(ABBAt[690:705]))
        a1_pointO3 = a1_pointO3[0][0] - 680
        resultQ3 = gmodel.fit(ABBAt[680:720], x=np.arange(40), a0=A03, a1=a1_pointO3, a2=1.3, a3=0, a4=0, a5=0)
        pQ3 = SimpleNamespace(**resultQ3.best_values)
        eQ3 = np.sqrt(np.diag(resultQ3.covar))
        A0_Q3.append(pQ3.a0)
        A1_Q3.append(pQ3.a1)
        A2_Q3.append(pQ3.a2*5.705793605517151e-05)
        Err_A0_Q3.append(eQ3[0])
        Err_A1_Q3.append(eQ3[1]*5.705793605517151e-05) #If you look for the gradient on wave at Q1 and Q3 it differs!
        Err_A2_Q3.append(eQ3[2]*5.705793605517151e-05)
        FWHM = pQ3.a2*5.705793605517151e-05*2*np.sqrt(2*math.log(2))
        FWHMQ3.append(FWHM)
        INTSQ3.append(pQ3.a0*FWHM)
        # if xx == 0:
        # # #     pass
        # #     # plt.plot(np.arange(1024), ABBAt)
        #     plt.figure()
        #     plt.plot(XX, ABBAt[0:300], 'bo')
        #     plt.plot(XX, resultQ1.best_fit)
                
    print('ABBA Set ' + str(xx) + ' completed!')

Q1_IntErr = []
Q3_IntErr = []

for o in range(len(INTSQ1)):
    Q1IntErr = INTSQ1[o]*np.sqrt((Err_A0_Q1[o]/A0_Q1[o])**2 + (Err_A2_Q1[o]/A2_Q1[o])**2)
    Q1_IntErr.append(Q1IntErr)

for o in range(len(INTSQ3)):
    Q3IntErr = INTSQ3[o]*np.sqrt((Err_A0_Q3[o]/A0_Q3[o])**2 + (Err_A2_Q3[o]/A2_Q3[o])**2)
    Q3_IntErr.append(Q3IntErr)
    
T_Int_Err = np.sqrt((np.nanmean(Q1_IntErr)**2) + (np.nanmean(Q3_IntErr)**2))
print('Total Q1 + Q3 Err is =' + str(T_Int_Err*10**6))   

INTSQ1_14_Set3 = np.transpose(np.reshape(INTSQ1, (2, 31)))
Errs_Int_Q1_Set3 = np.transpose(np.reshape(Q1_IntErr, (2, 31)))
INTSQ3_14_Set3 = np.transpose(np.reshape(INTSQ3, (2, 28)))
Errs_Int_Q3_Set3 = np.transpose(np.reshape(Q3_IntErr, (2, 28)))
Total_shift_Set3 = np.fliplr(np.reshape(A1_Q1, (2, 31)))
Err_shift_Set3 = np.fliplr(np.reshape(Err_A1_Q1, (2, 31)))

#%%
# Take on all lines together
#Quickly find the Q3,1 line
Q3_1_Pos = []
Q3_0_Pos = []
Q1_0_Pos = []
Q2_0_Pos = []
Q3_2_Pos = []

COADD_ABBA_TotalQ3i = f.zoom(COADD_ABBA_TotalQ3, zoom = (1, np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[705:745]))))
#print(np.shape(COADD_ABBA_TotalQ3i))
COADD_ABBA_TotalQ31i = f.zoom(COADD_ABBA_TotalQ3, zoom = (1, np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[675:715]))))
#print(np.shape(COADD_ABBA_TotalQ31i))
COADD_ABBA_TotalQ2i = f.zoom(COADD_ABBA_TotalQ2, zoom = (1, np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[427:467]))))
#print(np.shape(COADD_ABBA_TotalQ2i))
COADD_ABBA_TotalQ32i = f.zoom(COADD_ABBA_TotalQ32, zoom = (1, np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[836:876]))))
#print(np.shape(COADD_ABBA_TotalQ32i))

for o in range(25):
    data = COADD_ABBA_TotalQ31i[40+o,755:795]
    A03 = np.nanmax(data)
    a1_pointO3 = np.where(data == np.nanmax(data))
    a1_pointO3 = a1_pointO3[0][0]
    popt, pcov = curve_fit(gaussian_fit, data, np.arange(40), p0=(A03, a1_pointO3, 1.5, 0, 0, 0))
    # plt.figure()
    # plt.plot(np.arange(40), data)
    # plt.plot(np.arange(40), gaussian_fit(np.arange(40), *popt))
    # print(popt[1])
    Q3_1_Pos.append(popt[1]+755)

for o in range(23):
    data = COADD_ABBA_TotalQ3i[40+o,730:770]
    A03 = np.nanmax(data)
    a1_pointO3 = np.where(data == np.nanmax(data))
    a1_pointO3 = a1_pointO3[0][0]
    popt, pcov = curve_fit(gaussian_fit, data, np.arange(40), p0=(A03, a1_pointO3, 1.5, 0, 0, 0))
    # plt.figure()
    # plt.plot(np.arange(40), data)
    # plt.plot(np.arange(40), gaussian_fit(np.arange(40), *popt))
    Q3_0_Pos.append(popt[1]+730)

for o in range(25):
    data = COADD_ABBA_Total[31+o,125:165]
    A03 = np.nanmax(data)
    a1_pointO3 = np.where(data == np.nanmax(data))
    a1_pointO3 = a1_pointO3[0][0]
    popt, pcov = curve_fit(gaussian_fit, data, np.arange(40), p0=(A03, a1_pointO3, 1.5, 0, 0, 0))
    # plt.figure()
    # plt.plot(np.arange(40), data)
    # plt.plot(np.arange(40), gaussian_fit(np.arange(40), *popt))
    Q1_0_Pos.append(popt[1]+125)

for o in range(25):
    data = COADD_ABBA_TotalQ2i[32+o,445:485]
    A03 = np.nanmax(data)
    a1_pointO3 = np.where(data == np.nanmax(data))
    a1_pointO3 = a1_pointO3[0][0]
    popt, pcov = curve_fit(gaussian_fit, data, np.arange(40), p0=(A03, a1_pointO3, 1.5, 0, 0, 0))
    # plt.figure()
    # plt.plot(np.arange(40), data)
    # plt.plot(np.arange(40), gaussian_fit(np.arange(40), *popt))
    Q2_0_Pos.append(popt[1]+445)
    
for o in range(25):
    data = COADD_ABBA_TotalQ32i[37+o,917:957]
    A03 = np.nanmax(data)
    a1_pointO3 = np.where(data == np.nanmax(data))
    a1_pointO3 = a1_pointO3[0][0]
    popt, pcov = curve_fit(gaussian_fit, data, np.arange(40), p0=(A03, a1_pointO3, 1.5, 0, 0, 0))
    # plt.figure()
    # plt.plot(np.arange(40), data)
    # plt.plot(np.arange(40), gaussian_fit(np.arange(40), *popt))
    Q3_2_Pos.append(popt[1]+917)

#%% So now make a super line

shifted_21 = (np.nanmean(Q1_0_Pos)) - (np.nanmean(Q2_0_Pos))
#Q2_0 = f.zoom(f.shift(COADD_ABBA_TotalQ2, shift=(0, shifted_21)), zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[427:467]))))
Q2_0 = f.shift(f.zoom(COADD_ABBA_TotalQ2, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[427:467])))), shift = (0, shifted_21))

shifted_31 = (np.nanmean(Q1_0_Pos[0:23])) - (np.nanmean(Q3_0_Pos))
#Q3_0 = f.zoom(f.shift(COADD_ABBA_TotalQ3, shift=(0, shifted_31)), zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[675:715]))))
Q3_0 = f.shift(f.zoom(COADD_ABBA_TotalQ3, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[675:715])))), shift = (0, shifted_31+2))

shifted_311 = (np.nanmean(Q1_0_Pos)) - (np.nanmean(Q3_1_Pos))
Q3_1 = f.shift(f.zoom(COADD_ABBA_TotalQ3, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[705:745])))), shift = (0, shifted_311-3))

shifted_321 = (np.nanmean(Q1_0_Pos)) - (np.nanmean(Q3_2_Pos))
Q3_2 = f.shift(f.zoom(COADD_ABBA_TotalQ32, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[836:876])))), shift = (0, shifted_321))

plt.figure()
plt.imshow((COADD_ABBA_Total+Q2_0[0:75,0:1024]+Q3_0[0:75,0:1024]+Q3_1[0:75,0:1024]+Q3_2[0:75,0:1024])/5, cmap='gist_gray', vmax = 0.01, vmin = -0.01)

Combo_Qlines = [] #Almost

#%%
for xx in range(9):
    ABBA1 = (COADD_ABBA[xx*2] + COADD_ABBA[(xx*2)+1])/2
    ABBA2 = f.shift(f.zoom((COADD_ABBAQ2[xx*2] + COADD_ABBAQ2[(xx*2)+1])/2, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[427:467])))), shift = (0, shifted_21))
    ABBA3 = f.shift((COADD_ABBAQ3[xx*2] + COADD_ABBAQ3[(xx*2)+1])/2, shift = (0, shifted_31+52.5))
    ABBA31 = f.shift(f.zoom((COADD_ABBAQ3[xx*2] + COADD_ABBAQ3[(xx*2)+1])/2, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[705:745])))), shift = (0, shifted_311-3))
    ABBA32 = f.shift(f.zoom((COADD_ABBAQ32[xx*2] + COADD_ABBAQ32[(xx*2)+1])/2, zoom = (np.nanmean(np.gradient(wave[125:165]))/np.nanmean(np.gradient(wave[836:876])))), shift = (0, shifted_321)) 
    Combo_Qlines.append((ABBA1+ABBA3[0:75,0:1024])/2)
    plt.figure()
    plt.imshow(((ABBA1+ABBA3[0:75,0:1024])/2), cmap='gist_gray', vmax = 0.005, vmin = -0.005)

#%% So now we add the lines together before doing A1 analysis
# Now here we see how A1 has been moving around
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)

for xx in range(4):
    ABBAi = (Combo_Qlines[xx*2] + Combo_Qlines[1+(xx*2)])/2
    # ABBAi = acre(ABBAi, width = 2, verbose = False)
    for i in range(33): #28
        if xx == 0:
            ii = i + 29 #30
        elif xx == 1:
            ii = i + 28
        elif xx == 2:
            ii = i + 28
        else:
            ii = i + 29
        ABBAt = (ABBAi[ii,:])
        A01 = np.nanmax(ABBAt[140:155])
        a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[140:155]))
        a1_pointO1 = a1_pointO1[0][0] - 130
        #print(a1_pointO1, i)
        # print(a1_pointO1)
        resultQ1 = gmodel.fit(ABBAt[130:165], x=np.arange(35), a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
        pQ1 = SimpleNamespace(**resultQ1.best_values)
        eQ1 = np.sqrt(np.diag(resultQ1.covar))
        A0_Q1.append(pQ1.a0)
        A1_Q1.append(pQ1.a1)
        A2_Q1.append(pQ1.a2*np.nanmean(np.gradient(wave[147:149])))
        Err_A0_Q1.append(eQ1[0])
        #Err_A1_Q1.append(eQ1[1])
        Err_A1_Q1.append(eQ1[1]*np.nanmean(np.gradient(wave[147:149])))
        #print(pQ1.a1, eQ1[1], i)
        Err_A2_Q1.append(eQ1[2]*np.nanmean(np.gradient(wave[147:149])))
        FWHM = pQ1.a2*np.nanmean(np.gradient(wave[147:149]))*2*np.sqrt(2*math.log(2))
        FWHMQ1.append(FWHM)
        INTSQ1.append(pQ1.a0*FWHM) 
        # plt.figure()
        # plt.plot(np.arange(35), ABBAt[130:165])
        # plt.plot(np.arange(35), resultQ1.best_fit)
    # plt.figure()
    # plt.imshow(ABBAi[:, 130:165], cmap='gist_gray')
    # plt.plot(A1_Q1, np.arange(28, 61), color='b')


#%%
Total_lines = np.fliplr(np.reshape(INTSQ1, (4, 33)))
Total_shift = np.fliplr(np.reshape(A1_Q1, (4, 33)))

Vels_rough = []
Shiftings = []

for xx in range(4):
    Total_shifts = -1*(Total_shift[xx]-np.nanmean(Total_shift))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.965))*299792.458
    Vels_rough.append(Vels)

Q1_VelErr = []
Shiftings = np.reshape(Shiftings, (1, 132))

Errors_Vels = np.fliplr(np.reshape(Err_A1_Q1, (4, 33)))
Errors_Vels = np.reshape(Errors_Vels, (1, 132))

Vels = np.reshape(Vels_rough, (1, 132))

for o in range(len(INTSQ1)):
    Q1VelErr = (np.sqrt((Errors_Vels[:,o]/(Shiftings[:,o]))**2))
    Q1_VelErr.append(np.sqrt((Q1VelErr*Vels[:,o])**2))

Q1_VelErr = np.fliplr(np.reshape(Q1_VelErr, (4, 33)))

Longitudes = np.load('Longitudes5.npy')
Latitudes = np.load('Latitudes5.npy')

Longs = []
for o in range(33):
    Longs.append((Longitudes[o]+Longitudes[o+34]+Longitudes[o+1]+Longitudes[o+35])/4)

Lats = []
for o in range(33):
    Lats.append((Latitudes[o] + Latitudes[o+1] + Latitudes[o+34] + Latitudes[o+35])/4)

#%%
R_Uranus = 28559*1000 
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
Period_Uranus = []

for o in range(33):
    Period_Uranus2 = 2*np.pi*R_Uranus*np.cos(np.deg2rad(Lats[o]))/ (17.24*60*60)
    Period_Uranus.append(Period_Uranus2)
    Limb_velocity_pix_O4 = (Period_Uranus2/1000)
    Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

#So the speed of Uranus at surface is 
# Look into how to vary this across the slit so it can be used for 2016
Limb_velocity_pix_O4_test = []
Limb_velocity_pix_O4_test2 = []

Planet_diameter = []
for o in range(33):
    Planet_diameter.append((16.380257074056157*np.cos(np.deg2rad(Lats[o]))))

for a in range(17):
    Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((16-a)*np.pi/(2*Planet_diameter[a])))
for a in range(17):
    Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((16-a)*np.pi/(2*Planet_diameter[a])))
    
Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))

TIMES = ['10:45', '11:07', '11:29', '11:51']
for a in range(4):
    plt.figure()
    plt.plot(Longs, Vels_rough[a], color='b')
    plt.fill_between(Longs, y1 = Vels_rough[a]-Q1_VelErr[a], y2 = Q1_VelErr[a]+Vels_rough[a], alpha=0.5, color='b')
    plt.plot(Longs, Planet_rotation_O4[:], ls = '--', lw = 3, color = 'k')    
    plt.ylim(-4, 4)
    plt.tick_params(axis='both', which='major', labelsize=20)
    plt.tick_params(axis='both', which='major', labelsize=20)
    # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
    plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
    plt.ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
    plt.xlim(-90, 90) #0, 3601
    plt.grid(alpha = 0.7, linestyle = '--')
    plt.title('~20 min exposure around ' + str(TIMES[a]) + ' UTC with Keck/NIRSPEC on $13^{th}$ October 2014', fontsize=22, pad=25)

#%% Now need to average the Keck 

Total_Q_Lines = (Combo_Qlines[0] + Combo_Qlines[1] + Combo_Qlines[2] + Combo_Qlines[3] + Combo_Qlines[4] + Combo_Qlines[5] + Combo_Qlines[6] + Combo_Qlines[7])/8
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)

for xx in range(33):
    ABBAi = Total_Q_Lines
    # ABBAi = acre(ABBAi, width = 2, verbose = False)
    ABBAt = (ABBAi[28+xx,:])
    A01 = np.nanmax(ABBAt[140:155])
    a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[140:155]))
    a1_pointO1 = a1_pointO1[0][0] - 130
    #print(a1_pointO1, i)
    # print(a1_pointO1)
    resultQ1 = gmodel.fit(ABBAt[130:165], x=np.arange(35), a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
    pQ1 = SimpleNamespace(**resultQ1.best_values)
    eQ1 = np.sqrt(np.diag(resultQ1.covar))
    A0_Q1.append(pQ1.a0)
    A1_Q1.append(pQ1.a1)
    A2_Q1.append(pQ1.a2*np.nanmean(np.gradient(wave[147:149])))
    Err_A0_Q1.append(eQ1[0])
    #Err_A1_Q1.append(eQ1[1])
    Err_A1_Q1.append(eQ1[1]*np.nanmean(np.gradient(wave[147:149])))
    #print(pQ1.a1, eQ1[1], i)
    Err_A2_Q1.append(eQ1[2]*np.nanmean(np.gradient(wave[147:149])))
    FWHM = pQ1.a2*np.nanmean(np.gradient(wave[147:149]))*2*np.sqrt(2*math.log(2))
    FWHMQ1.append(FWHM)
    INTSQ1.append(pQ1.a0*FWHM)

Q1_IntErr_P = []
Q1_IntErr_N = []
Q1_IntErr = []

for o in range(33):
    QIntErr = INTSQ1[o]*np.sqrt((Err_A0_Q1[o]/A0_Q1[o])**2 + (Err_A2_Q1[o]/A2_Q1[o])**2)
    Q1_IntErr.append(QIntErr)
    Q1_IntErr_P.append(INTSQ1[o]+QIntErr)
    Q1_IntErr_N.append(INTSQ1[o]-QIntErr)

Vels_rough = []
Shiftings = []

for xx in range(33):
    Total_shifts = (A1_Q1[xx]-np.nanmean(A1_Q1))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.965))*299792.458
    Vels_rough.append(Vels)

Q1_VelErr = []

for o in range(len(INTSQ1)):
    Q1VelErr = (np.sqrt((Err_A1_Q1[o]/(Shiftings[o]))**2))
    Q1_VelErr.append(np.sqrt((Q1VelErr*Vels_rough[o])**2))

Period_Uranus = (2*np.pi*R_Uranus)/ (17.24*60*60)

Planet_rotation_O4_V1 = (np.linspace(16, -16, 33)*np.sin(np.deg2rad(90-27.35))*Period_Uranus)/(16.380257074056157*1000)

plt.figure()
plt.plot(Longs, Vels_rough, color='b')
plt.fill_between(Longs, y1 = np.array(Vels_rough[0:33])-np.array(Q1_VelErr[0:33]), y2 = np.array(Q1_VelErr[0:33])+np.array(Vels_rough[0:33]), alpha=0.5, color='b')
plt.plot(Longs, Planet_rotation_O4[:], ls = '--', lw = 3, color = 'k', label = 'Manual calculation Velocity') 
plt.plot(Longs, Planet_rotation_O4_V1, ls = '-.', lw = 3, color = 'k', alpha = 0.75, label= 'Sub-Earth-Latitude Velocity')   
plt.ylim(-4, 4)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
plt.xlim(-90, 90) #0, 3601
plt.grid(alpha = 0.7, linestyle = '--')
plt.title('~1hr exposure around 11:20 UTC with Keck/NIRSPEC on $13^{th}$ October 2014', fontsize=22, pad=25)
plt.legend(loc = 'upper center', fontsize = 25)
plt.vlines(72.61138080642526, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
plt.vlines(-72.61138080642526, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '11:20' +' UTC on $13^{th}$ October 2014', fontsize=35, pad=25)
ax.plot(Longs, Planet_rotation_O4_V1, color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Longs, Vels_rough, color='b', ls= '--', label='Combined $Q(1,0)$ and $Q(3,0)$ IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Longs, np.array(Vels_rough[0:33])-np.array(Q1_VelErr[0:33]), np.array(Vels_rough[0:33])+np.array(Q1_VelErr[0:33]), color='b', alpha=0.5)
ax2.plot(Longs, np.flip(INTSQ1)*10**6, color='r', label='IR Intensity', lw = 5)
#ax2.fill_between(Longs, np.flip(Q1_IntErr_N)*10**6, np.flip(Q1_IntErr_P)*10**6, color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=35, labelpad=15)
ax.set_ylabel('(ORF) $H_{3}^{+}$ Velocities (kms$^{-1}$)', fontsize=35, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=35)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 1)
ax.set_ylim(-4, 4)
ax.legend(loc='lower center', fontsize=35)
ax2.legend(loc='upper center', fontsize=35)
ax.grid(linestyle = '--')
#plt.vlines(46.76483059080472, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-46.76483059080472, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)

import scipy 

def lin_func(x, a , b):
    y = a* x + b
    return y

# Velocities
VelS_popt, VelS_pcov = scipy.optimize.curve_fit(lin_func, np.arange(2, 33, 1), Vels_rough[1:-1], p0 = [-0.3, 0], sigma = Q1_VelErr[1:-1])
Pos_parms = VelS_popt + np.sqrt(np.diag(VelS_pcov))
Neg_parms = VelS_popt - np.sqrt(np.diag(VelS_pcov))

Linefit = lin_func(np.arange(2, 33), *VelS_popt)
Linefit2 = lin_func(np.arange(2, 33), *Pos_parms)
Linefit3 = lin_func(np.arange(2, 33), *Neg_parms)

Min_Super_rotation = (1 + (Pos_parms[0]/((2.40775607*2)/31)))*-100 #(%)
Max_Super_rotation = (1 + (Neg_parms[0]/((2.40775607*2)/31)))*-100 #(%)

Sub_Earth_Velocity = np.linspace(2.5682731361673013, -2.5682731361673013, 33)
      
fig, ax = plt.subplots(figsize=(10,8))
ax.plot(np.arange(2, 33, 1), Vels_rough[1:-1], 'ko')
ax.plot(np.arange(2, 33, 1), Linefit, color = 'b', label = 'Best fit')
ax.plot(np.arange(2, 33, 1), Linefit2, color = 'r', label = 'Best fit + error')
ax.plot(np.arange(2, 33, 1), Linefit3, color = 'g', label = 'Best fit - error')
ax.plot(np.arange(2, 33, 1), Sub_Earth_Velocity[1:-1], 'k--', label = 'Planetary Rotation')
ax.tick_params(axis='both', which='major', labelsize=20)
#ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude (Pixel No.)', fontsize=25, labelpad=15)
ax.set_ylabel('(PRF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
#ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(1, 33) #0, 3601
ax.set_ylim(-4, 4)
#ax2.set_ylim(0, 0.75)
ax.set_title('~1 hr exposure around ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)

#%% Lets h3ppy over a total and break down from there
h3p = h3ppy.h3p(line_list_file = 'C:/Users/snowy/OneDrive/Documents/Python work/h3p_line_list_neale_1996_subset.txt')
wave = np.load('Wavelength_Order.npy')
Q1Ints = []
Q3Ints = []
Q1Max = []
Q3Max = []
Temps = []
T_err = []
Col_Densities = []
Col_err = []
T_emission = []
Tot_EM_Err = []
Shift = []
Shift_err = []
c = 0

def subdivide(wave, spec, middles, width = 20) : 
    ret = []
    for m in middles : 
        centre = np.abs(wave - m).argmin()
        for i in range(centre - width, centre + width) : 
            ret.append(spec[i])
    return np.array(ret)

# The H3+ line centeres contained withing this spectral band
centers = [3.953]
centers3 = [3.971, 3.986, 3.9945]
cpos = np.arange(1) * 41 + 20
subspecs = []
subspecsQ2 = []
subspecsQ3 = []
subspecsQ31 = []
subspecsQ32 = []

for xx in range(1):
    ABBAi = (COADD_ABBA[0] + COADD_ABBA[1] + COADD_ABBA[2] + COADD_ABBA[3] + COADD_ABBA[4] + COADD_ABBA[5])/6
    ABBAf = (COADD_ABBAQ3[0] + COADD_ABBAQ3[1] + COADD_ABBAQ3[2] + COADD_ABBAQ3[3] + COADD_ABBAQ3[4] + COADD_ABBAQ3[5])/6
    for i in range(32):
        subspec = subdivide(wave, ABBAi[28+i,:], centers)
        subspecQ3 = subdivide(wave, ABBAf[28+i,:], centers3)
        subwave = subdivide(wave, wave, centers)
        subwaveQ3 = subdivide(wave, wave, centers3)
        subspec = np.concatenate((subspec, subspecQ3))
        subwave = np.concatenate((subwave, subwaveQ3))
        h3p.set(wavelength = subwave, data = subspec, R = 16400, temperature = 550, density = 6*10**15)
        # Guess the density
        h3p.guess_density(verbose=False)
        # By generating a model at this point we can compare our initial guess to the data
        guess = h3p.model()
        # First let's just fit the background
        h3p.set(nbackground = 1)
        ptf = ['background_0']
        fit = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results(verbose=False)
        # Only fit a linear fit as the emission background is level and a quadratic for offset  
        # Remember, the variables are stored in the h3p object so it'll remember the results 
        # of the last fit
        h3p.set(noffset = 2)
        ptf = ['offset_0', 'offset_1']
        fit2 = h3p.fit(params_to_fit = ptf, verbose = False)
        Vars, errs = h3p.get_results(verbose=False)
        # Then let's fit the temperature, density and line-width
        ptf = ['temperature', 'density', 'sigma_0']
        fit3 = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results(verbose = False)
        if np.nanmean(fit3) == 0:
            h3p.set(nbackground = 1, wavelength = subwave, data = subspec, R = 20000, temperature = 600, density = 6*10**15)
            h3p.guess_density(verbose=False)
            guess = h3p.model()
            h3p.set(nbackground = 1)
            ptf = ['background_0']
            fit = h3p.fit(params_to_fit = ptf)
            Vars, errs = h3p.get_results(verbose=False)
            # Only fit a linear fit as the emission background is level and a quadratic for offset  
            # Remember, the variables are stored in the h3p object so it'll remember the results 
            # of the last fit
            h3p.set(noffset = 2)
            ptf = ['offset_0', 'offset_1']
            fit2 = h3p.fit(params_to_fit = ptf, verbose = False)
            Vars, errs = h3p.get_results(verbose=False)
            # Then let's fit the temperature, density and line-width
            ptf = ['temperature', 'density', 'sigma_0']
            fit3 = h3p.fit(params_to_fit = ptf)
            Vars, errs = h3p.get_results(verbose = False)
            print('Redone previous fit, fit now converges Pixel No. ' + str(i))
        Temps.append(Vars['temperature'])
        T_err.append(errs['temperature'])
        Col_Densities.append(Vars['density'])
        Col_err.append(errs['density']) 
        Shift.append(Vars['offset_0'])
        # Shift.append(Vars['offset_1'])
        Shift_err.append(errs['offset_0'])
        # fig, ax = plt.subplots()
        # ax.plot(np.arange(160), subspec, label='Raw Data')
        # ax.plot(np.arange(160), fit3, label = '$H_3^+$ spectrum fit')
        # ax.plot(np.arange(160), guess, label = 'Guessed model')

#%% Sort all and seeing edges!!!!!
plt.figure()
plt.plot(Longs[1:], np.flip(Temps), color='r')
plt.fill_between(Longs[1:], np.flip(Temps+np.array(T_err)), np.flip(Temps-np.array(T_err)), color='r', alpha = 0.5)
plt.ylim(350, 750)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('Temperatures of $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
plt.xlim(-90, 90) #0, 3601
plt.grid(alpha = 0.7, linestyle = '--')
plt.title('~1hr exposure around 11:20 UTC with Keck/NIRSPEC on $13^{th}$ October 2014', fontsize=22, pad=25)

plt.figure()
plt.plot(Longs[1:], np.flip(Col_Densities)/1e16, color='b')
plt.fill_between(Longs[1:], np.flip(Col_Densities+np.array(Col_err))/1e16, np.flip(Col_Densities-np.array(Col_err))/1e16, color='b', alpha = 0.5)
plt.ylim(0, 1.6)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
plt.xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
plt.ylabel('Column Densities of $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
plt.xlim(-90, 90) #0, 3601
plt.grid(alpha = 0.7, linestyle = '--')
plt.title('~1hr exposure around 11:20 UTC with Keck/NIRSPEC on $13^{th}$ October 2014', fontsize=22, pad=25)

#%%
fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '11:20' +' UTC on $13^{th}$ October 2014', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Longs[1:], np.flip(Col_Densities)/1e16, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Longs[1:], np.flip(Col_Densities+np.array(Col_err))/1e16, np.flip(Col_Densities-np.array(Col_err))/1e16, color='purple', alpha=0.5)
ax2.plot(Longs, np.flip(INTSQ1)*10**6, color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Longs, np.flip(Q1_IntErr_N)*10**6, np.flip(Q1_IntErr_P)*10**6, color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax2.set_ylabel('Ionospheric Column Density of $H_{3}^{+}$ $(x 10^{17}m^{-2})$', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(linestyle = '--')

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '11:20' +' UTC on $13^{th}$ October 2014', fontsize=35, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Longs[1:], np.flip(Temps), color='g', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Longs[1:], np.flip(Temps+np.array(T_err)), np.flip(Temps-np.array(T_err)), color='g', alpha=0.5)
ax2.plot(Longs[1:], np.flip(Col_Densities)/1e16, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
ax2.fill_between(Longs[1:], np.flip(Col_Densities+np.array(Col_err))/1e16, np.flip(Col_Densities-np.array(Col_err))/1e16, color='purple', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=35, labelpad=15)
ax.set_ylabel('$H_{3}^{+}$ Temperature $(K)$', fontsize=35, labelpad=15)
ax2.set_ylabel('$H_{3}^{+}$ Column Density $(x 10^{16}m^{-2})$', fontsize=35, labelpad=15)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 10)
ax.set_ylim(350, 850)
ax.legend(loc='upper left', fontsize=35)
ax2.legend(loc='upper right', fontsize=35)
ax.grid(linestyle = '--')

#%%
ABBAi = []
A0_Q1 = []
Err_A0_Q1 = []
A0_Q3 = []
Err_A0_Q3 = []
A1_Q1 = []
Err_A1_Q1 = []
A1_Q3 = []
Err_A1_Q3 = []
A2_Q1 = []
Err_A2_Q1 = []
A2_Q3 = []
Err_A2_Q3 = []
FWHMQ1 = []
FWHMQ3 = []
INTSQ1 = []
INTSQ3 = []
#from AvParams import AVALSAv
iiii = 0
gmodel = Model(gaussian_fit)

for a in range(8):
    ABBAi = Combo_Qlines[a]
    for xx in range(30):
        if a == 0:
            ABBAt = (ABBAi[28+xx,:])
        else:
            ABBAt = (ABBAi[29+xx,:])
        A01 = np.nanmax(ABBAt[140:155])
        a1_pointO1 = np.where(ABBAt == np.nanmax(ABBAt[140:155]))
        a1_pointO1 = a1_pointO1[0][0] - 130
        #print(a1_pointO1, i)
        # print(a1_pointO1)
        resultQ1 = gmodel.fit(ABBAt[130:165], x=np.arange(35), a0=A01, a1=a1_pointO1, a2=1.3, a3=0, a4=0, a5=0)
        pQ1 = SimpleNamespace(**resultQ1.best_values)
        eQ1 = np.sqrt(np.diag(resultQ1.covar))
        A0_Q1.append(pQ1.a0)
        A1_Q1.append(pQ1.a1)
        A2_Q1.append(pQ1.a2*np.nanmean(np.gradient(wave[147:149])))
        Err_A0_Q1.append(eQ1[0])
        #Err_A1_Q1.append(eQ1[1])
        Err_A1_Q1.append(eQ1[1]*np.nanmean(np.gradient(wave[147:149])))
        #print(pQ1.a1, eQ1[1], i)
        Err_A2_Q1.append(eQ1[2]*np.nanmean(np.gradient(wave[147:149])))
        FWHM = pQ1.a2*np.nanmean(np.gradient(wave[147:149]))*2*np.sqrt(2*math.log(2))
        FWHMQ1.append(FWHM)
        INTSQ1.append(pQ1.a0*FWHM)
    if a == 0:
        plt.figure()
        plt.imshow(ABBAi[28:28+30,130:165], cmap='gist_gray', vmax = np.nanmean(ABBAi) + 6*np.nanstd(ABBAi), vmin = np.nanmean(ABBAi) - 6*np.nanstd(ABBAi))
        plt.plot(A1_Q1[a*30:(a*30)+30], np.arange(30))
    else:
        plt.figure()
        plt.imshow(ABBAi[29:29+30,130:165], cmap='gist_gray', vmax = np.nanmean(ABBAi) + 6*np.nanstd(ABBAi), vmin = np.nanmean(ABBAi) - 6*np.nanstd(ABBAi))
        plt.plot(A1_Q1[a*30:(a*30)+30], np.arange(30))

#%%
Q1_IntErr_P = []
Q1_IntErr_N = []
Q1_IntErr = []

for o in range(len(INTSQ1)):
    QIntErr = INTSQ1[o]*np.sqrt((Err_A0_Q1[o]/A0_Q1[o])**2 + (Err_A2_Q1[o]/A2_Q1[o])**2)
    Q1_IntErr.append(QIntErr)
    Q1_IntErr_P.append(INTSQ1[o]+QIntErr)
    Q1_IntErr_N.append(INTSQ1[o]-QIntErr)
    
Vels_rough = []
Shiftings = []
a = 0

for xx in range(len(INTSQ1)):
    Total_shifts = (A1_Q1[xx]-np.nanmean(A1_Q1[a*30:(a*30)+30]))*np.nanmean(np.gradient(wave[145:150]))
    Shiftings.append(Total_shifts)
    Vels = (Total_shifts/(3.965))*299792.458
    Vels_rough.append(Vels)
    if (xx + 1) % 33 == 0:
        a += 1

Q1_VelErr = []

for o in range(len(INTSQ1)):
    Q1VelErr = (np.sqrt((Err_A1_Q1[o]/(Shiftings[o]))**2))
    Q1_VelErr.append(np.sqrt((Q1VelErr*Vels_rough[o])**2))

Period_Uranus = (2*np.pi*R_Uranus)/ (17.24*60*60)

Planet_rotation_O4_V1 = (np.linspace(14, -15, 30)*np.sin(np.deg2rad(90-27.35))*Period_Uranus)/(16.380257074056157*1000)
    
TIMES = ['10:38', '10:47', '10:55', '11:03', '11:11', '11:15', '11:24', '11:32']

for a in range(8):
    fig, ax = plt.subplots(figsize=(10,8))
    ax2 = ax.twinx()
    ax.set_title('~10 minute exposure at ' + str(TIMES[a]) +' UTC on $13^{th}$ October 2014', fontsize=22, pad=25)
    ax.plot(Longs[2:-1], Planet_rotation_O4_V1, color='k', ls = '--', label='Planetary Rotation')
    #ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
    #ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
    ax.plot(Longs[2:-1], Vels_rough[a*30:(a*30)+30], color='b', ls= '--', label='Combined $Q(1,0)$ and $Q(3,0)$ IR Velocities', lw = 5)
    #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
    #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
    ax.fill_between(Longs[2:-1], np.array(Vels_rough[a*30:(a*30)+30])-np.array(Q1_VelErr[a*30:(a*30)+30]), np.array(Vels_rough[a*30:(a*30)+30])+np.array(Q1_VelErr[a*30:(a*30)+30]), color='b', alpha=0.5)
    ax2.plot(Longs[2:-1], np.flip(INTSQ1[a*30:(a*30)+30])*10**6, color='r', label='IR Intensity', lw = 5)
    ax2.fill_between(Longs[2:-1], np.flip(Q1_IntErr_N[a*30:(a*30)+30])*10**6, np.flip(Q1_IntErr_P[a*30:(a*30)+30])*10**6, color='r', alpha=0.5)
    #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
    #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
    # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
    ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax2.tick_params(axis='both', which='major', labelsize=20)
    # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
    ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
    ax.set_ylabel('Ionospheric Column Density of $H_{3}^{+}$ $(x 10^{17}m^{-2})$', fontsize=25, labelpad=15)
    ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
    #ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    #ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
    ax.set_xlim(-90, 90) #0, 3601
    ax2.set_ylim(0, 0.75)
    ax.set_ylim(-4, 4)
    ax.legend(loc='lower left', fontsize=25)
    ax2.legend(loc='upper right', fontsize=25)
    ax.grid(linestyle = '--')