# -*- coding: utf-8 -*-
"""
Created on Mon Apr  6 12:16:02 2020

@author: Emma 'Amelia'
"""
# Here we will extract out the star files, the flat files and dark files to be used in the first stage of data reduction

import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits #import the relevant directories to read a fits file from your directory and plot it
from types import SimpleNamespace
from skimage import exposure
from matplotlib import cm
from gaussstar import gaussstar
from scipy.optimize import curve_fit

from lmfit import Model
Keck_Data = [] #First create lists in which the arrays of Keck data will go into
Keck_DataABBA = []
s = 0 #This will be used to create the ABBA pattern

#Find sky lines and match these
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
#%%
for n in range(72):
    if n == 0:
        Keck_Data_Total = Keck_Data[n]
    else:
        Keck_Data_Total += Keck_Data[n]

plt.figure()
plt.imshow(Keck_Data_Total, cmap='gist_gray')

#%%
data =np.loadtxt("mktrans_zm_30_10.dat") #import the data

#First lets fit the Q1, Q2, Q3, Q3,1 and Q3,2 if possible (Only Q3,1 possible)
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Waves = data.T[0]
Trans = data.T[1]

S_wave_O4 = np.where(Waves == 3.94400)
F_wave_O4 = np.where(Waves == 3.99700)
S_waveO4 = S_wave_O4[0][0]
F_waveO4 = F_wave_O4[0][0]
Waves_O4 = Waves[S_waveO4:F_waveO4]
Trans_O4 = Trans[S_waveO4:F_waveO4]

plt.figure()
plt.plot(Waves_O4, Trans_O4)

#%% Now do a calibration comparsion of both
# Start = np.nanmean(Keck_Data_Total[50:150,0:901], axis = 0)/(30*10000)
# p_lower, p_higher = np.percentile(Start, (0.01, 99.9))
# Linecut_Arc = exposure.rescale_intensity(Start, in_range=(p_lower, p_higher))

# Wave_O4 = np.load('Wavelength_Order.npy')
# fig, ax1 = plt.subplots()

# ax1.plot(Wave_O4[0:901], Start, color='k')
# ax1.set_ylim(2, 6)
# ax1.tick_params(axis='y')
# ax1.set_xlabel(r'Wavelength ($\mu$m)', fontsize=20)
# ax1.set_ylabel('Detector counts (counts$^{-1}$ x 10$^{4}$)', fontsize=20)
# labels = ['3.944', '3.950', '3.956', '3.962', '3.968', '3.974', '3.980', '3.986', '3.991', '3.997']
# ax1.set_xticks([Wave_O4[0], Wave_O4[100], Wave_O4[200], Wave_O4[300], Wave_O4[400], Wave_O4[500], Wave_O4[600], Wave_O4[700], Wave_O4[800], Wave_O4[900]], labels=labels, fontsize = 15)
# labels = ['2', '3', '4', '5', '6']
# ax1.set_yticks([2, 3, 4, 5, 6], labels=labels, fontsize = 15)
# ax1.grid()

# ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

# ax2.set_ylabel('Transmission', fontsize = 20)  # we already handled the x-label with ax1
# ax2.plot(Waves_O4, Trans_O4, color = 'gray', linestyle = '--')
# labels = ['0', '0.2', '0.4', '0.6', '0.8', '1.0']
# ax2.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1.0], labels=labels, fontsize = 15)

# plt.plot(Linecut_Arc , color='k')
# plt.xlabel(r'Wavelength ($\mu$m)', fontsize=20)
# plt.ylabel('Normalised Detector counts', fontsize=20)
# labels = ['3.944', '3.950', '3.956', '3.962', '3.968', '3.974', '3.980', '3.986', '3.991', '3.997']
# plt.xticks([0,100, 200, 300, 400, 500, 600, 700, 800, 900], labels=labels, fontsize = 15)
# plt.yticks(fontsize = 15)

#%%
#Lets check the four lines then
A1_L1 = []
A1_Err_L1 = []
gmodel = Model(gauss_fit)
for i in range(140):
    Sky_No_O4 = Keck_Data_Total[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[15:50]))
    a1_pointP = a1_pointP[0][0] - 10
    result_L1 = gmodel.fit(Sky_No_O4[10:60], x=np.arange(50), a0=np.nanmax(Sky_No_O4[15:50]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+10)
    A1_Err_L1.append(eL1[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation
def lin_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

A1_L1[1] = 32.3
A1_L1[2] = 32.3
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)

popt1, pcov1 = curve_fit(lin_func, xdata, ydata, p0=[9.07220639e-06, -1.93469812e-03, 0.126, 29.6])
# plt.figure(figsize = (12,10))
# plt.plot(np.arange(140), ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt1), 'g--')
# plt.ylabel(r'Wavelength Pixel Position (Pixel No.)', fontsize=25)
# plt.xlabel('Spatial Pixel Position across the Slit (Pixel No.)', fontsize=25)
# plt.xticks(np.arange(0, 141, 10), labels = np.arange(25, 166, 10), fontsize = 20)
# plt.yticks(fontsize = 20)
# plt.grid(linestyle='--', alpha = 0.5)

#%% Second line
A1_L1 = []
A1_Err_L1 = []
gmodel = Model(gauss_fit)
for i in range(140):
    Sky_No_O4 = Keck_Data_Total[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[105:120]))
    a1_pointP = a1_pointP[0][0] - 85
    result_L1 = gmodel.fit(Sky_No_O4[85:135], x=np.arange(50), a0=np.nanmax(Sky_No_O4[105:120]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+85)
    A1_Err_L1.append(eL1[1])
    
# A1_L1[1] = 32.3
# A1_L1[2] = 32.3
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)

popt2, pcov2 = curve_fit(lin_func, xdata, ydata, p0=[9.07220639e-06, -1.93469812e-03, 0.126, 29.6])
# plt.figure()
# plt.plot(np.arange(140), ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt2), 'g--')

#%% line 3
A1_L1 = []
A1_Err_L1 = []
gmodel = Model(gauss_fit)
for i in range(140):
    Sky_No_O4 = Keck_Data_Total[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[210:230]))
    a1_pointP = a1_pointP[0][0] - 195
    result_L1 = gmodel.fit(Sky_No_O4[195:245], x=np.arange(50), a0=np.nanmax(Sky_No_O4[210:230]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+195)
    A1_Err_L1.append(eL1[1])
    
A1_L1[5] = 223
A1_L1[6] = 223
A1_L1[8] = 223
A1_L1[11] = 223
A1_L1[18] = 222.9
A1_L1[19] = 222.9
A1_L1[20] = 222.9
A1_L1[26] = 222.8
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)

popt3, pcov3 = curve_fit(lin_func, xdata, ydata, p0=[9.07220639e-06, -1.93469812e-03, 0.126, 29.6])
# plt.figure()
# plt.plot(np.arange(140), ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt3), 'g--')

#%% line 4
A1_L1 = []
A1_Err_L1 = []
gmodel = Model(gauss_fit)
for i in range(140):
    Sky_No_O4 = Keck_Data_Total[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[330:345]))
    a1_pointP = a1_pointP[0][0] - 315
    result_L1 = gmodel.fit(Sky_No_O4[315:365], x=np.arange(50), a0=np.nanmax(Sky_No_O4[330:345]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+315)
    A1_Err_L1.append(eL1[1])
    
# A1_L1[5] = 223
# A1_L1[6] = 223
# A1_L1[8] = 223
# A1_L1[11] = 223
# A1_L1[18] = 222.9
# A1_L1[19] = 222.9
# A1_L1[20] = 222.9
# A1_L1[26] = 222.8
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)

popt4, pcov4 = curve_fit(lin_func, xdata, ydata, p0=[9.07220639e-06, -1.93469812e-03, 0.126, 29.6])
# plt.figure()
# plt.plot(np.arange(140), ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt4), 'g--')

#%% Line 5
A1_L1 = []
A1_Err_L1 = []
gmodel = Model(gauss_fit)
for i in range(140):
    Sky_No_O4 = Keck_Data_Total[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[625:635]))
    a1_pointP = a1_pointP[0][0] - 610
    result_L1 = gmodel.fit(Sky_No_O4[610:660], x=np.arange(50), a0=np.nanmax(Sky_No_O4[625:635]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+610)
    A1_Err_L1.append(eL1[1])
    
# A1_L1[5] = 223
# A1_L1[6] = 223
# A1_L1[8] = 223
# A1_L1[11] = 223
# A1_L1[18] = 222.9
# A1_L1[19] = 222.9
# A1_L1[20] = 222.9
# A1_L1[26] = 222.8
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)

popt5, pcov5 = curve_fit(lin_func, xdata, ydata, p0=[9.07220639e-06, -1.93469812e-03, 0.126, 29.6])
# plt.figure()
# plt.plot(np.arange(140), ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt5), 'g--')

#%%
import scipy.ndimage as f #This will allow the equivalent of fshift-ing in IDL

#Set up the shift factor(look to automate this)
#corrective_s_factor = (A_average + B_average)
A_average = popt1[0]
B_average = popt1[1]
C_average = popt1[2]
D_average = popt1[3]
Corrective_Array = []

for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 172:
        ii = i - 10
        Sky_row_O4 = Keck_Data_Total[i,:]
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average) - 32.26538491060073)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Keck_Data_Total, cmap='gist_gray')

A1_L1 = []
A1_Err_L1 = []
for i in range(140):
    Sky_No_O4 = shift_Sky_O4[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[15:50]))
    a1_pointP = a1_pointP[0][0] - 10
    result_L1 = gmodel.fit(Sky_No_O4[10:60], x=np.arange(50), a0=np.nanmax(Sky_No_O4[15:50]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+10)
    A1_Err_L1.append(eL1[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation

A1_L1[1] = 32.5
A1_L1[2] = 32.5
ydata = A1_L1
xdata = np.arange(140)
XPOS_L1 = np.nanmean(ydata)

popt1A, pcov1A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 31.917])
# plt.figure()
# plt.plot(np.arange(140), A1_L1, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt1A), 'g--')

#%%
Corrective_Array2 = []
A_average = popt1A[0]
B_average = popt1A[1]
C_average = popt1A[2]
D_average = popt1A[3]

for i in range(len(shift_Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 172:
        ii = i - 10
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average) - XPOS_L1)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)
         
#%%
A_average = popt2[0]
B_average = popt2[1]
C_average = popt2[2]
D_average = popt2[3]

for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = Keck_Data_Total[i,:]
        ii = i - 10
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average) - 112.11)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Keck_Data_Total, cmap='gist_gray')

A1_L1 = []
A1_Err_L1 = []
for i in range(140):
    Sky_No_O4 = shift_Sky_O4[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[105:120]))
    a1_pointP = a1_pointP[0][0] - 85
    result_L1 = gmodel.fit(Sky_No_O4[85:135], x=np.arange(50), a0=np.nanmax(Sky_No_O4[105:120]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+85)
    A1_Err_L1.append(eL1[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation

ydata = A1_L1
xdata = np.arange(140)
XPOS_L2 = np.nanmean(ydata)

popt2A, pcov2A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 100])
# plt.figure()
# plt.plot(np.arange(140), A1_L1, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt2A), 'g--')

#%%
A_average = popt2A[0]
B_average = popt2A[1]
C_average = popt2A[2]
D_average = popt2A[3]

for i in range(len(shift_Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 172:
        ii = i - 10
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L2)*-1)
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)
         
#%%
A_average = popt3[0]
B_average = popt3[1]
C_average = popt3[2]
D_average = popt3[3]

for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = Keck_Data_Total[i,:]
        ii = i - 10
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average) - 222.82)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Keck_Data_Total, cmap='gist_gray')

A1_L1 = []
A1_Err_L1 = []
for i in range(140):
    Sky_No_O4 = shift_Sky_O4[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[210:230]))
    a1_pointP = a1_pointP[0][0] - 195
    result_L1 = gmodel.fit(Sky_No_O4[195:245], x=np.arange(50), a0=np.nanmax(Sky_No_O4[210:230]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+195)
    A1_Err_L1.append(eL1[1])
    
A1_L1[5] = 223
A1_L1[6] = 223
A1_L1[8] = 223
A1_L1[11] = 223
A1_L1[18] = 222.9
A1_L1[19] = 222.9
A1_L1[20] = 222.9
A1_L1[26] = 222.8
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)
XPOS_L3 = np.nanmean(ydata)

popt3A, pcov3A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 31.917])
# plt.figure()
# plt.plot(np.arange(140), A1_L1, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt3A), 'g--')

#%%
A_average = popt3A[0]
B_average = popt3A[1]
C_average = popt3A[2]
D_average = popt3A[3]

for i in range(len(shift_Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = (((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average) + D_average - XPOS_L3)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)
#%%
A_average = popt4[0]
B_average = popt4[1]
C_average = popt4[2]
D_average = popt4[3]

for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = Keck_Data_Total[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 337.14)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Keck_Data_Total, cmap='gist_gray')

A1_L1 = []
A1_Err_L1 = []
for i in range(140):
    Sky_No_O4 = shift_Sky_O4[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[330:345]))
    a1_pointP = a1_pointP[0][0] - 315
    result_L1 = gmodel.fit(Sky_No_O4[315:365], x=np.arange(50), a0=np.nanmax(Sky_No_O4[330:345]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+315)
    A1_Err_L1.append(eL1[1])
    
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)
XPOS_L4 = np.nanmean(ydata)

popt4A, pcov4A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 31.917])
# plt.figure()
# plt.plot(np.arange(140), A1_L1, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt4A), 'g--')

#%%
A_average = popt4A[0]
B_average = popt4A[1]
C_average = popt4A[2]
D_average = popt4A[3]

for i in range(len(shift_Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L4)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)
    
#%%
A_average = popt5[0]
B_average = popt5[1]
C_average = popt5[2]
D_average = popt5[3]

for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = Keck_Data_Total[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 627.93)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Keck_Data_Total, cmap='gist_gray')

A1_L1 = []
A1_Err_L1 = []
for i in range(140):
    Sky_No_O4 = shift_Sky_O4[25+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[615:635]))
    a1_pointP = a1_pointP[0][0] - 605
    result_L1 = gmodel.fit(Sky_No_O4[605:655], x=np.arange(50), a0=np.nanmax(Sky_No_O4[615:635]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+605)
    A1_Err_L1.append(eL1[1])
    
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(140)
XPOS_L5 = np.nanmean(ydata)

popt5A, pcov5A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 31.917])
# plt.figure()
# plt.plot(np.arange(140), A1_L1, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt5A), 'g--')

#%%
A_average = popt5A[0]
B_average = popt5A[1]
C_average = popt5A[2]
D_average = popt5A[3]

for i in range(len(shift_Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 172:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L5)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)
#%%
Corrections = np.array(Corrective_Array)
Corrections2 = np.array(Corrective_Array2)
CorrectionsA = Corrections.reshape((5, 183))
CorrectionsB = Corrections2.reshape((5, 183))
X_Positions = np.asarray([np.round(XPOS_L1), np.round(XPOS_L2), np.round(XPOS_L3), np.round(XPOS_L4), np.round(XPOS_L5)])

def lin_func(x, a , b, c):
    y = a* x**2 + b*x + c
    return y

Corrections_Spectra = []
gmodel = Model(lin_func)
Margin2 = []

for m in range(183):
    poptT, pcovT = curve_fit(lin_func, X_Positions, CorrectionsA[:,m], p0=[0, 0, 2])
    Corrections_Spectra.append(lin_func(np.arange(1024), *poptT))
    # print(lin_func(2013, best_val.a, best_val.b, best_val.c))
    #Margin2.append(lin_func(2013, best_val.a, best_val.b, best_val.c, best_val.d))
    
plt.figure()
plt.imshow(Corrections_Spectra)
plt.colorbar()

Corrections_SpectraB = []

for m in range(183):
    poptTB, pcovTB = curve_fit(lin_func, X_Positions, CorrectionsB[:,m], p0=[0, 0, 0.5])
    Corrections_SpectraB.append(lin_func(np.arange(1024), *poptTB))
    
plt.figure()
plt.imshow(Corrections_SpectraB)
plt.colorbar()

Corrections_SpectraT = np.add(Corrections_SpectraB, Corrections_Spectra)

#%% Now we use the Correction Spectra to fix the spectra just for the Q1 spectra
for i in range(len(Keck_Data_Total[:,0])):
    if i == 0:
        Sky_row_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = Sky_row_O4
    elif i < 10:
        shifted_Sky_O4 = Keck_Data_Total[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
    elif i > 10 and i < 172:
        Sky_row_O4 = Keck_Data_Total[i,:]
        corrective_shift_factor = np.nanmean(Corrections_SpectraT[i][135:155])
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
    else:
         shifted_Sky_O4 = Keck_Data_Total[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))

np.save('Corrective_Array.npy', Corrections_SpectraT) #Corrrect the data and refind the wavelength range for the equations 
plt.figure()
plt.imshow(shift_Sky_O4, cmap='magma')
plt.figure()
plt.imshow(Keck_Data_Total, cmap='magma')

#%%
# Top = np.zeros((35,1024))
# Bottom = np.zeros((27, 1024))

# Data = np.concatenate((Top, shift_Sky_O4))
# Data = np.concatenate((Data, Bottom))

# plt.figure()
# plt.imshow(Data[:,0:1001], cmap='magma', vmax = np.nanmean(Data)+2*np.nanstd(Data), vmin = 1)
# plt.xlabel(r'Wavelength Axis (Pixel)', fontsize=25)
# #label_x = 0, 200, 400, 600, 800, 1000
# #plt.xticks(label_x, ('3.9445', '3.9565', '3.9685', '3.9804', '3.9924', '4.0044'), fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel('Spatial Axis (Pixel)', fontsize=25)
# #plt.hlines(60, xmin=0, xmax=149, color='w', linestyle='--')
# # cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
# # cbar.ax.tick_params(labelsize=15)
# # cbar.set_label('Spectral Pixel Adjustment (Pixel)', fontsize=25)

#%% Now we make a figure to show the corrective array
# Top = np.zeros((35,1024))
# Bottom = np.zeros((27, 1024))

# New_Fig = np.concatenate((Top, Corrections_SpectraT[0:173, :]))
# New_Fig = np.concatenate((New_Fig, Bottom))

# plt.figure(figsize=(10, 7.5))
# plt.imshow(New_Fig)
# #plt.title(r'First ABBA observation set of Uranus with NIRSPEC on $5^{th}$ September 2006', fontsize=23)
# plt.xlabel(r'Wavelength Axis (5.7851 x 10$^{-5}$ $\mu$m per pixel)', fontsize=25)
# #label_x = 0, 200, 400, 600, 800, 1000
# #plt.xticks(label_x, ('3.9445', '3.9565', '3.9685', '3.9804', '3.9924', '4.0044'), fontsize=15)
# plt.xticks(fontsize=15)
# plt.yticks(fontsize=15)
# plt.ylabel('Spatial Axis (0.144" per pixel)', fontsize=25)
# #plt.hlines(60, xmin=0, xmax=149, color='w', linestyle='--')
# cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
# cbar.ax.tick_params(labelsize=15)
# cbar.set_label('Spectral Pixel Adjustment (Pixel)', fontsize=25)

#%%

def lin_func(x, a , b, c):
    y = a* x**2 + b*x + c
    return y

Corrections_Spectra = []

#Skyline 1 Rough 32.5 = 3.9459797 #NEED TO DO WAVELNGTH CALS
#Skyline 2 Rough 241.3 = 3.9587773
#Skyline 3 Rough 250.7 = 3.9593299
#Skyline 4 Rough 628.1 = 3.9815994

Ratio1 = (3.9587773-3.945979)/(241.3-32.5)
Ratio2 = (3.9593299-3.9587773)/(250.7-241.3)
Ratio3 = (3.9815994-3.9593299)/(628.1-250.7)

plt.figure()
plt.plot([32.5, 241.3, 250.7, 628.1], [3.9459797, 3.9587773, 3.9593299, 3.9815994], 'bo')
popt, pcov = curve_fit(lin_func, [32.5, 241.3, 250.7, 628.1], [3.9459797, 3.9587773, 3.9593299, 3.9815994], p0=[2, 1, 20])
plt.plot(np.arange(1024), lin_func(np.arange(1024), *popt), 'g--')

np.save('Wavelength_Order.npy', lin_func(np.arange(1024), *popt))
# Wav_to_Pix = (Ratio1 + Ratio2 + Ratio3)/3

# Start_Wav_O1 = 3.9459797 - (Wav_to_Pix*32.5)
# End_Wav_O1 = Start_Wav_O1 + 1024*Wav_to_Pix

# print(Start_Wav_O1)
# print(End_Wav_O1)

#%%
Star1_Data = []
Star2_Data = []
Darks_Data = []
Flats_Data = []

for n in range(12): #We use this list to create a list which holds all the data from Order19
    num = n + 117
    if num < 121:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/order1/order1s0' + str(num) + '.fits.gz'
        image_datai = fits.getdata(image_filei, ext=0)
        Darks_Data.append(image_datai)       
    elif num >= 121 and num < 125:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/order1/order1s0' + str(num) + '.fits.gz'
        fits.open(image_filei)
        image_datai = fits.getdata(image_filei, ext=0)
        Flats_Data.append(image_datai)
    else:
        image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/order1/order1s0' + str(num) + '.fits.gz'
        hdul = fits.open(image_filei)
        hdr = hdul[0].header
        print(hdr['OBJECT'])
        image_datai = fits.getdata(image_filei, ext=0)
        Star1_Data.append(image_datai)

Gauss_Star_Factor = gaussstar()
s = 0
Darks = (Darks_Data[0] + Darks_Data[1] + Darks_Data[2] + Darks_Data[3])/4

True_Darks = Darks

Calibrate1 = (Flats_Data[s] - True_Darks)
Calibrate2 = (Flats_Data[s+1] - True_Darks)
Calibrate3 = (Flats_Data[s+2] - True_Darks)
Calibrate4 = (Flats_Data[s+3] - True_Darks)

Cal_Flats = (Calibrate1 + Calibrate2 + Calibrate3 + Calibrate4)/(4*10)

Cal_Darks = True_Darks

#Sort the view of the imshow
p_lower, p_higher = np.percentile(Cal_Darks, (0.01, 99.9))
array_to_view = exposure.rescale_intensity(Cal_Darks, in_range=(p_lower, p_higher))

fig = plt.figure(figsize=(12,8))
plt.imshow(array_to_view, cmap='gist_gray')
plt.xlabel(r'Wavelength Pixel Position (5.7851 x 10$^{-5}$ $\mu$m per pixel)', fontsize=20)
plt.ylabel('Spatial Pixel Position across the Slit (Pixel No.)', fontsize=20)
plt.title('Dark Frame', fontsize=25)
m = cm.ScalarMappable(cmap='gist_gray')
m.set_array(array_to_view)
cbar = fig.colorbar(m, ticks=[-1, -0.5, 0, 0.5, 1])
cbar.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1.0'], fontsize=12.5)
cbar.set_label('Normalisted Detector counts', fontsize=20)

p_lower, p_higher = np.percentile(Cal_Flats, (5, 95))
array_to_view = exposure.rescale_intensity(Cal_Flats, in_range=(p_lower, p_higher))

fig = plt.figure(figsize=(12,8))
plt.imshow(array_to_view, cmap='gist_gray')
plt.xlabel(r'Wavelength Pixel Position (5.7851 x 10$^{-5}$ $\mu$m per pixel)', fontsize=20)
plt.ylabel('Spatial Pixel Position across the Slit (Pixel No.)', fontsize=20)
plt.title('Flat Frame', fontsize=25)
m = cm.ScalarMappable(cmap='gist_gray')
m.set_array(array_to_view)
cbar = fig.colorbar(m, ticks=[0, 0.25, 0.5, 0.75, 1])
cbar.ax.set_yticklabels(['0', '0.25', '0.5', '0.75', '1.0'], fontsize=12.5)
cbar.set_label('Normalised Detector counts', fontsize=20)

Star1_DataA1 = (Star1_Data[0]/(10*Gauss_Star_Factor)) - Cal_Flats
Star1_DataB1 = (Star1_Data[1]/(10*Gauss_Star_Factor)) - Cal_Flats
Star1_DataB2 = (Star1_Data[2]/(10*Gauss_Star_Factor)) - Cal_Flats
Star1_DataA2 = (Star1_Data[3]/(10*Gauss_Star_Factor)) - Cal_Flats

Star1 = (Star1_DataA1 - Star1_DataB1 - Star1_DataB2 + Star1_DataA2)/40

plt.figure()
plt.imshow(Star1, cmap='gist_gray')
plt.axhline(84, color='r', ls='--')
plt.axhline(160, color='b', ls='--')

# plt.figure()
# plt.imshow(Star1, cmap='gist_gray')
# plt.title('Order 19 spectra of HD 215 143 observed on 5th September 2006', fontsize=25)
# plt.xlabel(r'Wavelength $\lambda$ (+- 0.0016 $\mu$m)', fontsize=20)
# label_x = 0, 200, 400, 600, 800, 1000
# plt.xticks(label_x, ('3.9448', '3.9564', '3.9680', '3.9795', '3.9911', '4.0003'), fontsize=10)
# plt.ylabel('Spatial Pixel No. across the slit (Pixel No.)', fontsize=20)
# plt.axvline(696.652257000841, color='r', ls='--')
# plt.axvline(142.22591552273104, color='r', ls='--')
# plt.axhline(82, color='r', ls='--')
# plt.axhline(158, color='b', ls='--')

#plt.figure()
#plt.imshow(Cal_Lamps, cmap='gist_gray')
#plt.xlabel('Wavelength Pixel Position', fontsize=20)
#plt.ylabel('Spatial Pixel Position across the Slit', fontsize=20)
#plt.title('Scaled Lamp Arcs across Order 19', fontsize=25)
#cbar = plt.colorbar() #Prints out an image in greyscale of the fits file
#cbar.set_label('Uncalibrated Intensity (CCD Counts per second)', fontsize=20)

#%%
#Now we will determine the correct fit for the line by fitting a gaussian over the star fits files and finding the peak position 
# Now match a gaussian profile across all 5 lines to find the intensity, and width of the slit
def gaussian_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Star1P = Star1[25:75,:]
Star1N = Star1[100:150,:]

# Star2P = Star1[50:120,:]
# Star2N = Star2[125:195,:]
# X = np.arange(70)

Star1PPeak = []
Star1NPeak = []
A0Star1P = []
A2Star1P = []
A0Star1N = []
A2Star1N = []

for i in range(1024):
    StarP = Star1P[:,i]
    StarN = Star1N[:,i]
    gmodel = Model(gaussian_fit)
    result1P = gmodel.fit(StarP, x=np.arange(50), a0=600, a1=22.5, a2=2.0, a3=0, a4=0, a5=0)
    result1N = gmodel.fit(StarN, x=np.arange(50), a0=-850, a1=31.5, a2=2.0, a3=0, a4=0, a5=0)
    p = SimpleNamespace(**result1P.best_values)
    n = SimpleNamespace(**result1N.best_values)
    Star1PPeak.append(p.a1+25)
    Star1NPeak.append(n.a1+100)
    A0Star1P.append(p.a0)
    A2Star1P.append(p.a2)
    A0Star1N.append(n.a0)
    A2Star1N.append(n.a2)
    if i < 510 and i > 500:
        plt.figure()
        plt.plot(StarP)
        plt.plot(StarN)

Peak1P = np.array(Star1PPeak)
Peak1N = np.array(Star1NPeak)
ABBASep = Peak1N - Peak1P
XI = np.arange(1024)

plt.figure()
plt.title(r'Peak pixel position for both A frame and B frame emissions of HD 218 639')
plt.xlabel(r'Wavelength $\lambda$ (+- 0.0016 $\mu$m)')
label_x = 0, 200, 400, 600, 800, 1000
plt.xticks(label_x, ('3.9448', '3.9564', '3.9680', '3.9795', '3.9911', '4.0003'), fontsize=10)
plt.ylabel(r'Spatial Pixel No. across the slit (Pixel No.)')
plt.plot(XI, Peak1P, color='r')
plt.plot(XI, Peak1N, color='b')
plt.plot(XI, ABBASep, color='g')
plt.ylim(0,150)
plt.xlim(0,1024)
#%%

#Now find the difference between the AR and DEC
PyThag1 = np.sqrt((37.03979-37.04174)**2 + (8.46013-8.46285)**2)
Sep_Ratio = PyThag1*3600
print('Arc seconds seperation for Star 1 is ' + str(np.mean(Sep_Ratio)) + ' +- ' + str(np.std(Sep_Ratio)) + ' Arc seconds')

Spat_Ratio = (PyThag1*3600)/ABBASep
print('The spatial scale for Star 1 is ' + str(np.mean((Spat_Ratio))) + ' +- ' + str(np.std((Spat_Ratio))) + ' Arc second per pixel')

# #Now we will do the exact same for Star 2 HD 215 143 as we did for Star 1 HD 218 639
# Star2PPeak = []
# Star2NPeak = []
# A0Star2P = []
# A2Star2P = []
# A0Star2N = []
# A2Star2N = []

# for i in range(1024):
#     StarP = Star2P[:,i]
#     StarN = Star2N[:,i]
#     gmodel = Model(gaussian_fit)
#     result2P = gmodel.fit(StarP, x=X, a0=54.09817216036432, a1=31.9, a2=2.6, a3=0, a4=0, a5=0)
#     result2N = gmodel.fit(StarN, x=X, a0=-54.09817216036432, a1=33.7, a2=2.65, a3=0, a4=0, a5=0)
#     p = SimpleNamespace(**result2P.best_values)
#     n = SimpleNamespace(**result2N.best_values)
#     Star2PPeak.append(p.a1+50)
#     Star2NPeak.append(n.a1+125)
#     A0Star2P.append(p.a0)
#     A2Star2P.append(p.a2)
#     A0Star2N.append(n.a0)
#     A2Star2N.append(n.a2)

# Peak2P = np.array(Star2PPeak)
# Peak2N = np.array(Star2NPeak)
# ABBA2Sep = Peak2N - Peak2P
# XI = np.arange(1024)

# plt.figure()
# plt.title(r'Peak pixel position for both A frame and B frame emissions of HD 215 143', fontsize=20)
# plt.xlabel(r'Wavelength $\lambda$ (+- 0.0016 $\mu$m)', fontsize=15)
# label_x = 0, 200, 400, 600, 800, 1000
# plt.xticks(label_x, ('3.9448', '3.9564', '3.9680', '3.9795', '3.9911', '4.0003'), fontsize=10)
# plt.ylabel(r'Spatial Pixel No. across the slit (Pixel No.)', fontsize=15)
# plt.plot(XI, Peak2P, color='r')
# plt.plot(XI, Peak2N, color='b')
# plt.plot(XI, ABBA2Sep, color='g')
# plt.ylim(0,228)
# plt.xlim(0,1024)

# #Now find the difference between the AR and DEC
# Pythag2 = np.sqrt((-6.963025+6.962605)**2 + (340.81135-340.807805)**2)
# Sep_Ratio2 = Pythag2*3600
# print('Arc seconds seperation for Star 2 is ' + str(np.mean(Sep_Ratio2)) + ' +- ' + str(np.std(Sep_Ratio2)) + ' Arc seconds')

# Spat_Ratio2 = (Pythag2*3600)/ABBA2Sep
# print('The spatial scale for Star 2 is ' + str(np.mean(Spat_Ratio2)) + ' +- ' + str(np.std(Spat_Ratio2)) + ' Arc second per pixel')

#%%
#Now that the Spatial Ratio's are known between Star's we will now work out their Observed Flux to compare it against the expected.
A0Star1P = np.array(A0Star1P)
A2Star1P = np.array(A2Star1P)
A0Star1N = np.array(A0Star1N)
A2Star1N = np.array(A2Star1N)

FWHM_Star1P = A2Star1P*2*np.sqrt(2*np.log(2))
FWHM_Star1N = A2Star1N*2*np.sqrt(2*np.log(2))

Flux_Star1 = (((A0Star1P*FWHM_Star1P)+(-1*A0Star1N*FWHM_Star1N))/2)

#Now we need to find the spectral scaling across each order to find the flux calibration

#From here we need to find the expected Intensity of each Star to find the ratio used for our Data files, as found in Rosie's Thesis
#This is between 3.94 and 4.00 micrometers so the atmospheric window is most likely in the L' filter
F_HR718 = 4.07*(10**-10)
F_Star = 7.3*(10**-11)
M_lambda_HR718 = 4.392 # at 2.2 micrometers
M_lambda_star = 3.4

F_A0_HR718 = F_HR718*(10**(-0.4*M_lambda_HR718))
F_A0_star = F_Star*(10**(-0.4*M_lambda_star))

hc = (6.63*10**-34)*(2.99*10**8)
kb = 1.38*10**-23
T_HR718 = 10630
T_star = 10000
Lambda_AW = 2.2 # In micrometers
Lambda_AW1 = 3.975
Lambda = lin_func(np.arange(1024), *popt)

Fbb_HR718 = F_A0_HR718*((Lambda_AW/Lambda)**5)*(np.exp((hc/kb)/(Lambda_AW*T_HR718))/(np.exp((hc/kb)/(Lambda*T_HR718))))
Fbb_star = F_A0_star*((Lambda_AW/Lambda)**5)*(np.exp((hc/kb)/(Lambda_AW*T_star))/(np.exp((hc/kb)/(Lambda*T_star))))

Flux_Star1 = (Flux_Star1)/Lambda

#Now we can compare the lines against the expected values for the flux of these stars
fig, ax = plt.subplots()
ax.plot(Lambda, Flux_Star1, color='r', label='HR 718 observed')
#ax.plot(Lambda, Flux_Star2, color='k', label='HD 215 143 observed')
ax.set_xlabel(r'Wavelength  $\lambda$ ($\mu$m)', fontsize=15)
ax.set_ylabel(r'Uncalibrated Intensity (CCD countss$^{-1}$$\mu$m$^{-1}$)', fontsize=15)
ax.set_title(r'Observed and Expected Flux of HD 215 143 across Order 19', fontsize=20)
ax.set_ylim(0, 100)
ax.set_xlim(Lambda[0], Lambda[1023])
ax.grid()
ax2=ax.twinx()
ax2.plot(Lambda, Fbb_HR718, color='k', label='HR 718 expected')
# ax2.plot(Lambda, Fbb_HD215, color='b', label='HD 215 143 expected')
# ax2.plot(Lambda, Fbb_HR1578, color='g', label='HR 1578 expected')
# ax2.plot(Lambda, Fbb_star, color='r', label='Henrik example')
#ax2.set_ylabel(r'Blackbody Flux (Wm$^{-2}$$\mu$m$^{-1}$)', fontsize=15)
#ax2.set_ylim(0, 40*(10**-13))
ax.legend(loc='upper left')
ax2.legend(loc='upper right')

#Finally we will work out the Intensity to CCD counts ratio for order 19 which can be used on the data from order 19
Fc_HR718 = 2*Fbb_HR718/Flux_Star1
# Fc_HD215 = Fbb_HD215/Flux_Star2

# plt.figure()
# plt.plot(Lambda, Fc_HD218, color='r', label='HD 218 639 calibration spectrum')
# plt.plot(Lambda, Fc_HD215, color='k', label='HD 215 143 calibration spectrum')
# plt.xlabel(r'Wavelength $\lambda$ ($\mu$m)', fontsize=15)
# plt.ylabel(r'Flux (Wm$^{-2}$$\mu$m$^{-1}$Counts$^{-1}$)', fontsize=15)
# plt.title(r'Flux Calibration Spectre HD 215 143 aross Order 19', fontsize=20)
# plt.legend(loc='upper right')
# plt.xlim(Lambda[0], Lambda[1023])
# plt.grid()
# plt.ylim(0.18*(10**-15), 0.3*(10**-15))