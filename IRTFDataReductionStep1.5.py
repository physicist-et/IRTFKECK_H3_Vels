# -*- coding: utf-8 -*-
"""
Created on Mon Jun 14 11:52:28 2021

@author: snowy
"""
#Quick Data reduction for 2016 Data to find spatial resolution
from astropy.io import fits
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import math
from scipy.optimize import curve_fit

#First lets fit the Q1, Q2, Q3, Q3,1 and Q3,2 if possible (Only Q3,1 possible)
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Sky_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/sky.order.4.fits.gz'
Sky_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/sky.order.5.fits.gz'
Wave_info3 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.3.fits.gz'
Wave_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.4.fits.gz'
Wave_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/wavelength.order.5.fits.gz'

Wavelength_O3 = fits.getdata(Wave_info3, ext=0)
Wavelength_O4 = fits.getdata(Wave_info4, ext=0)
Wavelength_O5 = fits.getdata(Wave_info5, ext=0)

Sky_O4 = fits.getdata(Sky_info4, ext=0)
Sky_O5 = fits.getdata(Sky_info5, ext=0)

plt.imshow(Sky_O4, cmap='gist_gray', vmax = 75000, vmin = -75000)
plt.figure()
plt.imshow(Sky_O5, cmap='gist_gray', vmax = 75000, vmin = -75000)

#%%First plot the emission and transmission lines for Earth's atmosphere

data =np.loadtxt("mktrans_zm_30_10.dat") #import the data

Waves = data.T[0]
Trans = data.T[1]

S_wave_O4 = np.where(Waves == 3.93290)
F_wave_O4 = np.where(Waves == 3.96550)
S_waveO4 = S_wave_O4[0][0]
F_waveO4 = F_wave_O4[0][0]
Waves_O4 = Waves[S_waveO4:F_waveO4]
Trans_O4 = Trans[S_waveO4:F_waveO4]

plt.figure()
plt.plot(Waves_O4, Trans_O4,)
plt.xlabel(r'Wavelength ($\mu$m)')
plt.ylabel(r'Transmission Value')
plt.ylim(0, 1.0) #show the plot of Wavelength vs Transmission
plt.xlim(Wavelength_O4[0], Wavelength_O4[-1])

m= 25
Tot_Sky_O4 = np.zeros(2048)

for i in range(50):
    T_Sky_O4 = Sky_O4[m+i,:]
    Tot_Sky_O4 = np.add(Tot_Sky_O4, T_Sky_O4)

plt.figure()
#plt.plot(Wavelength_O4, (Tot_Sky_O4)/50)
plt.plot(np.arange(2048), (Tot_Sky_O4)/50)
#plt.xlabel(r'Wavelength ($\mu$m)')
#plt.xlim(3.93290, 3.96550)

# Now for the second order

Waves = data.T[0]
Trans = data.T[1]
S_wave_O5 = np.where(Waves == 3.96310)
F_wave_O5 = np.where(Waves == 3.99620)
S_waveO5 = S_wave_O5[0][0]
F_waveO5 = F_wave_O5[0][0]
Waves_O5 = Waves[S_waveO5:F_waveO5]
Trans_O5 = Trans[S_waveO5:F_waveO5]

plt.figure()
plt.plot(Waves_O5, Trans_O5)
plt.xlabel(r'Wavelength ($\mu$m)')
plt.ylabel(r'Transmission Value')
plt.ylim(0, 1.0) #show the plot of Wavelength vs Transmission
plt.xlim(Wavelength_O5[0], Wavelength_O5[-1])

m= 25
Tot_Sky_O5 = np.zeros(2048)

for i in range(50):
    T_Sky_O5 = Sky_O5[m+i,:]
    Tot_Sky_O5 = np.add(Tot_Sky_O5, T_Sky_O5)

plt.figure()
plt.plot(Wavelength_O5, (Tot_Sky_O5)/50)
plt.xlabel(r'Wavelength ($\mu$m)')
plt.xlim(3.96310, 3.99620)

#%% Need to find the angle of the star lines on both orders
#Lets start with the easier order O4

#Pick out 3 lines one from the fron, middle and back
#Front Need to go through all rows and find the centre of the gaussian curve and see if its moved!

Emission_O4 = Trans_O4[874:924]

gmodel = Model(gauss_fit)
result_L1 = gmodel.fit(Emission_O4, x=Waves_O4[874:924], a0=-0.15, a1=3.9509, a2=0.00018, a3=0, a4=0, a5=0)
pL1 = SimpleNamespace(**result_L1.best_values)
Line_1 = pL1.a1

result_L1A = gmodel.fit(Tot_Sky_O4[1065:1095], x=np.arange(30), a0=4300000, a1=20, a2=2, a3=0, a4=0, a5=0)
pL1A = SimpleNamespace(**result_L1A.best_values)
Line_1A = pL1A.a1 + 1065

Emission_O4 = Trans_O4[1210:1260]

result_L5 = gmodel.fit(Emission_O4, x=Waves_O4[1210:1260], a0=np.nanmin(Emission_O4), a1=3.957618, a2=0.00162683, a3=0.89, a4=0, a5=0)
pL5 = SimpleNamespace(**result_L5.best_values)
Line_5 = pL5.a1

result_L5A = gmodel.fit(Tot_Sky_O4[1720:1760], x=np.arange(40), a0=55000, a1=18.7041288, a2=2.64, a3=0, a4=0, a5=0)
pL5A = SimpleNamespace(**result_L5A.best_values)
Line_5A = pL5A.a1 + 1720

Emission_O4 = Trans_O4[1125:1175]

result_L2 = gmodel.fit(Emission_O4, x=Waves_O4[1125:1175], a0=np.nanmin(Emission_O4), a1=3.95592, a2=0.00162683, a3=0, a4=0, a5=0)
pL2 = SimpleNamespace(**result_L2.best_values)
Line_2 = pL2.a1

result_L2A = gmodel.fit(Tot_Sky_O4[1390:1420], x=np.arange(30), a0=14000, a1=10, a2=2, a3=0, a4=0, a5=0)
pL2A = SimpleNamespace(**result_L2A.best_values)
Line_2A = pL2A.a1 + 1390


#%% Here's the first line

A1_L1 = []
A1_Err_L1 = []
for i in range(75):
    Sky_No_O4 = Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[24:30]))
    a1_pointP = a1_pointP[0][0] - 10
    result_L1 = gmodel.fit(Sky_No_O4[10:60], x=np.arange(50), a0=np.nanmax(Sky_No_O4[24:30]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+10)
    A1_Err_L1.append(eL1[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation
def lin_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

A1_L1[5] = 26.5
A1_L1[28] = 25.6
A1_L1[70] = 24.3
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(75)

popt1, pcov1 = curve_fit(lin_func, xdata, ydata, p0=[0, 0, -0.032, 26.5])
plt.figure()
plt.plot(np.arange(75), ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt1), 'g--')

#%%Here's the middle line 500

A1_L2 = []
A1_Err_L2 = []
for i in range(75):
    Sky_No_O4 = Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[665:685]))
    a1_pointP = a1_pointP[0][0] - 655
    result_L1 = gmodel.fit(Sky_No_O4[655:705], x=np.arange(50), a0=np.nanmax(Sky_No_O4[665:685]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L2.append(pL1.a1+655)
    A1_Err_L2.append(eL1[1])
    
ydata = A1_L2
xdata = np.arange(75)

popt2, pcov2 = curve_fit(lin_func, xdata, ydata, p0=[0, -0.00049, -0.028, 675])
plt.figure()
plt.plot(np.arange(75), A1_L2, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt2), 'g--')

#%%Here's the end line
A1_L4 = []
A1_Err_L4 = []
for i in range(75):
    if i == 24 or i == 74:
        A1_L4.append(pL1.a1+1170)
        A1_Err_L4.append(eL1[1])
    else:
        Sky_No_O4 = Sky_O4[10+i,:]
        a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1180:1200]))
        a1_pointP = a1_pointP[0][0] - 1170
        result_L1 = gmodel.fit(Sky_No_O4[1170:1220], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1180:1200]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        pL1 = SimpleNamespace(**result_L1.best_values)
        if result_L1.covar is None:
            eL1 = 0,0,0
            A1_L4.append(pL1.a1+1170)
            A1_Err_L4.append(eL1[1])            
        else:
            eL1 = np.sqrt(np.diag(result_L1.covar))
            A1_L4.append(pL1.a1+1170)
            A1_Err_L4.append(eL1[1])
    
A1_L4[2] = 1189.7
A1_L4[15] = 1188.9
A1_L4[53] = 1187
A1_L4[46] = 1187.6
ydata = A1_L4
xdata = np.arange(75)

popt4, pcov4 = curve_fit(lin_func, xdata, ydata, p0=[0, 6.71647395e-04, -0.108, 1191])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt4), 'g--')

#%% And one on the other side of Q1
A1_L5 = []
A1_Err_L5 = []
for i in range(75):
    if i == 32 or i == 36 or i == 39:
        A1_L5.append(pL1.a1+1270)
        A1_Err_L5.append(eL1[1])
    else:
        Sky_No_O4 = Sky_O4[10+i,:]
        a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1290:1303]))
        a1_pointP = a1_pointP[0][0] - 1270
        result_L1 = gmodel.fit(Sky_No_O4[1270:1320], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1290:1303]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        pL1 = SimpleNamespace(**result_L1.best_values)
        if result_L1.covar is None:
            eL1 = 0,0,0
            A1_L5.append(pL1.a1+1270)
            A1_Err_L5.append(eL1[1])
        else:
            eL1 = np.sqrt(np.diag(result_L1.covar))
            A1_L5.append(pL1.a1+1270)
            A1_Err_L5.append(eL1[1])

A1_L5[21] = 1295.6
A1_L5[22] = 1295.8
A1_L5[31] = 1295

ydata = A1_L5
xdata = np.arange(75)

popt5, pcov5 = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 1405])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt5), 'g--')

#%% And another later line will be helpful here!
A1_L3 = []
A1_Err_L3 = []
for i in range(75):
    Sky_No_O4 = Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1506:1520]))
    if a1_pointP[0][0] < 1500:
        a1_pointP = 1514 - 1485
    else:
        a1_pointP = a1_pointP[0][0] - 1485
    result_L3 = gmodel.fit(Sky_No_O4[1485:1535], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1506:1526]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL3 = SimpleNamespace(**result_L3.best_values)
    eL3 = np.sqrt(np.diag(result_L3.covar))
    A1_L3.append(pL3.a1+1485)
    A1_Err_L3.append(eL3[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation
def lin_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

A1_L3[31] = 1513.6
A1_L3[41] = 1513.2
# A1_L1[70] = 24.3
ydata = A1_L3
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(75)

popt3, pcov3 = curve_fit(lin_func, xdata, ydata, p0=[0, 0, -0.032, 26.5])
plt.figure()
plt.plot(np.arange(75), ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt3), 'g--')

#%% Another line if possible
A1_L6 = []
A1_Err_L6 = []
for i in range(75):
    Sky_No_O4 = Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1729:1749]))
    if a1_pointP[0][0] < 1700:
        a1_pointP = 1739 - 1710
    else:
        a1_pointP = a1_pointP[0][0] - 1710
    result_L6 = gmodel.fit(Sky_No_O4[1710:1760], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1729:1749]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL6 = SimpleNamespace(**result_L6.best_values)
    eL6 = np.sqrt(np.diag(result_L6.covar))
    A1_L6.append(pL6.a1+1710)
    A1_Err_L6.append(eL6[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation
def lin_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

A1_L6[23] = 1739.7
A1_L6[32] = 1739.3
A1_L6[33] = 1739.3
ydata = A1_L6
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(75)

popt6, pcov6 = curve_fit(lin_func, xdata, ydata, p0=[0, 0, -0.032, 26.5])
plt.figure()
plt.plot(np.arange(75), ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt6), 'g--')

#%%Now we know what we want to shift the spectra by lets fidget with the Sky_order

import scipy.ndimage as f #This will allow the equivalent of fshift-ing in IDL

#Set up the shift factor(look to automate this)
#corrective_s_factor = (A_average + B_average)
A_average = popt1[0]
B_average = popt1[1]
C_average = popt1[2]
D_average = popt1[3]
Corrective_Array = []

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 85:
        ii = i - 10
        Sky_row_O4 = Sky_O4[i,:]
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 25.404971212701255)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Sky_O4[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
# plt.figure()
# plt.imshow(Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)

A1_L1 = []
A1_Err_L1 = []
for i in range(75):
    # if i == 73 and i == 74:
    #     A1_L1.append(19.377984948979858+10)
    #     A1_Err_L1.append(eL1[1])
    # else:
        Sky_No_O4 = shift_Sky_O4[10+i,:]
        a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[10:40]))
        a1_pointP = a1_pointP[0][0] - 10
        result_L1 = gmodel.fit(Sky_No_O4[10:60], x=np.arange(50), a0=400000, a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        pL1 = SimpleNamespace(**result_L1.best_values)
        eL1 = np.sqrt(np.diag(result_L1.covar))
        A1_L1.append(pL1.a1+10)
        A1_Err_L1.append(eL1[1])
    
#Now to find out the curve in y and x this line has!
#Define a linear equation and a quadratic equation

A1_L1[0] = 25.3
A1_L1[5] = 25.3
A1_L1[70] = 25.4
ydata = A1_L1
xdata = np.arange(75)
XPOS_L1 = np.nanmean(ydata)

popt1A, pcov1A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 31.917])
plt.figure()
plt.plot(np.arange(75), A1_L1, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt1A), 'g--')

#%%
Corrective_Array2 = []
A_average = popt1A[0]
B_average = popt1A[1]
C_average = popt1A[2]
D_average = popt1A[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 95:
        ii = i - 10
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L1)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)

#%%
#corrective_s_factor = (A_average + B_average)
A_average = popt2[0]
B_average = popt2[1]
C_average = popt2[2]
D_average = popt2[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 678.0391483419393)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Sky_O4[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
# plt.figure()
# plt.imshow(Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
A1_L2 = []
A1_Err_L2 = []
for i in range(75):
    Sky_No_O4 = shift_Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[665:685]))
    a1_pointP = a1_pointP[0][0] - 655
    result_L1 = gmodel.fit(Sky_No_O4[655:705], x=np.arange(50), a0=np.nanmax(Sky_No_O4[665:685]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L2.append(pL1.a1+655)
    A1_Err_L2.append(eL1[1])
    
A1_L2[18] = 678.2
A1_L2[0] = 678.0

ydata = A1_L2
xdata = np.arange(75)
XPOS_L2 = np.nanmean(ydata)

popt2A, pcov2A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 580])
plt.figure()
plt.plot(np.arange(75), A1_L2, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt2A), 'g--')

#%%
A_average = popt2A[0]
B_average = popt2A[1]
C_average = popt2A[2]
D_average = popt2A[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L2)*-1
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

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 1187.9681265525958)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Sky_O4[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)
         
A1_L4 = []
A1_Err_L4 = []
for i in range(75):
        Sky_No_O4 = shift_Sky_O4[10+i,:]
        a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1181:1200]))
        a1_pointP = a1_pointP[0][0] - 1170
        result_L1 = gmodel.fit(Sky_No_O4[1170:1220], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1181:1200]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
        pL1 = SimpleNamespace(**result_L1.best_values)
        if result_L1.covar is None:
            eL1 = 0,0,0
            A1_L4.append(pL1.a1+1170)
            A1_Err_L4.append(eL1[1])
        else:
            eL1 = np.sqrt(np.diag(result_L1.covar))
            A1_L4.append(pL1.a1+1170)
            A1_Err_L4.append(eL1[1])
    
A1_L4[2] = 1188.7
A1_L4[53] = 1188.7
A1_L4[15] = 1188.7
ydata = A1_L4
xdata = np.arange(75)
XPOS_L4 = np.nanmean(ydata)

popt4A, pcov4A = curve_fit(lin_func, xdata, ydata, p0=[0, 6.71647395e-04, -0.108, 1191])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt4A), 'g--')

#%%
A_average = popt4A[0]
B_average = popt4A[1]
C_average = popt4A[2]
D_average = popt4A[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 95:
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
# A_average = popt5[0]
# B_average = popt5[1]
# C_average = popt5[2]

# for i in range(len(Sky_O4[:,0])):
#     if i == 0:
#         Sky_row_O4 = Sky_O4[i,:]
#         shift_Sky_O4 = Sky_row_O4
#         Corrective_Array.append(0)
#     elif i < 10:
#         shifted_Sky_O4 = Sky_O4[i,:]
#         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
#         Corrective_Array.append(0)
#     elif i > 10 and i < 95:
#         Sky_row_O4 = Sky_O4[i,:]
#         corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average)*-1
#         shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
#         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
#         Corrective_Array.append(corrective_shift_factor)
#     else:
#           shifted_Sky_O4 = Sky_O4[i,:]
#           shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
#           Corrective_Array.append(0)

# # plt.figure()
# # plt.imshow(shift_Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
# # plt.figure()
# # plt.imshow(Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)

# A1_L5 = []
# A1_Err_L5 = []
# for i in range(75):
#         Sky_No_O4 = shift_Sky_O4[10+i,:]
#         a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1293:1303]))
#         a1_pointP = a1_pointP[0][0] - 1270
#         result_L1 = gmodel.fit(Sky_No_O4[1270:1320], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1293:1303]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
#         pL1 = SimpleNamespace(**result_L1.best_values)
#         if result_L1.covar is None:
#             eL1 = 0,0,0
#         else:
#             eL1 = np.sqrt(np.diag(result_L1.covar))
#         A1_L5.append(pL1.a1+1270)
#         A1_Err_L5.append(eL1[1])

# A1_L5[8] = 1297.6
# A1_L5[47] = 1297.6
# A1_L5[21] = 1297.6
# A1_L5[22] = 1297.6
# A1_L5[24] = 1297.6

# ydata = A1_L5
# xdata = np.arange(75)
# XPOS_L5 = np.nanmean(ydata)

# popt5A, pcov5A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 1405])
# plt.figure()
# plt.plot(xdata, ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt5A), 'g--')

#%%
# A_average = popt5A[0]
# B_average = popt5A[1]
# C_average = popt5A[2]

# for i in range(len(Sky_O4[:,0])):
#     if i == 0:
#         Sky_row_O4 = shift_Sky_O4[i,:]
#         shift_Sky_O4B = Sky_row_O4
#         Corrective_Array2.append(0)
#     elif i < 10:
#         shifted_Sky_O4 = shift_Sky_O4[i,:]
#         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
#         Corrective_Array2.append(0)
#     elif i > 10 and i < 95:
#         Sky_row_O4 = shift_Sky_O4[i,:]
#         corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average)*-1
#         shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
#         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
#         Corrective_Array2.append(corrective_shift_factor)
#     else:
#           shifted_Sky_O4 = shift_Sky_O4[i,:]
#           shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
#           Corrective_Array2.append(0)

#%%
A_average = popt3[0]
B_average = popt3[1]
C_average = popt3[2]
D_average = popt3[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 1513.7053821168158)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
          shifted_Sky_O4 = Sky_O4[i,:]
          shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
          Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
# plt.figure()
# plt.imshow(Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)

A1_L1 = []
A1_Err_L1 = []
for i in range(75):
    Sky_No_O4 = shift_Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1508:1521]))
    if a1_pointP[0][0] < 1500:
        a1_pointP = 1515 - 1485
    else:
        a1_pointP = a1_pointP[0][0] - 1485
    result_L1 = gmodel.fit(Sky_No_O4[1485:1535], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1508:1521]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+1485)
    A1_Err_L1.append(eL1[1])
    
#A1_L1[31] = 1513.6
#A1_L1[41] = 1513.2
# A1_L1[70] = 24.3
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(75)
XPOS_L3 = np.nanmean(ydata)

popt3A, pcov3A = curve_fit(lin_func, xdata, ydata, p0=[0, 0, -0.032, 26.5])
plt.figure()
plt.plot(np.arange(75), ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt3A), 'g--')

#%%
A_average = popt3A[0]
B_average = popt3A[1]
C_average = popt3A[2]
D_average = popt3A[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L3)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
          shifted_Sky_O4 = shift_Sky_O4[i,:]
          shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
          Corrective_Array2.append(0)

#%%
A_average = popt6[0]
B_average = popt6[1]
C_average = popt6[2]
D_average = popt6[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - 1739.3471398206007)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Sky_O4[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

A1_L1 = []
A1_Err_L1 = []
for i in range(75):
    Sky_No_O4 = shift_Sky_O4[10+i,:]
    a1_pointP = np.where(Sky_No_O4 == np.nanmax(Sky_No_O4[1729:1749]))
    if a1_pointP[0][0] < 1700:
        a1_pointP = 1739 - 1710
    else:
        a1_pointP = a1_pointP[0][0] - 1710
    result_L1 = gmodel.fit(Sky_No_O4[1710:1760], x=np.arange(50), a0=np.nanmax(Sky_No_O4[1729:1749]), a1=a1_pointP, a2=2, a3=0, a4=0, a5=0)
    pL1 = SimpleNamespace(**result_L1.best_values)
    eL1 = np.sqrt(np.diag(result_L1.covar))
    A1_L1.append(pL1.a1+1710)
    A1_Err_L1.append(eL1[1])
    
A1_L1[43] = 1740
A1_L1[32] = 1740
A1_L1[33] = 1740
A1_L1[34] = 1740
A1_L1[49] = 1740
A1_L1[50] = 1740
# A1_L1[41] = 1513.2
# A1_L1[70] = 24.3
ydata = A1_L1
#ydata = np.clip(ydata, 26, 27)
xdata = np.arange(75)
XPOS_L6 = np.nanmean(A1_L1)

popt6A, pcov6A = curve_fit(lin_func, xdata, ydata, p0=[0, 0, -0.032, 26.5])
plt.figure()
plt.plot(np.arange(75), ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt6A), 'g--')

#%% I think there's still a second element to to sort so we'll make a second corrective array
A_average = popt6A[0]
B_average = popt6A[1]
C_average = popt6A[2]
D_average = popt6A[3]

for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 10:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i > 10 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        ii = i - 10
        corrective_shift_factor = ((ii**3)*A_average + (ii**2)*B_average + (ii)*C_average + D_average - XPOS_L6)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
        Corrective_Array2.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = shift_Sky_O4[i,:]
         shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4))
         Corrective_Array2.append(0)

#%% Now to make multiple
Corrections = np.array(Corrective_Array)
Corrections2 = np.array(Corrective_Array2)
CorrectionsA = Corrections.reshape((5, 110))
CorrectionsB = Corrections2.reshape((5, 110))
X_Positions = np.asarray([np.round(XPOS_L1), np.round(XPOS_L2), np.round(XPOS_L4), np.round(XPOS_L3), np.round(XPOS_L6)])

def lin_func(x, a , b, c):
    y = a* x**2 + b*x + c
    return y

Corrections_Spectra = []
gmodel = Model(lin_func)
Margin2 = []

for m in range(110):
    poptT, pcovT = curve_fit(lin_func, X_Positions, CorrectionsA[:,m], p0=[0, 0, 2])
    #plt.figure()
    #plt.plot(X_Positions, CorrectionsA[:,m+20], 'bo')
    #best_val = SimpleNamespace(**result.best_values)
    #plt.plot(np.arange(2048),lin_func(np.arange(2048), best_val.a, best_val.b, best_val.c), 'g--')
    #plt.plot(X_Positions, result.init_fit, 'r--')
    #plt.plot(X_Positions, result.best_fit, 'g--')
    #print(best_val.a, best_val.b, best_val.c, best_val.d)
    Corrections_Spectra.append(lin_func(np.arange(2048), *poptT))
    # print(lin_func(2013, best_val.a, best_val.b, best_val.c))
    #Margin2.append(lin_func(2013, best_val.a, best_val.b, best_val.c, best_val.d))
    
plt.figure()
plt.imshow(Corrections_Spectra)
plt.colorbar()

Corrections_SpectraB = []

for m in range(110):
    poptTB, pcovTB = curve_fit(lin_func, X_Positions, CorrectionsB[:,m], p0=[0, 0, 0.5])
    # plt.figure()
    # plt.plot(X_Positions, CorrectionsA[:,m], 'bo')
    # plt.plot(np.arange(2048), lin_func(np.arange(2048), *poptT), 'g--')
    Corrections_SpectraB.append(lin_func(np.arange(2048), *poptTB))
    
plt.figure()
plt.imshow(Corrections_SpectraB)
plt.colorbar()

Corrections_SpectraT = np.add(Corrections_SpectraB, Corrections_Spectra)

#%% Now we use the Correction Spectra to fix the spectra just for the Q1 spectra
for i in range(len(Sky_O4[:,0])):
    if i == 0:
        Sky_row_O4 = Sky_O4[i,:]
        shift_Sky_O4 = Sky_row_O4
    elif i < 10:
        shifted_Sky_O4 = Sky_O4[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
    elif i > 10 and i < 95:
        Sky_row_O4 = Sky_O4[i,:]
        corrective_shift_factor = np.nanmean(Corrections_SpectraT[i][1205:1225])
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
    else:
         shifted_Sky_O4 = Sky_O4[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))

np.save('Corrective_Array_O4', Corrections_SpectraT)
plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)
plt.figure()
plt.imshow(Sky_O4, cmap='gist_gray', vmax = 600000, vmin = -600000)

#%%Now to find out the wavelength on the O4 order ## Going to need more lines here!!!!!!!!!!!!!!!!!!!
Waves = data.T[0]
Trans = data.T[1]
S_wave_O4 = np.where(Waves == 3.93270)
F_wave_O4 = np.where(Waves == 3.96560)
S_waveO4 = S_wave_O4[0][0]
F_waveO5 = F_wave_O4[0][0]
Waves_O4 = Waves[S_waveO4:F_waveO4]
Trans_O4 = -1*Trans[S_waveO4:F_waveO4]
gmodel = Model(gauss_fit)

a1_pointP = np.where(Waves_O4[10:60] == np.nanmax(Waves_O4[10:60]))
a1_pointP = a1_pointP[0][0]
a1_pointP = Waves_O4[a1_pointP]
result_L1_O4 = gmodel.fit(Trans_O4[10:60], x=Waves_O4[10:60], a0=-1*np.nanmax(Trans_O4[10:60]), a1=a1_pointP, a2=0.0015, a3=-0.90, a4=0, a5=0)
pL1_O4 = SimpleNamespace(**result_L1_O4.best_values)
eL1_O4 = np.sqrt(np.diag(result_L1_O4.covar))
A1_L1_O4 = pL1_O4.a1
A1_Err_L1_O4 = eL1_O4[1]

print(A1_L1_O4, A1_Err_L1_O4)
#print(np.nanmean(A1_L1))
plt.figure()
plt.plot(Waves_O4[10:60], Trans_O4[10:60], 'b')
plt.plot(Waves_O4[10:60], result_L1_O4.best_fit, 'g-', label='best fit')
plt.close()

#Next line (2 is too close to a second line)
result_L2_O4 = gmodel.fit(Trans_O4[550:590], x=Waves_O4[550:590], a0=0.15, a1=3.94436, a2=0.0016, a3=-0.9, a4=0, a5=0)
pL2_O4 = SimpleNamespace(**result_L2_O4.best_values)
eL2_O4 = np.sqrt(np.diag(result_L2_O4.covar))
A1_L2_O4 = pL2_O4.a1
A1_Err_L2_O4 = eL2_O4[1]

print(A1_L2_O4, A1_Err_L2_O4)
print(np.nanmean(A1_L2))
plt.figure()
plt.plot(Waves_O4[550:590], Trans_O4[550:590], 'b')
plt.plot(Waves_O4[550:590], result_L2_O4.best_fit, 'g-', label='best fit')
plt.close()

#Next line 
result_L4_O4 = gmodel.fit(Trans_O4[970:1020], x=Waves_O4[970:1020], a0=0.1, a1=3.9526, a2=0.002, a3=-0.89, a4=0, a5=0)
pL4_O4 = SimpleNamespace(**result_L4_O4.best_values)
eL4_O4 = np.sqrt(np.diag(result_L4_O4.covar))
A1_L4_O4 = pL4_O4.a1
A1_Err_L4_O4 = eL4_O4[1]

print(A1_L4_O4, A1_Err_L4_O4)
print(np.nanmean(A1_L4))
plt.figure()
plt.plot(Waves_O4[970:1020], Trans_O4[970:1020], 'b')
plt.plot(Waves_O4[970:1020], result_L4_O4.best_fit, 'g-', label='best fit')
plt.close()

# result_L5_O4 = gmodel.fit(Trans_O4[1055:1105], x=Waves_O4[1055:1105], a0=0.1, a1=3.9542, a2=0.002, a3=-0.89, a4=0, a5=0)
# pL5_O4 = SimpleNamespace(**result_L5_O4.best_values)
# eL5_O4 = np.sqrt(np.diag(result_L5_O4.covar))
# A1_L5_O4 = pL5_O4.a1
# A1_Err_L5_O4 = eL5_O4[1]
# print(A1_L5_O4, A1_Err_L5_O4)
# print(np.nanmean(A1_L5))

result_L3_O4 = gmodel.fit(Trans_O4[1215:1265], x=Waves_O4[1215:1265], a0=0.1, a1=3.9576, a2=0.002, a3=-0.89, a4=0, a5=0)
pL3_O4 = SimpleNamespace(**result_L3_O4.best_values)
eL3_O4 = np.sqrt(np.diag(result_L3_O4.covar))
A1_L3_O4 = pL3_O4.a1
A1_Err_L3_O4 = eL3_O4[1]
print(A1_L3_O4, A1_Err_L3_O4)
print(np.nanmean(A1_L3))

result_L6_O4 = gmodel.fit(Trans_O4[1385:1435], x=Waves_O4[1385:1435], a0=0.1, a1=3.96105, a2=0.002, a3=-0.89, a4=0, a5=0)
pL6_O4 = SimpleNamespace(**result_L6_O4.best_values)
eL6_O4 = np.sqrt(np.diag(result_L6_O4.covar))
A1_L6_O4 = pL6_O4.a1
A1_Err_L6_O4 = eL6_O4[1]
print(A1_L6_O4, A1_Err_L6_O4)
print(np.nanmean(A1_L6))

#%%
#Calculate the ratio
Waves = [A1_L1_O4, A1_L2_O4, A1_L4_O4, A1_L3_O4, A1_L6_O4]
PixWav = X_Positions

def poly_func(x, a , b, c):
    y = a* x**2 + b*x + c
    return y

poptW, pcovW = curve_fit(poly_func, X_Positions, Waves, p0=[0, 0, 3.96])
plt.figure()
plt.plot(X_Positions, Waves, 'bo')
    #best_val = SimpleNamespace(**result.best_values)
plt.plot(np.arange(2048),poly_func(np.arange(2048), *poptW), 'g--')

#Lets find Q1's expected position!
# poptWA, poptWA = curve_fit(poly_func, Waves, X_Positions,  p0=[0, 0, 100])
# print(poly_func(3.95295, *poptWA ))
# X_Q1 = (-1*poptW[1] + np.sqrt((poptW[1])**2 - 4*poptW[0]*(poptW[2]-3.95295)))/2*poptW[0]
# print(X_Q1)
# X_Q1 = (-1*poptW[1] - np.sqrt((poptW[1])**2 - 4*poptW[0]*(poptW[2]-3.95295)))/2*poptW[0]
# print(X_Q1)

#%% Now one final try to see if there's a linear relationship
