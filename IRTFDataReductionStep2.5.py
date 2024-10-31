# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 15:18:55 2022

@author: snowy
"""

from astropy.io import fits
from lmfit import Model
import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace
import math
from lmfit.models import GaussianModel, ExponentialModel
from scipy.optimize import curve_fit

#As Order 5 is difficult to do with sky lines we switch to using ThAr lamp lines

filename = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/CR_GCAT_080901A_lines_thar.fits'
ThAr_lines = fits.getdata(filename)

ThAr_lines_Wav = []
ThAr_lines_Emission = []
for a in range(len(ThAr_lines)):
    ThAr_lines_Wav.append(ThAr_lines[a][0])
    ThAr_lines_Emission.append(ThAr_lines[a][1])

plt.figure()
plt.plot(ThAr_lines_Wav, ThAr_lines_Emission)

#This is the emission spectra for the ThAr lamp hence we bring up the arc lines spectra again and 1) measure the position of the lines and two assign a wavelength to pixels

filename = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/arcs.order.5.fits.gz'
Arc_linesO5 = fits.getdata(filename)

plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

#%% 
#Line 1 = 193
#Line 2 = 230
#Line 3 = 1234
#Line 4 = 1709
#Line 5 = 1713

Total_Arc_lines = np.zeros(2048)

#First lets make a big plot
for a in range(110):
    Total_Arc_lines += Arc_linesO5[a,:]
    
plt.figure()
plt.plot(np.arange(2048), Total_Arc_lines)
#Lets start at the end with two lines

# #%% So lets focus on 1660 to 1760

# ydata = Arc_linesO5[:,1660:1760]
# Line4A1 = []
# Line5A1 = []

# for a in range(75):
#     gauss1 = GaussianModel(prefix='g1_')
#     pars = gauss1.guess(ydata[a+10,:], x=np.arange(100))
    
#     pars['g1_center'].set(value=49, min=46, max=53)
#     pars['g1_sigma'].set(value=3, min=0)
#     pars['g1_amplitude'].set(value=35000, min=1000)

#     gauss2 = GaussianModel(prefix='g2_')
#     pars.update(gauss2.make_params())

#     pars['g2_center'].set(value=56, min=53, max=59)
#     pars['g2_sigma'].set(value=3, min=0)
#     pars['g2_amplitude'].set(value=35000, min=1000)
    
#     mod = gauss1 + gauss2
    
#     init = mod.eval(pars, x=np.arange(100))
#     out = mod.fit(ydata[a+10,:], pars, x=np.arange(100))
#     Lines_Vals = out.best_values
#     Line4A1.append(Lines_Vals['g1_center']+1660)
#     Line5A1.append(Lines_Vals['g2_center']+1660)
#     print('Loop ' + str(a))

#%%
#Now focus on 1180 to 1280
#First lets fit the Q1, Q2, Q3, Q3,1 and Q3,2 if possible (Only Q3,1 possible)
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Line3A1 = []
gmodel = Model(gauss_fit)

for a in range(75):
    Arcs_data = Arc_linesO5[a+10,:]
    a1_point = np.where(Arcs_data[1228:1240] == np.nanmax(Arcs_data[1228:1240]))
    a1_point = a1_point[0][0] + 28
    result_L3 = gmodel.fit(Arcs_data[1200:1250], x=np.arange(50), a0=np.nanmax(Arcs_data[1228:1240]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL3A = SimpleNamespace(**result_L3.best_values)
    Line3A1.append(pL3A.a1+1200)

#%% Now focus on 200 to 260
Line2A1 = []

for a in range(75):
    Arcs_data = Arc_linesO5[a+10,:]    
    a1_point = np.where(Arcs_data[220:240] == np.nanmax(Arcs_data[220:240]))
    a1_point = a1_point[0][0] + 20
    result_L2 = gmodel.fit(Arcs_data[200:260], x=np.arange(60), a0=np.nanmax(Arcs_data[200:260]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL2A = SimpleNamespace(**result_L2.best_values)
    Line2A1.append(pL2A.a1+200)
    
#%% Now focus on 160 to 220
Line1A1 = []

for a in range(75):
    Arcs_data = Arc_linesO5[a+10,:]
    a1_point = np.where(Arcs_data[180:200] == np.nanmax(Arcs_data[180:200]))
    a1_point = a1_point[0][0] + 20
    result_L1 = gmodel.fit(Arcs_data[160:220], x=np.arange(60), a0=np.nanmax(Arcs_data[160:220]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL1A = SimpleNamespace(**result_L1.best_values)
    Line1A1.append(pL1A.a1+160)
    
Line4A1 = []
for a in range(75):
    if a == 14 or a == 37 or a == 41:
        pass
    else:
        Arcs_data = Arc_linesO5[a+10,:]
        a1_point = np.where(Arcs_data[1680:1730] == np.nanmax(Arcs_data[1680:1730]))
        a1_point = a1_point[0][0]
        result_L4 = gmodel.fit(Arcs_data[1680:1730], x=np.arange(50), a0=np.nanmax(Arcs_data[1680:1730]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
        pL4A = SimpleNamespace(**result_L4.best_values)
        Line4A1.append(pL4A.a1+1680)
    # if a > 20 and a < 50:
    #     plt.figure()
    #     plt.plot(np.arange(50), Arcs_data[1680:1730])
    #     plt.plot(np.arange(50), result_L4.best_fit)
    
Line5A1 = []
for a in range(75):
        Arcs_data = Arc_linesO5[a+10,:]
        a1_point = np.where(Arcs_data[685:705] == np.nanmax(Arcs_data[685:705]))
        a1_point = a1_point[0][0] + 15
        result_L5 = gmodel.fit(Arcs_data[670:720], x=np.arange(50), a0=np.nanmax(Arcs_data[685:705]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
        pL5A = SimpleNamespace(**result_L5.best_values)
        Line5A1.append(pL5A.a1+670)
    
#%% Now we have a rough idea of how skewed everything is we now need to fshift it back into position
def lin_func(x, a , b, c, d):
    y = a* x**3 + b*x**2 + c*x + d
    return y

# plt.figure()
# plt.plot(np.arange(75), Line3A1)

# plt.figure()
# plt.plot(np.arange(75), Line2A1)

# plt.figure()
# plt.plot(np.arange(75), Line1A1)

#Now to work out the quadratic equations of these curves and take the a and b values over
Line3A1[38] = 1233.9
Line3A1[39] = 1233.8
Line3A1[40] = 1233.8
Line3A1[41] = 1233.8
ydata = Line3A1
xdata = np.arange(75)

popt3, pcov3 = curve_fit(lin_func, xdata, ydata, p0=[0, 2.3e-4, -6.3e-2, 1.23588853e3])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt3), 'g--')

Line2A1[18] = 230.8
Line2A1[19] = 230.7
ydata = Line2A1
xdata = np.arange(75)

popt2, pcov2 = curve_fit(lin_func, xdata, ydata, p0=[0, 1.15e-4, -4.7e-2, 2.32e2])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt2), 'g--')

Line1A1[73] = 191.2
ydata = Line1A1
xdata = np.arange(75)

popt1, pcov1 = curve_fit(lin_func, xdata, ydata, p0=[0, 1.5762e-4, -4.86e-2, 1.94e2])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt1), 'g--')

Line4A1[35] = 1711
ydata = Line4A1
xdata = np.concatenate((np.arange(14), np.arange(15,37), np.arange(38, 41), np.arange(42, 75)))

popt4, pcov4 = curve_fit(lin_func, xdata, ydata, p0=[0, 1.5762e-4, -4.86e-2, 1.94e2])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt4), 'g--')

Line5A1[40] = 695.1
Line5A1[42] = 695.2
ydata = Line5A1
xdata = np.arange(75)

popt5, pcov5 = curve_fit(lin_func, xdata, ydata, p0=[0, 1.5762e-4, -4.86e-2, 1.94e2])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt5), 'g--')

# A_shift = (popt1[0]+popt2[0]+popt3[0])/3
# B_shift = (popt1[1]+popt2[1]+popt3[1])/3

# print(A_shift)
# print(B_shift)
#%%
import scipy.ndimage as f #This will allow the equivalent of fshift-ing in IDL

#Set up the shift factor(look to automate this)

#corrective_s_factor = (A_average + B_average)
A_average = popt1[0]
B_average = popt1[1]
C_average = popt1[2]
D_average = popt1[3]
Corrective_Array = []

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - 192.4693895462319)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

Line1A1 = []

for a in range(75):
    Arcs_data = shift_Sky_O4[a+10,:]
    a1_point = np.where(Arcs_data[180:200] == np.nanmax(Arcs_data[180:200]))
    a1_point = a1_point[0][0] + 20
    result_L1 = gmodel.fit(Arcs_data[160:220], x=np.arange(60), a0=np.nanmax(Arcs_data[160:220]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL1A = SimpleNamespace(**result_L1.best_values)
    Line1A1.append(pL1A.a1+160)

Line1A1[73] = 192.8
ydata = Line1A1
xdata = np.arange(75)
XPOS_L1 = np.nanmean(ydata)

popt1A, pcov1A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 170])
plt.figure()
plt.plot(np.arange(75), Line1A1, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt1A), 'g--')

#%%
Corrective_Array2 = []
A_average = popt1A[0]
B_average = popt1A[1]
C_average = popt1A[2]
D_average = popt1A[3]

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 5:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - XPOS_L1)*-1
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

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - 230.05115184986138)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

Line2A1 = []

for a in range(75):
    Arcs_data = shift_Sky_O4[a+10,:]    
    a1_point = np.where(Arcs_data[220:240] == np.nanmax(Arcs_data[220:240]))
    a1_point = a1_point[0][0] + 20
    result_L2 = gmodel.fit(Arcs_data[200:260], x=np.arange(60), a0=np.nanmax(Arcs_data[200:260]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL2A = SimpleNamespace(**result_L2.best_values)
    Line2A1.append(pL2A.a1+200)

Line2A1[18] = 230.4
Line2A1[19] = 230.4
ydata = Line2A1
xdata = np.arange(75)
XPOS_L2 = np.nanmean(ydata)

popt2A, pcov2A = curve_fit(lin_func, xdata, ydata, p0=[0, 0.00011, -0.028, 170])
plt.figure()
plt.plot(np.arange(75), Line2A1, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt2A), 'g--')

#%%
A_average = popt2A[0]
B_average = popt2A[1]
C_average = popt2A[2]
D_average = popt2A[3]

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 5:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - XPOS_L2)*-1
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

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - 1234.020616362938)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

Line3A1 = []

for a in range(75):
    Arcs_data = shift_Sky_O4[a+10,:]
    a1_point = np.where(Arcs_data[1228:1240] == np.nanmax(Arcs_data[1228:1240]))
    a1_point = a1_point[0][0] + 28
    result_L3 = gmodel.fit(Arcs_data[1200:1250], x=np.arange(50), a0=np.nanmax(Arcs_data[1228:1240]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
    pL3A = SimpleNamespace(**result_L3.best_values)
    Line3A1.append(pL3A.a1+1200)

Line3A1[20] = 1234.5
Line3A1[38] = 1234.5
Line3A1[39] = 1234.5
Line3A1[40] = 1234.5
Line3A1[41] = 1234.5
Line3A1[71] = 1234.3
ydata = Line3A1
xdata = np.arange(75)
XPOS_L3 = np.nanmean(ydata)

popt3A, pcov3A = curve_fit(lin_func, xdata, ydata, p0=[0, 2.3e-4, -6.3e-2, 1.23588853e3])
plt.figure()
plt.plot(xdata, ydata, 'bo')
#plt.plot(xdata, Att1, 'g--')
#Fitdata = lin_func(A)
plt.plot(xdata, lin_func(xdata, *popt3A), 'g--')

#%%
A_average = popt3A[0]
B_average = popt3A[1]
C_average = popt3A[2]
D_average = popt3A[3]

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 5:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - XPOS_L3)*-1
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

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
        Corrective_Array.append(0)
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
        Corrective_Array.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - 1711.0559336957474)*-1
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
        Corrective_Array.append(corrective_shift_factor)
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
         Corrective_Array.append(0)

plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

Line4A1 = []
for a in range(75):
    if a == 19:
        Line4A1.append(17.322738846444+1680)
    else:
        Arcs_data = shift_Sky_O4[a+10,:]
        a1_point = np.where(Arcs_data[1700:1720] == np.nanmax(Arcs_data[1700:1720]))
        a1_point = a1_point[0][0]
        result_L4 = gmodel.fit(Arcs_data[1680:1730], x=np.arange(50), a0=np.nanmax(Arcs_data[1700:1720]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
        pL4A = SimpleNamespace(**result_L4.best_values)
        Line4A1.append(pL4A.a1+1680)

Line4A1[41] = 1696.9
Line4A1[43] = 1696.8
ydata = Line4A1
xdata = np.arange(75)
XPOS_L4 = np.nanmean(ydata)

popt4A, pcov4A = curve_fit(lin_func, xdata, ydata, p0=[0, 2.3e-4, -6.3e-2, 1.23588853e3])
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

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = Sky_row_O4
        Corrective_Array2.append(0)
    elif i < 5:
        shifted_Sky_O4 = shift_Sky_O4[i,:]
        shift_Sky_O4B = np.vstack((shift_Sky_O4B, shifted_Sky_O4)) 
        Corrective_Array2.append(0)
    elif i >= 5 and i < 95:
        Sky_row_O4 = shift_Sky_O4[i,:]
        corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average + D_average - XPOS_L4)*-1
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

# for i in range(len(Arc_linesO5[:,0])):
#     if i == 0:
#         Sky_row_O4 = Arc_linesO5[i,:]
#         shift_Sky_O4 = Sky_row_O4
#         Corrective_Array.append(0)
#     elif i < 5:
#         shifted_Sky_O4 = Arc_linesO5[i,:]
#         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
#         Corrective_Array.append(0)
#     elif i >= 5 and i < 95:
#         Sky_row_O4 = Arc_linesO5[i,:]
#         corrective_shift_factor = ((i**3)*A_average + (i**2)*B_average + (i)*C_average)*-1
#         shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
#         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
#         Corrective_Array.append(corrective_shift_factor)
#     else:
#          shifted_Sky_O4 = Arc_linesO5[i,:]
#          shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
#          Corrective_Array.append(0)

# plt.figure()
# plt.imshow(shift_Sky_O4, cmap='gist_gray')
# plt.figure()
# plt.imshow(Arc_linesO5, cmap='gist_gray')

# Line5A1 = []
# for a in range(75):
#     if a == 8 or a == 11:
#         Line5A1.append(698.5)
#     else:
#         Arcs_data = shift_Sky_O4[a+10,:]
#         a1_point = np.where(Arcs_data[690:700] == np.nanmax(Arcs_data[690:700]))
#         a1_point = a1_point[0][0] + 15
#         result_L5 = gmodel.fit(Arcs_data[670:720], x=np.arange(50), a0=np.nanmax(Arcs_data[690:700]), a1=a1_point, a2=2, a3=0, a4=0, a5=0)
#         pL5A = SimpleNamespace(**result_L5.best_values)
#         Line5A1.append(pL5A.a1+670)

# ydata = Line5A1
# xdata = np.arange(75)
# XPOS_L5 = np.nanmean(ydata)

# popt5A, pcov5A = curve_fit(lin_func, xdata, ydata, p0=[0, 2.3e-4, -6.3e-2, 1.23588853e3])
# plt.figure()
# plt.plot(xdata, ydata, 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# plt.plot(xdata, lin_func(xdata, *popt5A), 'g--')

#%%
Corrections = np.array(Corrective_Array)
Corrections2 = np.array(Corrective_Array2)
CorrectionsA = Corrections.reshape((4, 110))
CorrectionsB = Corrections2.reshape((4, 110))
X_Positions = np.asarray([round(XPOS_L1), round(XPOS_L2), round(XPOS_L3), round(XPOS_L4)])

def lin_func(x, a , b, c):
    y = a* x**2 + b*x + c
    return y

Corrections_Spectra = []

for m in range(110):
    poptT, pcovT = curve_fit(lin_func, X_Positions, CorrectionsA[:,m], p0=[0, 0, 25.7])
    # plt.figure()
    # plt.plot(X_Positions, CorrectionsA[:,m], 'bo')
    # plt.plot(np.arange(2048), lin_func(np.arange(2048), *poptT), 'g--')
    Corrections_Spectra.append(lin_func(np.arange(2048), *poptT))
    
plt.figure()
plt.imshow(Corrections_Spectra)
plt.colorbar()

Corrections_SpectraB = []

for m in range(110):
    poptTB, pcovTB = curve_fit(lin_func, X_Positions, CorrectionsB[:,m], p0=[0, 0, 25.7])
    # plt.figure()
    # plt.plot(X_Positions, CorrectionsB[:,m+10], 'bo')
    # plt.plot(np.arange(2048), lin_func(np.arange(2048), *poptTB), 'g--')
    Corrections_SpectraB.append(lin_func(np.arange(2048), *poptTB))
    
plt.figure()
plt.imshow(Corrections_SpectraB)
plt.colorbar()

Corrections_SpectraT = np.add(Corrections_SpectraB, Corrections_Spectra)

#%%
for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = np.nanmean(Corrections_SpectraT[i][1350:1370])
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))

np.save('Corrective_Array_O5', Corrections_SpectraT) #For Q3
plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')

#Now make the figure for arclines
from skimage import exposure
Start = np.nanmean(Arc_linesO5[10:85, :], axis = 0)
p_lower, p_higher = np.percentile(Start, (0.01, 99.9))
Linecut_Arc = exposure.rescale_intensity(Start, in_range=(p_lower, p_higher))

plt.figure(figsize=(12,7))
plt.plot(Linecut_Arc , color='k')
plt.xlabel(r'Wavelength Axis (1.6286 x 10$^{-5}$ $\mu$m per pixel)', fontsize=20)
plt.ylabel('Normalised Detector counts', fontsize=20)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

for i in range(len(Arc_linesO5[:,0])):
    if i == 0:
        Sky_row_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = Sky_row_O4
    elif i < 5:
        shifted_Sky_O4 = Arc_linesO5[i,:]
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) 
    elif i >= 5 and i < 95:
        Sky_row_O4 = Arc_linesO5[i,:]
        corrective_shift_factor = np.nanmean(Corrections_SpectraT[i][1445:1465])
        shifted_Sky_O4 = f.shift(Sky_row_O4, shift=[corrective_shift_factor], mode='wrap')
        shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4))
    else:
         shifted_Sky_O4 = Arc_linesO5[i,:]
         shift_Sky_O4 = np.vstack((shift_Sky_O4, shifted_Sky_O4)) #For Q3,1

plt.figure()
plt.imshow(shift_Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_linesO5, cmap='gist_gray')
#%% Now we work out the arc lines and calculate the wavelength

gmodel = Model(gauss_fit)
A1_L1_O4 = 3.9663502
A1_L2_O4 = 3.96699509
A1_L3_O4 = 3.98353121
A1_L4_O4 = 3.9909353
Waves = [A1_L1_O4, A1_L2_O4, A1_L3_O4, A1_L4_O4]

poptW, pcovW = curve_fit(lin_func, X_Positions, Waves, p0=[0, 1.73130990e-05, 3.967])

plt.figure()
plt.plot(X_Positions, Waves, 'bo')
plt.plot(np.arange(2048), lin_func(np.arange(2048), *poptW), 'g--')
WavelengthsO5 = lin_func(np.arange(2048), *poptW)

