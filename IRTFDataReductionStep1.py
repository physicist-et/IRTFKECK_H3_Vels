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

#First lets fit the Q1, Q2, Q3, Q3,1 and Q3,2 if possible (Only Q3,1 possible)
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

Sky_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.4.fits.gz'
Sky_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/sky.order.5.fits.gz'
Wave_info3 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.3.fits.gz'
Wave_info4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.4.fits.gz'
Wave_info5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.5.fits.gz'
Wave_info6 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.6.fits.gz'
Wave_info7 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration/wavelength.order.7.fits.gz'

Wavelength_O3 = fits.getdata(Wave_info3, ext=0)
Wavelength_O4 = fits.getdata(Wave_info4, ext=0)
Wavelength_O5 = fits.getdata(Wave_info5, ext=0)

Sky_O4 = fits.getdata(Sky_info4, ext=0)
Sky_O5 = fits.getdata(Sky_info5, ext=0)

plt.imshow(Sky_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Sky_O5, cmap='gist_gray')

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
plt.xlim(3.93290, 3.96550)

m= 25
Tot_Sky_O4 = np.zeros(2048)

for i in range(50):
    T_Sky_O4 = Sky_O4[m+i,:]
    Tot_Sky_O4 = np.add(Tot_Sky_O4, T_Sky_O4)

plt.figure()
plt.plot(np.arange(2048), (Tot_Sky_O4)/50)

#%% Match 5 lines across both orders to get a better wavelength

#Line 1 (First), Line 2 (Fifth), Line 3 (Eigth), Line 4 (Tweleth), Line 5 (Sixteen)

Emission_O4 = Trans_O4[0:50]

gmodel = Model(gauss_fit)
result_L1 = gmodel.fit(Emission_O4, x=Waves_O4[0:50], a0=-0.50, a1=3.93336328, a2=0.00018, a3=0, a4=0, a5=0)
pL1 = SimpleNamespace(**result_L1.best_values)
Line_1 = pL1.a1

result_L1A = gmodel.fit(Tot_Sky_O4[10:40], x=np.arange(30), a0=4300000, a1=15.0689762, a2=2, a3=0, a4=0, a5=0)
pL1A = SimpleNamespace(**result_L1A.best_values)
Line_1A = pL1A.a1 + 10

Emission_O4 = Trans_O4[1210:1260]

result_L5 = gmodel.fit(Emission_O4, x=Waves_O4[1210:1260], a0=-0.07230124, a1=3.957618, a2=0.00162683, a3=0.89, a4=0, a5=0)
pL5 = SimpleNamespace(**result_L5.best_values)
Line_5 = pL5.a1

result_L5A = gmodel.fit(Tot_Sky_O4[1720:1760], x=np.arange(40), a0=55000, a1=18.7041288, a2=2.64, a3=0, a4=0, a5=0)
pL5A = SimpleNamespace(**result_L5A.best_values)
Line_5A = pL5A.a1 + 1720

Ratio = (Line_5 - Line_1)/(Line_5A - Line_1A)
Start_W4 = Line_1 - Line_1A*Ratio
print(Start_W4)
Finish_W4 = Start_W4 + 2048*Ratio
print(Finish_W4)

# =============================================================================
# print(result_L5A.fit_report())
# 
# plt.figure()
# plt.plot(np.arange(40), Tot_Sky_O4[1720:1760], 'bo')
# plt.plot(np.arange(40), result_L5A.init_fit, 'k--', label='initial fit')
# plt.plot(np.arange(40), result_L5A.best_fit, 'r-', label='best fit')
# plt.legend(loc='best')
# plt.show()
# =============================================================================

# =============================================================================
# plt.figure()
# plt.plot(Waves_O4[1210:1260], Emission_O4, 'bo')
# plt.plot(Waves_O4[1210:1260], result_L5.init_fit, 'k--', label='initial fit')
# plt.plot(Waves_O4[1210:1260], result_L5.best_fit, 'r-', label='best fit')
# plt.legend(loc='best')
# plt.show()
# =============================================================================

#%% Now plot for the Fifth order

S_wave_O5 = np.where(Waves == 3.96314)
F_wave_O5 = np.where(Waves == 3.99624)
S_waveO5 = S_wave_O5[0][0]
F_waveO5 = F_wave_O5[0][0]
Waves_O5 = Waves[S_waveO5:F_waveO5]
Trans_O5 = Trans[S_waveO5:F_waveO5]

plt.figure()
plt.plot(Waves_O5, Trans_O5)
plt.xlabel(r'Wavelength ($\mu$m)')
plt.ylabel(r'Transmission Value')
plt.ylim(0, 1.0) #show the plot of Wavelength vs Transmission
plt.xlim(3.96313, 3.99624)

m= 20
Tot_Sky_O5 = np.zeros(2048)

for i in range(60):
    T_Sky_O5 = Sky_O5[m+i,:]
    Tot_Sky_O5 = np.add(Tot_Sky_O5, T_Sky_O5)

plt.figure()
plt.plot(np.arange(2048), (Tot_Sky_O5)/60)

#%% Now we look at order 4 and order 5 arc lines

Arc_lineo4 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/arcs.order.4.fits.gz'
Arc_lineo5 = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Calibration2/arcs.order.5.fits.gz'

Arc_O4 = fits.getdata(Arc_lineo4, ext=0)
Arc_O5 = fits.getdata(Arc_lineo5, ext=0)

plt.imshow(Arc_O4, cmap='gist_gray')
plt.figure()
plt.imshow(Arc_O5, cmap='gist_gray')
