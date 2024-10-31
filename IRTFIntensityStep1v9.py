# -*- coding: utf-8 -*-
"""
Created on Fri Jun 16 11:10:30 2023

@author: snowy
"""

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

CombinedImage = np.load('ABBA_Sets_111016.npy')
# Combine sets together to get the intensity and velocity

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
#%%
from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

#Lets clean the image a bit more
CleaningImage = CombinedImage[8]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
CleaningImage[31:34,1209:1211] = np.nanmean(CombinedImage[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[9]
Wavelength_O4 = np.load('Wavelength_O4.npy')
wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225]))

dataI = (CleaningImage + CleaningImage2)/2
plt.figure()
plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
#%% Now lets extract the line and use it for a for the los velocities
image = dataI[:40, 1205:1226]
image = np.fliplr(np.rot90(image, k = 1))
import matplotlib.ticker as tkr
wave = np.load('Wavelength_O4.npy')

plt.figure(figsize = (15,8))
plt.imshow(np.array(image)*1000, cmap='gist_heat', vmax = 11.1, vmin = 0)
plt.title('a) Average $H_{3}$$^{+}$ emission line from iSHELL', fontsize = 30, pad = 20)
plt.xticks(np.arange(0, 40, 5), labels=['50', '55', '60', '65', '70', '75', '80', '85'], fontsize = 20)
plt.xlabel('Spatial Axis row position (0.17")', fontsize = 25)
cbar = plt.colorbar(format=tkr.FormatStrFormatter('%.2f'), shrink = 0.75)
cbar.ax.tick_params(labelsize = 20)
cbar.ax.set_ylabel('Spectral Radiance (m$Wm^{-2}sr^{-1}$$\mu$m$^{-1}$)', fontsize = 25)
labels = str("{:.5f}".format(wave[1225])), str("{:.5f}".format(wave[1220])), str("{:.5f}".format(wave[1215])), str("{:.5f}".format(wave[1210])), str("{:.5f}".format(wave[1205]))
plt.yticks([0, 5, 10, 15, 20], labels = labels, fontsize=20)
plt.xticks(fontsize = 20)
plt.ylabel('Wavelength ($\mu$m)', fontsize = 25, labelpad = 10)
plt.xlabel('Spatial axis row position (pixels)', fontsize = 25)

#%%
dataI[dataI < -0.004] = 0
for m in range(25): #31
    datai = dataI[7+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
    a1_pointP = a1_pointP[0][0]
    ex_resultQPAS = gmodel.fit(datai, x=np.arange(200), a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
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

for o in range(25):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(25):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[0:o+24]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)

Velocity_QAS = []
lamda_H3_O4 = 3.965
lamda_H3_O5 = 3.98558
c = 299792458 
for n in range(25):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(25):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(25):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(25):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

# Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

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
# for o in range(26):
#     Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2))/ (17.24*60*60)
#     Limb_velocity_pix_O4 = (Period_Uranus/1000)
#     Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

# #So the speed of Uranus at surface is 
# # Look into how to vary this across the slit so it can be used for 2016
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []

# Planet_diameter = []
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []
# for o in range(26):
#     Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2)))/(2*np.pi))
# for a in range(13):
#     Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
# for a in range(13):
#     Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))

#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
# Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
# Planet_rotation_O4b = Planet_rotation_O4*1.1
Period_Uranus = (28559*2*np.pi)/(17.24*60*60)
Planet_rotation_O4 = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesS2.npy')
Approx_Longs = []

for o in range(25):
    Approx_Longs.append((Longs[o]+Longs[o+1])/2)

fig, ax = plt.subplots(figsize=(14,10))
ax2 = ax.twinx()
#ax.set_title('~1 hr exposure at ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=25, pad=25)
ax.set_title('b) iSHELL Intensity and LOS velocity', fontsize = 30, pad = 20)
for o in range(10):
    Planet_rotation_O4 = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-(o*10))))/(13.098297869601417)
    if o > 0:
        ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), lw = 2, color='k', ls = '--', alpha = 0.5)
    else:
        ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), lw = 3, color='g', ls = '--')
Planet_rotation_O4 = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), lw = 3, color='k', ls = '--', label = 'Planetary Rotation' )
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:25]), color='b', label='Average ion velocity', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
#ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:25]), np.flip(Velocity_Err_QAS_P[0:25]), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQAS[0:25]), color='r', ls= '--', label='IR Intensity', lw = 5, alpha = 0.5)
#ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N[0:25]), np.flip(INTS_Err_QAS_P[0:25]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
#ax.set_xticks(np.arange(2, 28, 5), labels=['50', '55', '60', '65', '70', '75'])
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('ORF $H_{3}^{+}$ Velocities (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower center', fontsize=25)
#ax2.legend(loc='upper center', fontsize=35)
ax.grid(linestyle = '--')
#plt.vlines(55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)


Ints_1158 = ex_INTSQAS
Errs_Ints_1158_N = INTS_Err_QAS_N
Errs_Ints_1158_P = INTS_Err_QAS_P

#%%
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

CleaningImage = CombinedImage[2]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
# CleaningImage[31:34,1209:1211] = np.nanmean(CombinedImage[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[3]
Wavelength_O4 = np.load('Wavelength_O4.npy')
wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225]))

dataI = (CleaningImage + CleaningImage2)/2
plt.figure()
plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
dataI[dataI < -0.004] = 0
for m in range(25): #31
    datai = dataI[7+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
    a1_pointP = a1_pointP[0][0]
    ex_resultQPAS = gmodel.fit(datai, x=np.arange(200), a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
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

for o in range(25):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(25):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[0:o+24]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)

Velocity_QAS = []
lamda_H3_O4 = 3.965
lamda_H3_O5 = 3.98558
c = 299792458 
for n in range(25):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(25):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(25):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(25):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

# Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

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
# for o in range(26):
#     Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2))/ (17.24*60*60)
#     Limb_velocity_pix_O4 = (Period_Uranus/1000)
#     Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

# #So the speed of Uranus at surface is 
# # Look into how to vary this across the slit so it can be used for 2016
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []

# Planet_diameter = []
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []
# for o in range(26):
#     Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2)))/(2*np.pi))
# for a in range(13):
#     Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
# for a in range(13):
#     Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))

#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
# Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
# Planet_rotation_O4b = Planet_rotation_O4*1.1
Period_Uranus = (28559*2*np.pi)/(17.24*60*60)
Planet_rotation_O4 = (np.arange(-12, 13, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesS2.npy')
Approx_Longs = []

for o in range(25):
    Approx_Longs.append((Longs[o]+Longs[o+1])/2)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:59' +' UTC on $11^{th}$ October 2016', fontsize=35, pad=25)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:25]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:25]), np.flip(Velocity_Err_QAS_P[0:25]), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQAS[0:25]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N[0:25]), np.flip(INTS_Err_QAS_P[0:25]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=35, labelpad=15)
ax.set_ylabel('ORF $H_{3}^{+}$ Velocities (kms$^{-1}$)', fontsize=35, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=35)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower center', fontsize=35)
ax2.legend(loc='upper center', fontsize=35)
ax.grid(linestyle = '--')
#plt.vlines(55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)


Ints_0859 = ex_INTSQAS
Errs_Ints_0859_N = INTS_Err_QAS_N
Errs_Ints_0859_P = INTS_Err_QAS_P

#%% Now lets do Ints and Temps/ CDs 
Temps_O859 = np.load('0859_Temps.npy')
Err_Temps_O859 = np.load('Err_0859_Temps.npy')
CDs_O859 = np.load('0859_CDs.npy')
Err_CDs_O859 = np.load('Err_0859_CDs.npy')
Temps_1158 = np.load('1158_Temps.npy')
Err_Temps_1158 = np.load('Err_1158_Temps.npy')
CDs_1158 = np.load('1158_CDs.npy')
Err_CDs_1158 = np.load('Err_1158_CDs.npy')

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:59' +' UTC on $11^{th}$ October 2016', fontsize=35, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:24], np.flip(Temps_O859), color='g', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.errorbar(Approx_Longs[1:24], np.flip(Temps_O859), yerr = np.flip(Err_Temps_O859), xerr = 0, color = 'g', ecolor='g', lw = 5, capthick = 2.5, label='$H_{3}^{+}$ Temperature')
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:24], np.flip(Temps_O859 + Err_Temps_O859), np.flip(Temps_O859 - Err_Temps_O859), color='g', alpha=0.5)
ax2.plot(Approx_Longs[1:24], np.flip(CDs_O859)/1e16, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
ax2.fill_between(Approx_Longs[1:24], np.flip(CDs_O859 + Err_CDs_O859)/1e16, np.flip(CDs_O859 - Err_CDs_O859)/1e16, color='purple', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=35, labelpad=15)
ax.set_ylabel('$H_{3}^{+}$ Temperature (K)', fontsize=35, labelpad=15)
ax2.set_ylabel(r'$H_{3}^{+}$ Column Density $(x 10^{16}m^{-2})$', fontsize=35, labelpad=15)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(350, 850)
ax2.set_ylim(0, 10)
ax.legend(loc='upper left', fontsize=35)
ax2.legend(loc='upper right', fontsize=35)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=35, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Temps_1158), color='g', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Temps_1158 + Err_Temps_1158), np.flip(Temps_1158 - Err_Temps_1158), color='g', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(CDs_1158)/1e16, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(CDs_1158 + Err_CDs_1158)/1e16, np.flip(CDs_1158 - Err_CDs_1158)/1e16, color='purple', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=30)
ax.tick_params(axis='both', which='major', labelsize=30)
ax2.tick_params(axis='both', which='major', labelsize=30)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=35, labelpad=15)
ax.set_ylabel('$H_{3}^{+}$ Temperature (K)', fontsize=35, labelpad=15)
ax2.set_ylabel(r'$H_{3}^{+}$ Column Density $(x 10^{16}m^{-2})$', fontsize=35, labelpad=15)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(350, 850)
ax2.set_ylim(0, 10)
ax.legend(loc='upper left', fontsize=35)
ax2.legend(loc='upper right', fontsize=35)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:24], np.flip(CDs_O859)/1e17, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:24], np.flip(CDs_O859 + Err_CDs_O859)/1e17, np.flip(CDs_O859 - Err_CDs_O859)/1e17, color='purple', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQAS[0:25]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N[0:25]), np.flip(INTS_Err_QAS_P[0:25]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel(r'Ionospheric Column Density of $H_{3}^{+}$ $(x 10^{17}m^{-2})$', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '11:57' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(CDs_1158)/1e17, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(CDs_1158 + Err_CDs_1158)/1e17, np.flip(CDs_1158 - Err_CDs_1158)/1e17, color='purple', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(Ints_1158), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(Errs_Ints_1158_N), np.flip(Errs_Ints_1158_P), color='r', alpha=0.5)
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
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

#%% Lets sort the Equator set and then North
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

from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

#Lets clean the image a bit more
CleaningImage = CombinedImage[0]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
# CleaningImage[31:34,1209:1211] = np.nanmean(CombinedImage[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[1]
Wavelength_O4 = np.load('Wavelength_O4.npy')
wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225]))

dataI = (CleaningImage + CleaningImage2)/2
plt.figure()
plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
dataI[dataI < -0.004] = 0
for m in range(27): #31
    datai = dataI[6+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
    a1_pointP = a1_pointP[0][0]
    ex_resultQPAS = gmodel.fit(datai, x=np.arange(200), a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
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

for o in range(27):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(27):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[0:o+26]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)

Velocity_QAS = []
lamda_H3_O4 = 3.965
lamda_H3_O5 = 3.98558
c = 299792458 
for n in range(27):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(27):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(27):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(27):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

# Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

# Limb_velocity_pix_O4_test = []
# for a in range(26):
#     aa = 12.5 - a
#     Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin(aa*np.pi/24))

# Planet_rotation_O4 = Limb_velocity_pix_O4_test[1:25]

R_Uranus = 28559*1000 #m (2*28559*np.pi*np.cos(np.deg2rad(Lats)))
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

Lats = np.load('LatitudesE.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
# for o in range(26):
#     Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2))/ (17.24*60*60)
#     Limb_velocity_pix_O4 = (Period_Uranus/1000)
#     Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

# #So the speed of Uranus at surface is 
# # Look into how to vary this across the slit so it can be used for 2016
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []

# Planet_diameter = []
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []
# for o in range(26):
#     Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2)))/(2*np.pi))
# for a in range(13):
#     Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
# for a in range(13):
#     Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))

#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
# Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
# Planet_rotation_O4b = Planet_rotation_O4*1.1
Period_Uranus = (28559*2*np.pi)/(17.24*60*60)
Planet_rotation_O4 = (np.arange(-13, 14, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesE.npy')
Approx_Longs = []

for o in range(27):
    Approx_Longs.append((Longs[o]))

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:00' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QAS), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N), np.flip(Velocity_Err_QAS_P), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQAS), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N), np.flip(INTS_Err_QAS_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('ORF $H_{3}^{+}$ Velocities (kms$^{-1}$)', fontsize=35, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower center', fontsize=25)
ax2.legend(loc='upper center', fontsize=25)
ax.grid(linestyle = '--')
#plt.vlines(55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)

Ints_0800 = ex_INTSQAS
Errs_Ints_0800_N = INTS_Err_QAS_N
Errs_Ints_0800_P = INTS_Err_QAS_P

Temps_0800 = np.load('Combo_Set1_Temps.npy')
Err_Temps_0800 = np.load('Combo_Set1_Temps_Err.npy')
CDs_0800 = np.load('Combo_Set1_CD.npy')
Err_CDs_0800 = np.load('Combo_Set1_CD_Err.npy')

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:00' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(CDs_0800)/1e17, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(CDs_0800 + Err_CDs_0800)/1e17, np.flip(CDs_0800 - Err_CDs_0800)/1e17, color='purple', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(Ints_0800), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N), np.flip(INTS_Err_QAS_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel(r'Ionospheric Column Density of $H_{3}^{+}$ $(x 10^{17}m^{-2})$', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '08:00' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(Temps_0800), color='g', ls= '--', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(Temps_0800 + Err_Temps_0800), np.flip(Temps_0800 - Err_Temps_0800), color='g', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(Ints_0800), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(Errs_Ints_0800_N), np.flip(Errs_Ints_0800_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('Ionospheric Temperature from $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(350, 850)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)
#%%
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

from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

#Lets clean the image a bit more
CleaningImage = CombinedImage[6]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
# CleaningImage[31:34,1209:1211] = np.nanmean(CombinedImage[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[7]
Wavelength_O4 = np.load('Wavelength_O4.npy')
wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225]))

dataI = (CleaningImage + CleaningImage2)/2
plt.figure()
plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
dataI[dataI < -0.004] = 0
for m in range(27): #31
    datai = dataI[5+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
    a1_pointP = a1_pointP[0][0]
    ex_resultQPAS = gmodel.fit(datai, x=np.arange(200), a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
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

for o in range(27):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(27):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[0:o+26]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)

Velocity_QAS = []
lamda_H3_O4 = 3.965
lamda_H3_O5 = 3.98558
c = 299792458 
for n in range(27):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(27):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(27):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(27):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

# Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

# Limb_velocity_pix_O4_test = []
# for a in range(26):
#     aa = 12.5 - a
#     Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin(aa*np.pi/24))

# Planet_rotation_O4 = Limb_velocity_pix_O4_test[1:25]

R_Uranus = 28559*1000 #m (2*28559*np.pi*np.cos(np.deg2rad(Lats)))
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

Lats = np.load('LatitudesE.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []


#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
# Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
# Planet_rotation_O4b = Planet_rotation_O4*1.1
Period_Uranus = (28559*2*np.pi)/(17.24*60*60)
Planet_rotation_O4 = (np.arange(-13, 14, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesE.npy')
Approx_Longs = []

for o in range(27):
    Approx_Longs.append((Longs[o]))

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs, np.flip(Velocity_QAS), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N), np.flip(Velocity_Err_QAS_P), color='b', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(ex_INTSQAS), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N), np.flip(INTS_Err_QAS_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower center', fontsize=25)
ax2.legend(loc='upper center', fontsize=25)
ax.grid(linestyle = '--')
#plt.vlines(55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)

Ints_1058 = ex_INTSQAS
Errs_Ints_1058_N = INTS_Err_QAS_N
Errs_Ints_1058_P = INTS_Err_QAS_P

Temps_1058 = np.load('Combo_Set4_Temps.npy')
Err_Temps_1058 = np.load('Combo_Set4_Temps_Err.npy')
CDs_1058 = np.load('Combo_Set4_CD.npy')
Err_CDs_1058 = np.load('Combo_Set4_CD_Err.npy')

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(CDs_1058)/1e17, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(CDs_1058 + Err_CDs_1058)/1e17, np.flip(CDs_1058 - Err_CDs_1058)/1e17, color='purple', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(Ints_1058), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N), np.flip(INTS_Err_QAS_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel(r'Ionospheric Column Density of $H_{3}^{+}$ $(x 10^{17}m^{-2})$', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '10:58' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(Temps_1058), color='g', ls= '--', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(Temps_1058 + Err_Temps_1058), np.flip(Temps_1058 - Err_Temps_1058), color='g', alpha=0.5)
ax2.plot(Approx_Longs, np.flip(Ints_1058), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs, np.flip(INTS_Err_QAS_N), np.flip(INTS_Err_QAS_P), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('Ionospheric Temperature from $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(350, 850)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

#%%
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

from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

#Lets clean the image a bit more
CleaningImage = CombinedImage[4]
# CleaningImage[CleaningImage == 0] = Mean/11
#CleaningImage[10:13,1209:1213] = np.nanmean(CombinedImage[8])
# CleaningImage[31:34,1209:1211] = np.nanmean(CombinedImage[8])
# CleaningImage[10:13,1209:1212] = np.nanmean(Q1_Data_Capture[8])
# CleaningImage[38:41,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[38,1209] = np.nanmean(CombinedImage[8])
# # CleaningImage[40,1209] = np.nanmean(CombinedImage[8])
# CleaningImage[39,1208] = np.nanmean(CombinedImage[8])

CleaningImage2 = CombinedImage[5]
Wavelength_O4 = np.load('Wavelength_O4.npy')
wav_pixel_ratio = np.nanmean(np.gradient(Wavelength_O4[1220:1225]))

dataI = (CleaningImage + CleaningImage2)/2
plt.figure()
plt.imshow(dataI, cmap='gist_gray', vmax = 0.01, vmin = -0.01)
dataI[dataI < -0.004] = 0
for m in range(27): #31
    datai = dataI[5+m, 1110:1310] #4
            #datai = datai+(acre(CombinedImage[8][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
            #datai = datai+(acre(CombinedImage[9][15+m, 1100:1300], width = 2, verbose = False)*AntiMask[15+m, 1100:1300])
        # datai = f.shift(datai, shift=-1*Offset_Recenter[24-m], mode='wrap')      
    gmodel = Model(gauss_fit)
    a1_pointP = np.where(datai == np.nanmax(datai[95:115]))
    a1_pointP = a1_pointP[0][0]
    ex_resultQPAS = gmodel.fit(datai, x=np.arange(200), a0=np.nanmax(datai[95:115]), a1=a1_pointP, a2=1.8, a3=0, a4=0, a5=0)
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

for o in range(27):
    QIntErr = ex_INTSQAS[o]*np.sqrt((ex_Err_A0_QAS[o]/ex_A0_QAS[o])**2 + (ex_Err_A2_QAS[o]/ex_A2_QAS[o])**2)
    QAS_IntErr.append(QIntErr)
  
    #%%
A1_QAS = []

o = 0
for n in range(27):
    A1 = ((ex_A1_QAS[n]-np.nanmean(ex_A1_QAS[0:o+27]))*-wav_pixel_ratio)
    #A1 = ((ex_A1_QAS[n]-1222))*-wav_pixel_ratio
    A1_QAS.append(A1)

Velocity_QAS = []
lamda_H3_O4 = 3.965
lamda_H3_O5 = 3.98558
c = 299792458 
for n in range(27):
    V = ((A1_QAS[n]/lamda_H3_O4)*c)
    Velocity_QAS.append(V/1000)
    
Velocity_QAS_Err = []

# Err1 = Uranus_width_pix_O4*(4.4111989937692037e-07/0.1818900559493812)
# Err3 = Uranus_width_pix_O5*(4.3670243056407377e-07/0.18006857011149846)

for n in range(27):
    Err = Velocity_QAS[n]*np.sqrt((ex_Err_A1_QAS[n]/A1_QAS[n])**2)
    Velocity_QAS_Err.append(np.sqrt(Err**2))
    
# Velocities_QAS = np.reshape(Velocity_QAS, (2, 25))
# Velocities_QAS = np.fliplr(np.flipud(Velocities_QAS))

Velocity_Err_QAS_N = []
Velocity_Err_QAS_P = []

for n in range(27):
    Velocity_Err_QAS_N.append(Velocity_QAS[n]-Velocity_QAS_Err[n])
    Velocity_Err_QAS_P.append(Velocity_QAS[n]+Velocity_QAS_Err[n])
    
INTS_Err_QAS_N = []
INTS_Err_QAS_P = []

for n in range(27):
    INTS_Err_QAS_N.append(ex_INTSQAS[n]-QAS_IntErr[n])
    INTS_Err_QAS_P.append(ex_INTSQAS[n]+QAS_IntErr[n])

# Limb_velocity_pix_O4_90 = Limb_velocity_pix_O4 + (0.1*Limb_velocity_pix_O4)

# Limb_velocity_pix_O4_test = []
# for a in range(26):
#     aa = 12.5 - a
#     Limb_velocity_pix_O4_test.append(Limb_velocity_pix_O4*np.sin(aa*np.pi/24))

# Planet_rotation_O4 = Limb_velocity_pix_O4_test[1:25]

R_Uranus = 28559*1000 #m (2*28559*np.pi*np.cos(np.deg2rad(Lats)))
Time_Uranus = 17.24 #hours
Time_s_Uranus = 17.24*60*60

Lats = np.load('LatitudesN.npy')
#Calculate the latitude for each pixel
Vels_Limb_Uranus = []
# for o in range(26):
#     Period_Uranus = 2*np.pi*R_Uranus*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2))/ (17.24*60*60)
#     Limb_velocity_pix_O4 = (Period_Uranus/1000)
#     Vels_Limb_Uranus.append(Limb_velocity_pix_O4)

# #So the speed of Uranus at surface is 
# # Look into how to vary this across the slit so it can be used for 2016
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []

# Planet_diameter = []
# Limb_velocity_pix_O4_test = []
# Limb_velocity_pix_O4_test2 = []
# for o in range(26):
#     Planet_diameter.append((2*np.pi*13.098297869601417*np.cos(np.deg2rad((Lats[o]+Lats[o+26])/2)))/(2*np.pi))
# for a in range(13):
#     Limb_velocity_pix_O4_test.append(Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))
# for a in range(13):
#     Limb_velocity_pix_O4_test2.append(-1*Vels_Limb_Uranus[a]*np.sin((12-a)*np.pi/(2*Planet_diameter[a])))

#%%
#Uranus diameter was 3.719" and the pixel scale is 4 = 0.1818900559493812 and 5 = 0.18006857011149846
Uranus_width_pix_O4 = 3.719/0.1818900559493812
Uranus_width_pix_O5 = 3.719/0.18006857011149846
  
# Planet_rotation_O4 = np.append(Limb_velocity_pix_O4_test, np.flip(Limb_velocity_pix_O4_test2[:-1]))
# Planet_rotation_O4b = Planet_rotation_O4*1.1
Period_Uranus = (28559*2*np.pi)/(17.24*60*60)
Planet_rotation_O4 = (np.arange(-13, 14, 1)*Period_Uranus*np.sin(np.deg2rad(90-35.78)))/(13.098297869601417)
# NORTH_SLIT = np.load('NORTH_SLIT_CONFIG.npy')
# EQUATOR_SLIT = np.load('EQUATOR_SLIT_CONFIG.npy')
Longs = np.load('LongitudesN.npy')
Approx_Longs = []

for o in range(27):
    Approx_Longs.append((Longs[o]))

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
ax.plot(Approx_Longs[1:-1], np.flip(Planet_rotation_O4[1:-1]), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(Velocity_QAS[1:-1]), color='b', ls= '--', label='Combined Q(1,0) and Q(3,0) IR Velocities', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(Velocity_Err_QAS_N[1:-1]), np.flip(Velocity_Err_QAS_P[1:-1]), color='b', alpha=0.5)
ax2.plot(Approx_Longs[1:-1], np.flip(ex_INTSQAS[1:-1]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:-1], np.flip(INTS_Err_QAS_N[1:-1]), np.flip(INTS_Err_QAS_P[1:-1]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('(ORF) LOS Velocities from $H_{3}^{+}$ (kms$^{-1}$)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax.set_ylim(-4, 4)
ax2.set_ylim(0, 0.75)
ax.legend(loc='lower center', fontsize=25)
ax2.legend(loc='upper center', fontsize=25)
ax.grid(linestyle = '--')
#plt.vlines(55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)
#plt.vlines(-55.00599839387289, -4, 4, colors='k', ls = '--', lw = 4, alpha = 0.5)

Ints_0959 = ex_INTSQAS
Errs_Ints_0959_N = INTS_Err_QAS_N
Errs_Ints_0959_P = INTS_Err_QAS_P

Temps_0959 = np.load('Combo_Set3_Temps.npy')
Err_Temps_0959 = np.load('Combo_Set3_Temps_Err.npy')
CDs_0959 = np.load('Combo_Set3_CD.npy')
Err_CDs_0959 = np.load('Combo_Set3_CD_Err.npy')

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(CDs_0959)/1e17, color='purple', ls= '--', label='$H_{3}^{+}$ Ion Density', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(CDs_0959 + Err_CDs_0959)/1e17, np.flip(CDs_0959 - Err_CDs_0959)/1e17, color='purple', alpha=0.5)
ax2.plot(Approx_Longs[1:-1], np.flip(ex_INTSQAS[1:-1]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:-1], np.flip(INTS_Err_QAS_N[1:-1]), np.flip(INTS_Err_QAS_P[1:-1]), color='r', alpha=0.5)
#ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
#ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
ax.tick_params(axis='both', which='major', labelsize=20)
ax2.tick_params(axis='both', which='major', labelsize=20)
# #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
ax.set_ylabel('Ionospheric Temperature from $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$)', fontsize=25)
#ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#ax.hlines(np.linspace(-4, 4, 17), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
ax.set_xlim(-90, 90) #0, 3601
ax2.set_ylim(0, 0.75)
ax.set_ylim(0, 1)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

fig, ax = plt.subplots(figsize=(10,8))
ax2 = ax.twinx()
ax.set_title('~1 hr exposure at ' + '09:59' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
#ax.plot(Approx_Longs, np.flip(Planet_rotation_O4), color='k', ls = '--', label='Planetary Rotation')
#ax.plot(Approx_Longs, Planet_rotation_O4b, color='green', ls = '--', label='Planetary Rotation + 250 m/s')
#ax.plot(Approx_Longitudes, Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
ax.plot(Approx_Longs[1:-1], np.flip(Temps_0959), color='g', ls= '--', label='$H_{3}^{+}$ Temperature', lw = 5)
#ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
#ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
ax.fill_between(Approx_Longs[1:-1], np.flip(Temps_0959 + Err_Temps_0959), np.flip(Temps_0959 - Err_Temps_0959), color='g', alpha=0.5)
ax2.plot(Approx_Longs[1:-1], np.flip(Ints_0959[1:-1]), color='r', label='IR Intensity', lw = 5)
ax2.fill_between(Approx_Longs[1:-1], np.flip(INTS_Err_QAS_N[1:-1]), np.flip(INTS_Err_QAS_P[1:-1]), color='r', alpha=0.5)
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
ax.set_ylim(350, 850)
ax.legend(loc='lower left', fontsize=25)
ax2.legend(loc='upper right', fontsize=25)
ax.grid(ls = 'dotted', alpha = 0.75)

#%% And now the Temps and Col Densities

