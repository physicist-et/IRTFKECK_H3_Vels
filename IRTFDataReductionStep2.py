# -*- coding: utf-8 -*-
"""
Created on Wed Jun 16 14:12:06 2021

@author: snowy
"""
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits #import the relevant directories to read a fits file from your directory and plot it
import h3ppy
from ACRE_tss import acre
h3p = h3ppy.h3p()

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
hdu = fits.open(image_file1)
hdr = hdu[0].header
#rint(hdr)

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
    
Total_IRTF_Data = []

for n in range(45):
    if n == 0:
        Total_IRTF_Data = All_IRTF_Data[0]
    else:
        Total_IRTF_Data = np.add(Total_IRTF_Data, All_IRTF_Data[n])

Total_IRTF_Data_Alt = Total_IRTF_Data/45
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt > 0.025] = 0
Total_IRTF_Data_Alt[Total_IRTF_Data_Alt < -0.025] = 0
T_IRTFData = np.flipud(Total_IRTF_Data_Alt)/2

#%% Now to get h3ppy involved (provide a model so we know what to compare to)
wave = h3p.wavegen((np.nanmin(Wavelength_O4))+0.0006052099999995455, (np.nanmax(Wavelength_O4)+0.0006052099999995455), 2048)
wave2 = h3p.wavegen((np.nanmin(Wavelength_O5)+0.0003660963848233223), (np.nanmax(Wavelength_O5)+0.00043096861236335826), 2048)
wave3 = h3p.wavegen((np.nanmin(Wavelength_O4)+0.0006052099999995455), (np.nanmax(Wavelength_O5)+0.00032096861236335826), 3932)
model = h3p.model(density = 8e15, temperature = 450, R = 75000, wavelength = wave)
model2 = h3p.model(density = 8e15, temperature = 450, R = 75000, wavelength = wave2)
model3 = h3p.model(density = 8e15, temperature = 450, R = 75000, wavelength = wave3)

noise = np.random.normal(size = model.size) * np.max(model) / 50
pretend_data = model + noise

noise2 = np.random.normal(size = model2.size) * np.max(model2) / 50
pretend_data2 = model2 + noise2

h3p.set(density = 8e15, temperature = 450, data = pretend_data, wavelength = wave)
# Let h3ppy make an initial guess at the density
h3p.guess_density()
params_to_fit = ['density', 'sigma_0']
# Fit temperature and density to the pretend data
fit = h3p.fit(params_to_fit)

h3p.set(density = 8e15, temperature = 450, data = pretend_data2, wavelength = wave2)
# Let h3ppy make an initial guess at the density
h3p.guess_density()
# Fit temperature and density to the pretend data
fit2 = h3p.fit()

fig, ax = plt.subplots()
ax.plot(wave, pretend_data * 1e3, 'o', label = 'IRTF iSHELL 2016 spectrum')
ax.plot(wave, fit * 1e3, label = 'h3ppy fit to data')
ax.legend(frameon = False)

fig, ax = plt.subplots()
ax.plot(wave, pretend_data2 * 1e3, 'o', label = 'IRTF iSHELL 2016 spectrum')
ax.plot(wave, fit2 * 1e3, label = 'h3ppy fit to data')
ax.legend(frameon = False)

#%% And now we start fitting
b = 0

psQ1 = 550 + 30
nsQ1 = 550 + 71
psQ3 = 440 + 30
nsQ3 = 440 + 71

# This function sub-divides data centered on a list of wavelengths
def subdivide(wave, spec, middles, width = 20) : 
    ret = []
    for m in middles : 
        centre = np.abs(wave - m).argmin()
        for i in range(centre - width, centre + width) : 
            ret.append(spec[i])
    return np.array(ret)

O1_IRTF_T = T_IRTFData[550:660,:]
O2_IRTF_T = T_IRTFData[440:550,:]
O3_IRTF_T = T_IRTFData[330:440,:]

P_O1_IRTF = O1_IRTF_T[29:53,:]
N_O1_IRTF = O1_IRTF_T[70:94,:]

Tot_O1_IRTF = np.zeros(2048)

P_O2_IRTF = O2_IRTF_T[29:53,:]
N_O2_IRTF = O2_IRTF_T[70:94,:]

Tot_O2_IRTF = np.zeros(2048)

for i in range(24):
    T_O1_IRTF = (P_O1_IRTF[i,:] - N_O1_IRTF[i,:])/2
    Tot_O1_IRTF = np.add(Tot_O1_IRTF, T_O1_IRTF)
    
Tot_O1_IRTF = Tot_O1_IRTF/24

for i in range(24):
    T_O2_IRTF = (P_O2_IRTF[i,:] - N_O2_IRTF[i,:])/2
    Tot_O2_IRTF = np.add(Tot_O2_IRTF, T_O2_IRTF)
    
Tot_O2_IRTF = Tot_O2_IRTF/24

# The H3+ line centeres contained withing this spectral band
centers = [3.952995962732753]
centers2 = [3.971064543098212]
centers3 = [3.985528767996527]
centers4 = [3.9870189233052002]
cpos = np.arange(4) * 41 + 20

Tot_O1_IRTF = np.concatenate((np.zeros(10), Tot_O1_IRTF))
Tot_O2_IRTF = np.concatenate((np.zeros(10), Tot_O2_IRTF))
Tot_O2_IRTF = np.concatenate((Tot_O2_IRTF, np.zeros(10)))

# Create sub-arrays, focusing on where the H3+ lines are
subspecQ1 = subdivide(wave, Tot_O1_IRTF[0:2048], centers)
subspecQ2 = subdivide(wave2, Tot_O2_IRTF[11:2059], centers2)
subspecQ3 = subdivide(wave2, Tot_O2_IRTF[8:2056], centers3)
subspecQ4 = subdivide(wave2, Tot_O2_IRTF[13:2061], centers4)
subwave = subdivide(wave, wave, centers)
subwave2 = subdivide(wave2, wave2, centers2)
subwave3 = subdivide(wave2, wave2, centers3)
subwave4 = subdivide(wave2, wave2, centers2)

subwave = np.concatenate((subwave, subwave2, subwave3, subwave4))
subspec = np.concatenate((subspecQ1, subspecQ2, subspecQ3, subspecQ4))

# Set the wavelength and the data
h3p.set(wavelength = subwave, data = subspec, R = 75000)

# Create a x scale for plotting 
xx      = range(len(subspec))

# Guess the density and proceed with a five parameter fit
h3p.guess_density()
fit = h3p.fit()
vars, errs = h3p.get_results()

centers = [3.952995962732753, 3.971064543098212, 3.985528767996527, 3.9870189233052002]
# Plot the fit
title = 'IRTF iSHELL H$_3^+$ spectrum of the total emission from Uranus October 2016'
fig, ax = plt.subplots()
ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'), xticks = cpos, title=title)
ax.set_xticklabels(centers)
ax.legend(frameon = False)
plt.tight_layout()

# =============================================================================
# for a in range(11):
#     data = np.flipud(All_IRTF_Data[b] + All_IRTF_Data[b+1] + All_IRTF_Data[b+2] + All_IRTF_Data[b+3])/4
#     data[data > 0.05] = 0
#     data[data < -0.05] = 0
#     data = acre(data, width = 5, verbose = False)
#     b += 4
#     for m in range(24):
#         datai = data[psQ1+m, :]
#         dataii = data[nsQ1+m, :]
#         dataiii = data[psQ3+m, :]
#         dataiv = data[nsQ3+m,:]
#         dataQ1 = (datai + dataii)/2
#         dataQ3 = (dataiii + dataiv)/2
#         # The H3+ line centeres contained withing this spectral band
#         centers = [3.952995962732753, 3.971064543098212, 3.985528767996527, 3.9870189233052002]
#         cpos = np.arange(4) * 41 + 20
# 
#         # Create sub-arrays, focusing on where the H3+ lines are
#         subspec = subdivide(wave, spec, centers)
#         subwave = subdivide(wave, wave, centers)
# 
#         # Set the wavelength and the data
#         h3p.set(wavelength = subwave, data = subspec, R = 20000)
# 
#         # Create a x scale for plotting 
#         xx      = range(len(subspec))
# 
#         # Guess the density and proceed with a five parameter fit
#         h3p.guess_density()
#         fit = h3p.fit()
#         vars, errs = h3p.get_results()
# 
#         # Plot the fit
#         fig, ax = plt.subplots()
#         ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
#         ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
#         ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'), xticks = cpos, title=title)
#         ax.set_xticklabels(centers)
#         ax.legend(frameon = False)
#         plt.tight_layout()
#         
# =============================================================================
