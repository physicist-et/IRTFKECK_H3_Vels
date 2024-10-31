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
    # if x == 33:
    #     IRTFData = acre(IRTFData, width = 15, verbose = False)
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
AB_Combined_Order4 = []
AB_Combined_Order5 = []
AB_Combined_Order5B = []

for x in range(45):
    IRTFData_Q1 = -1*Combined_Q1[x] + IRTF_Data_Q1[:,:,x]
    AB_Combined_Order4.append(IRTFData_Q1)

for x in range(45):
    IRTFData_Q3 = -1*Combined_Q3[x] + IRTF_Data_Q3[:,:,x]
    AB_Combined_Order5.append(IRTFData_Q3)
    
for x in range(45):
    IRTFData_Q31 = -1*Combined_Q31[x] + IRTF_Data_Q31[:,:,x]
    AB_Combined_Order5B.append(IRTFData_Q31)    
    
#%% Lets attempt to fit H3+ onto this (which didn' work so lets just focus Q1 and Q3 seperately so we use the Keck example from h3ppy)
from lmfit import Model
def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
    z = (x - a1) / a2
    y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
    return y

def subdivide(wave, spec, middles, width = 20) : 
    ret = []
    for m in middles : 
        centre = np.abs(wave - m).argmin()
        for i in range(centre - width, centre + width) : 
            ret.append(spec[i])
    return np.array(ret)

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

# centers = [3.953, 3.971, 3.986, 3.9945]
# cpos = np.arange(4) * 41 + 20

# # Create sub-arrays, focusing on where the H3+ lines are
# subspec = subdivide(wave, spec, centers)
# subwave = subdivide(wave, wave, centers)

# # Set the wavelength and the data
# h3p.set(wavelength = subwave, data = subspec, R = 20000)

# # Create a x scale for plotting 
# xx      = range(len(subspec))

# # Guess the density and proceed with a five parameter fit
# h3p.guess_density()
# fit = h3p.fit()
# vars, errs = h3p.get_results()

# # Plot the fit
# fig, ax = plt.subplots()
# ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
# ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
# ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'), xticks = cpos, title=title)
# ax.set_xticklabels(centers)
# ax.legend(frameon = False)
# plt.tight_layout()

Tots_Q1 = np.zeros((40, 2048))
Tots_Q3 = np.zeros((40, 2048))

for x in range(11):
    dataI = (AB_Combined_Order4[o]+AB_Combined_Order4[o+1]+AB_Combined_Order4[o+2]+AB_Combined_Order4[o+3])/8
    dataI[dataI < -0.004] = 0
    Q1_Data_Capture.append(dataI[49:89, :])
    Tots_Q1 += dataI[49:89, :]
    dataI = (AB_Combined_Order5[o]+AB_Combined_Order5[o+1]+AB_Combined_Order5[o+2]+AB_Combined_Order5[o+3])/8
    dataI[dataI < -0.004] = 0
    Q3_Data_Capture.append(dataI[49:89, :])
    Tots_Q3 += dataI[49:89, :]
    o += 4

# #%% Now conduct h3ppy on this
# center1 = [3.95295]
# center2 = [3.98558]
# #cpos = np.arange(2) * 41 + 20
# CleaningImage = Q1_Data_Capture[8]
# CleaningImage[31:34,1209:1211] = np.nanmean(Q1_Data_Capture[8])

# CleaningImage2 = Q1_Data_Capture[9]
# CleaningImage3 = Q3_Data_Capture[8]
# CleaningImage3[6,1325:1355] = np.nanmean(Q3_Data_Capture[8])
# CleaningImage4 = Q3_Data_Capture[9]
# CleaningImage4[6,1325:1355] = np.nanmean(Q3_Data_Capture[9])
wave = np.load('Wavelength_O4.npy')
wave2 = np.load('Wavelength_O5.npy')
# Temps_AURORA = []
# Err_Temps_AURORA = []
# CD_AURORA = []
# Err_CD_AURORA = []
# Tempo_Temp = np.load('Combo_Temps_Set5.npy')
# Tempo_Density = np.load('Combo_CD_Set5.npy')

# for x in range(1):
#     dataI = (CleaningImage + CleaningImage2)/2
#     dataII = (CleaningImage3 + CleaningImage4)/2
#     subwaveQ1 = wave[1197:1238]
#     subwaveQ3 = wave2[1342:1383]
#     subwave = np.concatenate((subwaveQ1, subwaveQ3))
#     for m in range(25):
#         datai = dataI[6+m, :] 
#         data2i = dataI[6+m, :] #datai = dataI[8+m, :] 
#         dataii= dataII[6+m, :] 
#         data2ii = dataII[6+m, :]
#         #Now to centraalise the data
#         a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ1 = datai[a1_pointP-20:a1_pointP+21]
#         a1_pointP = np.where(data2i == np.nanmax(data2i[1205:1225]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ1a = data2i[a1_pointP-20:a1_pointP+21]
#         a1_pointP = np.where(dataii == np.nanmax(dataii[1346:1367]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ3 = dataii[a1_pointP-20:a1_pointP+21]
#         a1_pointP = np.where(data2ii == np.nanmax(data2ii[1346:1367]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ3a = data2ii[a1_pointP-20:a1_pointP+21]
#         subspecQ1 = (subspecQ1 + subspecQ1a)/2
#         subspecQ3 = (subspecQ3 + subspecQ3a)/2
#         subspec = np.concatenate((subspecQ1, subspecQ3))
#         if m == 0:
#             h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 800, density = 3*10**15)
#         elif m == 23 or 24:
#             h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 500, density = 7*10**15)
#         else:
#             h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = Tempo_Temp[m], density = Tempo_Density[m])
#         # Guess the density and proceed with a five parameter fit
#         h3p.guess_density(verbose=False) #Keep an append of the Temperaturs and Density so they can be used to guess in the next fittings
#         h3p.set(nbackground = 1)
#         ptf = ['background_0']
#         fit = h3p.fit(params_to_fit = ptf)
#         Vars, errs = h3p.get_results(verbose=False)
#         h3p.set(noffset = 2)
#         ptf = ['offset_0', 'offset_1']
#         fit2 = h3p.fit(params_to_fit = ptf, verbose = False)           
#         Vars, errs = h3p.get_results(verbose=False)
#         ptf = ['temperature', 'density', 'sigma_0']
#         fit3 = h3p.fit(params_to_fit = ptf)
#         Vars, errs = h3p.get_results()
#         Temps_AURORA.append(Vars['temperature'])
#         Err_Temps_AURORA.append(errs['temperature'])
#         CD_AURORA.append(Vars['density'])
#         Err_CD_AURORA.append(errs['density'])
#         xx = range(len(subspec))
#         # # Plot the fit
#         # fig, ax = plt.subplots()
#         # ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
#         # ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
#         # ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'))
#         # #ax.set_xticklabels(np.concatenate((center1, center2)))
#         # ax.legend(frameon = False)
#         # plt.tight_layout()

# plt.figure()
# plt.plot(np.arange(25), Temps_AURORA)
# plt.errorbar(np.arange(25), Temps_AURORA, yerr=Err_Temps_AURORA)

# plt.figure()
# plt.plot(np.arange(25), CD_AURORA)
# plt.errorbar(np.arange(25), CD_AURORA, yerr=Err_CD_AURORA)
# plt.ylim(0, 1e17)

# Aurora_Temps = []
# Aurora_CDs = []
# Aurora_Temps.append(Temps_AURORA)
# Aurora_CDs.append(CD_AURORA)

# Err_Aurora_Temps = []
# Err_Aurora_CDs = []
# Err_Aurora_Temps.append(Err_Temps_AURORA)
# Err_Aurora_CDs.append(Err_CD_AURORA)
o = 0
for x in range(11):
    dataI = (AB_Combined_Order4[o]+AB_Combined_Order4[o+1]+AB_Combined_Order4[o+2]+AB_Combined_Order4[o+3])/8
    dataI[dataI < -0.004] = 0
    Q1_Data_Capture.append(dataI[49:89, :])
    Tots_Q1 += dataI[49:89, :]
    dataI = (AB_Combined_Order5[o]+AB_Combined_Order5[o+1]+AB_Combined_Order5[o+2]+AB_Combined_Order5[o+3])/8
    dataI[dataI < -0.004] = 0
    Q3_Data_Capture.append(dataI[49:89, :])
    Tots_Q3 += dataI[49:89, :]
    o += 4

#%% Now Set 4
Temps_AURORA = []
Err_Temps_AURORA = []
CD_AURORA = []
Err_CD_AURORA = []
dgsgdhgfjTempo_Temp = np.load('Temps_Set4.npy')
Tempo_Density = np.load('CD_Set4.npy')
#%% Now conduct h3ppy on this
center1 = [3.95295]
center2 = [3.98558]
#cpos = np.arange(2) * 41 + 20

for x in range(1):
    dataI = (Q1_Data_Capture[6] + Q1_Data_Capture[7])/2
    dataII = (Q3_Data_Capture[6] + Q3_Data_Capture[7])/2 #Widen up the fittings to 20 and 21 for better fitting
    # dataI[28:31, 1191:1194] = np.nanmean(dataI)
    subwaveQ1 = wave[1177:1258]
    subwaveQ3 = wave2[1312:1393]
    for m in range(25):
        datai = dataI[6+m, :]
        dataii= dataII[6+m, :]       
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0]
        subspecQ1 = datai[a1_pointP-40:a1_pointP+41]
        a1_pointP = np.where(dataii == np.nanmax(dataii[1346:1367]))
        a1_pointP = a1_pointP[0][0]
        subspecQ3 = dataii[a1_pointP-40:a1_pointP+41]
        subspec = np.concatenate((subspecQ1, subspecQ3))
        subwave = np.concatenate((subwaveQ1, subwaveQ3))
        h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 445, density = 2*10**16)
        model = h3p.model(density = 2e16, temperature  = 445, R = 64000, wavelength = subwave)
        if m == 0:
            h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 400, density = 10**16)
            model = h3p.model(density = 1e16, temperature  = 400, R = 64000, wavelength = subwave)
        if m == 9 or m == 18 or m == 19:
            h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 500, density = 10**16)
            model = h3p.model(density = 1e16, temperature  = 500, R = 64000, wavelength = subwave)
        if m == 24 or m == 1:
            h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 700, density = 1.5*10**16)
            model = h3p.model(density = 1.5e16, temperature  = 700, R = 64000, wavelength = subwave)
        # Guess the density and proceed with a five parameter fit
        h3p.guess_density(verbose=False) #Keep an append of the Temperaturs and Density so they can be used to guess in the next fittings
        h3p.set(nbackground = 1)
        ptf = ['background_0']
        fit = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results(verbose=False)
        h3p.set(noffset = 2)
        ptf = ['offset_0', 'offset_1']
        fit2 = h3p.fit(params_to_fit = ptf, verbose = False)           
        Vars, errs = h3p.get_results(verbose=False)
        ptf = ['temperature', 'density', 'sigma_0']
        fit3 = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results()
        Temps_AURORA.append(Vars['temperature'])
        Err_Temps_AURORA.append(errs['temperature'])
        CD_AURORA.append(Vars['density'])
        Err_CD_AURORA.append(errs['density'])
        xx = range(len(subspec))
        # Plot the fit
        fig, ax = plt.subplots()
        ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
        ax.plot(xx, model * 1e3, 'o-', label = 'Model')
        # ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
        ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'))
        #ax.set_xticklabels(np.concatenate((center1, center2)))
        ax.legend(frameon = False)
        plt.tight_layout()

#%%
plt.figure()
plt.plot(np.arange(25), Temps_AURORA)
plt.errorbar(np.arange(25), Temps_AURORA, yerr=Err_Temps_AURORA)

plt.figure()
plt.plot(np.arange(25), CD_AURORA)
plt.errorbar(np.arange(25), CD_AURORA, yerr=Err_CD_AURORA)
plt.ylim(0, 1e17)

#%% Set 3 #Finish this bits
Temps_AURORA = []
Err_Temps_AURORA = []
CD_AURORA = []
Err_CD_AURORA = []
Tempo_Temp = np.load('Temps_Set3.npy')
Tempo_Density = np.load('CD_Set3.npy')

for x in range(1):
    dataI = (Q1_Data_Capture[4] + Q1_Data_Capture[5])/2
    dataII = (Q3_Data_Capture[4] + Q3_Data_Capture[5])/2
    for m in range(25):
        datai = dataI[6+m, :]
        dataii= dataII[6+m, :]       
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0]
        subspecQ1 = datai[a1_pointP-40:a1_pointP+41]
        a1_pointP = np.where(dataii == np.nanmax(dataii[1346:1367]))
        a1_pointP = a1_pointP[0][0]
        subspecQ3 = dataii[a1_pointP-40:a1_pointP+41]
        subspec = np.concatenate((subspecQ1, subspecQ3))
        h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = Tempo_Temp[m-1], density = Tempo_Density[m-1])
        if m == 0:
             h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 700, density = 5*10**15)
        # Guess the density and proceed with a five parameter fit
        h3p.guess_density(verbose=False) #Keep an append of the Temperaturs and Density so they can be used to guess in the next fittings
        h3p.set(nbackground = 1)
        ptf = ['background_0']
        fit = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results(verbose=False)
        h3p.set(noffset = 2)
        ptf = ['offset_0', 'offset_1']
        fit2 = h3p.fit(params_to_fit = ptf, verbose = False)           
        Vars, errs = h3p.get_results(verbose=False)
        ptf = ['temperature', 'density', 'sigma_0']
        fit3 = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results()
        Temps_AURORA.append(Vars['temperature'])
        Err_Temps_AURORA.append(errs['temperature'])
        CD_AURORA.append(Vars['density'])
        Err_CD_AURORA.append(errs['density'])
        xx = range(len(subspec))
        # Plot the fit
        # fig, ax = plt.subplots()
        # ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
        # ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
        # ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'))
        # #ax.set_xticklabels(np.concatenate((center1, center2)))
        # ax.legend(frameon = False)
        # plt.tight_layout()

plt.figure()
plt.plot(np.arange(25), Temps_AURORA)
plt.errorbar(np.arange(25), Temps_AURORA, yerr=Err_Temps_AURORA)

plt.figure()
plt.plot(np.arange(25), CD_AURORA)
plt.errorbar(np.arange(25), CD_AURORA, yerr=Err_CD_AURORA)
plt.ylim(0, 1e17)

#%%Set 2
Temps_AURORA = []
Err_Temps_AURORA = []
CD_AURORA = []
Err_CD_AURORA = []
Tempo_Temp = np.load('Combo_Temps_Set2.npy')
Tempo_Density = np.load('Combo_CD_Set2.npy')

for x in range(1):
    dataI = (Q1_Data_Capture[2] + Q1_Data_Capture[3])/2
    dataII = (Q3_Data_Capture[2] + Q3_Data_Capture[3])/2
    #dataI[8:15,1190:1193] = np.nanmean(dataI)
    for m in range(23):
        datai = dataI[8+m, :]
        data2i = dataI[8+m, :]
        dataii = dataII[8+m, :]
        data2ii = dataII[8+m, :]
        a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
        a1_pointP = a1_pointP[0][0]
        subspecQ1 = datai[a1_pointP-40:a1_pointP+41]
        a1_pointP = np.where(data2i == np.nanmax(data2i[1205:1225]))
        a1_pointP = a1_pointP[0][0]
        subspecQ1a = data2i[a1_pointP-40:a1_pointP+41]
        a1_pointP = np.where(dataii == np.nanmax(dataii[1346:1367]))
        a1_pointP = a1_pointP[0][0]
        subspecQ3 = dataii[a1_pointP-40:a1_pointP+41]
        a1_pointP = np.where(data2ii == np.nanmax(data2ii[1346:1367]))
        a1_pointP = a1_pointP[0][0]
        subspecQ3a = data2ii[a1_pointP-40:a1_pointP+41]
        subspecQ1 = (subspecQ1 + subspecQ1a)/2
        subspecQ3 = (subspecQ3 + subspecQ3a)/2
        subspec = np.concatenate((subspecQ1, subspecQ3))
        h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = Tempo_Temp[m-1], density = Tempo_Density[m-1])
        if m == 0 or m == 22:
            h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 700, density = 5*10**15)
        # Guess the density and proceed with a five parameter fit
        h3p.guess_density(verbose=False) #Keep an append of the Temperaturs and Density so they can be used to guess in the next fittings
        h3p.set(nbackground = 1)
        ptf = ['background_0']
        fit = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results(verbose=False)
        h3p.set(noffset = 2)
        ptf = ['offset_0', 'offset_1']
        fit2 = h3p.fit(params_to_fit = ptf, verbose = False)           
        Vars, errs = h3p.get_results(verbose=False)
        ptf = ['temperature', 'density', 'sigma_0']
        fit3 = h3p.fit(params_to_fit = ptf)
        Vars, errs = h3p.get_results()
        Temps_AURORA.append(Vars['temperature'])
        Err_Temps_AURORA.append(errs['temperature'])
        CD_AURORA.append(Vars['density'])
        Err_CD_AURORA.append(errs['density'])
        xx = range(len(subspec))
        # Plot the fit
        # fig, ax = plt.subplots()
        # ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
        # ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
        # ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'))
        # #ax.set_xticklabels(np.concatenate((center1, center2)))
        # ax.legend(frameon = False)
        # plt.tight_layout()

plt.figure()
plt.plot(np.arange(23), Temps_AURORA)
plt.errorbar(np.arange(23), Temps_AURORA, yerr=Err_Temps_AURORA)

plt.figure()
plt.plot(np.arange(23), CD_AURORA)
plt.errorbar(np.arange(23), CD_AURORA, yerr=Err_CD_AURORA)
plt.ylim(0, 1e17)
# #%%
# from scipy.optimize import curve_fit
# Longs = np.load('LongitudesS.npy')
# South_Pos = []
# for o in range(25):
#     South_Pos.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)

# New_Approx_Longs = South_Pos[1:24]

# Perr_Temps_AURORA = []
# Nerr_Temps_AURORA = []
# Perr_CD_AURORA = []
# Nerr_CD_AURORA= []

# for n in range(23):
#     Perr_Temps_AURORA.append(Temps_AURORA[n] + Err_Temps_AURORA[n])
#     Nerr_Temps_AURORA.append(Temps_AURORA[n] - Err_Temps_AURORA[n])
#     Perr_CD_AURORA.append(CD_AURORA[n] + Err_CD_AURORA[n])
#     Nerr_CD_AURORA.append(CD_AURORA[n] - Err_CD_AURORA[n])
    
# Perr_Temps_AURORA1 = []
# Nerr_Temps_AURORA1 = []
# Perr_CD_AURORA1 = []
# Nerr_CD_AURORA1= []

# for n in range(25):
#     Perr_Temps_AURORA1.append(Aurora_Temps[0][n] + Err_Aurora_Temps[0][n])
#     Nerr_Temps_AURORA1.append(Aurora_Temps[0][n] - Err_Aurora_Temps[0][n])
#     Perr_CD_AURORA1.append(Aurora_CDs[0][n] + Err_Aurora_CDs[0][n])
#     Nerr_CD_AURORA1.append(Aurora_CDs[0][n] - Err_Aurora_CDs[0][n])
    
# fig, ax = plt.subplots(figsize=(10,8))
# #ax2 = ax.twinx()
# # ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(New_Approx_Longs[0:23], np.flip(Temps_AURORA), color='b', label='8:59 UTC')
# ax.fill_between(New_Approx_Longs[0:23], np.flip(Perr_Temps_AURORA), np.flip(Nerr_Temps_AURORA), color='b', alpha=0.5)
# ax.plot(South_Pos, np.flip(Aurora_Temps[0]), color='r', label='11:57 UTC')
# ax.fill_between(South_Pos, np.flip(Perr_Temps_AURORA1), np.flip(Nerr_Temps_AURORA1), color='r', alpha=0.5)
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# #ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# #ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# #ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('Ionospheric Temperature from $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
# #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# #ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(350, 750)
# #ax2.set_ylim(0, 0.75)
# ax.legend(loc='upper right', fontsize=25)
# ax.grid(alpha = 0.5)
# #ax2.legend(loc='upper right', fontsize=25) 

# fig, ax = plt.subplots(figsize=(10,8))
# #ax2 = ax.twinx()
# # ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(New_Approx_Longs[0:23], np.flip(CD_AURORA)/1e16, color='b', label='8:59 UTC')
# ax.fill_between(New_Approx_Longs[0:23], np.flip(Perr_CD_AURORA)/1e16, np.flip(Nerr_CD_AURORA)/1e16, color='b', alpha=0.5)
# ax.plot(South_Pos, np.flip(Aurora_CDs[0])/1e16, color='r', label='11:57 UTC')
# ax.fill_between(South_Pos, np.flip(Perr_CD_AURORA1)/1e16, np.flip(Nerr_CD_AURORA1)/1e16, color='r', alpha=0.5)
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# #ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# #ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# #ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('Ionospheric Column Density of $H_{3}^{+}$ (x 10$^{16}$ m$^{-2}$)', fontsize=25, labelpad=15)
# #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# #ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(0.1, 9)
# #ax2.set_ylim(0, 0.75)
# ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ax.grid(alpha = 0.7, linestyle = '--')
# #ax2.legend(loc='upper right', fontsize=25) 
# #%%
# Temps_AURORA = []
# Err_Temps_AURORA = []
# CD_AURORA = []
# Err_CD_AURORA = []
# Tempo_Temp = np.load('Combo_Temps_Set1.npy')
# Tempo_Density = np.load('Combo_CD_Set1.npy')

# for x in range(1):
#     dataI = (Q1_Data_Capture[0] + Q1_Data_Capture[1])/2
#     dataII = (Q3_Data_Capture[0] + Q3_Data_Capture[1])/2
#     #dataI[8:15,1190:1193] = np.nanmean(dataI)
#     for m in range(25):
#         datai = dataI[6+m, :]
#         data2i = dataI[6+m, :]
#         dataii = dataII[6+m, :]
#         data2ii = dataII[6+m, :]
#         a1_pointP = np.where(datai == np.nanmax(datai[1205:1225]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ1 = datai[a1_pointP-40:a1_pointP+41]
#         a1_pointP = np.where(data2i == np.nanmax(data2i[1205:1225]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ1a = data2i[a1_pointP-40:a1_pointP+41]
#         a1_pointP = np.where(dataii == np.nanmax(dataii[1346:1367]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ3 = dataii[a1_pointP-40:a1_pointP+41]
#         a1_pointP = np.where(data2ii == np.nanmax(data2ii[1346:1367]))
#         a1_pointP = a1_pointP[0][0]
#         subspecQ3a = data2ii[a1_pointP-40:a1_pointP+41]
#         subspecQ1 = (subspecQ1 + subspecQ1a)/2
#         subspecQ3 = (subspecQ3 + subspecQ3a)/2
#         subspec = np.concatenate((subspecQ1, subspecQ3))
#         h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = Tempo_Temp[m-1], density = Tempo_Density[m-1])
#         if m == 0 or m == 23 or m == 24:
#             h3p.set(wavelength = subwave, data = subspec, R = 64000, temperature = 700, density = 5*10**15)
#         # Guess the density and proceed with a five parameter fit
#         h3p.guess_density(verbose=False) #Keep an append of the Temperaturs and Density so they can be used to guess in the next fittings
#         h3p.set(nbackground = 1)
#         ptf = ['background_0']
#         fit = h3p.fit(params_to_fit = ptf)
#         Vars, errs = h3p.get_results(verbose=False)
#         h3p.set(noffset = 2)
#         ptf = ['offset_0', 'offset_1']
#         fit2 = h3p.fit(params_to_fit = ptf, verbose = False)           
#         Vars, errs = h3p.get_results(verbose=False)
#         ptf = ['temperature', 'density', 'sigma_0']
#         fit3 = h3p.fit(params_to_fit = ptf)
#         Vars, errs = h3p.get_results()
#         Temps_AURORA.append(Vars['temperature'])
#         Err_Temps_AURORA.append(errs['temperature'])
#         CD_AURORA.append(Vars['density'])
#         Err_CD_AURORA.append(errs['density'])
#         xx = range(len(subspec))
#         # Plot the fit
#         # fig, ax = plt.subplots()
#         # ax.plot(xx, subspec * 1e3, '.', label = 'Observation')
#         # ax.plot(xx, fit * 1e3, label = 'h3ppy H$_3^+$ fit')
#         # ax.set(xlabel = h3p.xlabel(), ylabel = h3p.ylabel(prefix = 'm'))
#         # #ax.set_xticklabels(np.concatenate((center1, center2)))
#         # ax.legend(frameon = False)
#         # plt.tight_layout()

# plt.figure()
# plt.plot(np.arange(26), Temps_AURORA)
# plt.errorbar(np.arange(26), Temps_AURORA, yerr=Err_Temps_AURORA)

# plt.figure()
# plt.plot(np.arange(26), CD_AURORA)
# plt.errorbar(np.arange(26), CD_AURORA, yerr=Err_CD_AURORA)
# plt.ylim(0, 1e17)

# #%%
# from scipy.optimize import curve_fit
# Longs = np.load('LongitudesS.npy')
# South_Pos = []
# for o in range(25):
#     South_Pos.append((Longs[o]+Longs[o+1]+Longs[o+26]+Longs[o+27])/4)

# New_Approx_Longs = South_Pos

# Perr_Temps_AURORA = []
# Nerr_Temps_AURORA = []
# Perr_CD_AURORA = []
# Nerr_CD_AURORA= []

# for n in range(25):
#     Perr_Temps_AURORA.append(Temps_AURORA[n] + Err_Temps_AURORA[n])
#     Nerr_Temps_AURORA.append(Temps_AURORA[n] - Err_Temps_AURORA[n])
#     Perr_CD_AURORA.append(CD_AURORA[n] + Err_CD_AURORA[n])
#     Nerr_CD_AURORA.append(CD_AURORA[n] - Err_CD_AURORA[n])
    
# Perr_Temps_AURORA1 = []
# Nerr_Temps_AURORA1 = []
# Perr_CD_AURORA1 = []
# Nerr_CD_AURORA1= []

# for n in range(25):
#     Perr_Temps_AURORA1.append(Aurora_Temps[0][n] + Err_Aurora_Temps[0][n])
#     Nerr_Temps_AURORA1.append(Aurora_Temps[0][n] - Err_Aurora_Temps[0][n])
#     Perr_CD_AURORA1.append(Aurora_CDs[0][n] + Err_Aurora_CDs[0][n])
#     Nerr_CD_AURORA1.append(Aurora_CDs[0][n] - Err_Aurora_CDs[0][n])
    
# fig, ax = plt.subplots(figsize=(10,8))
# #ax2 = ax.twinx()
# # ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(New_Approx_Longs, np.flip(Temps_AURORA), color='b', label='8:59 UTC')
# ax.fill_between(New_Approx_Longs, np.flip(Perr_Temps_AURORA), np.flip(Nerr_Temps_AURORA), color='b', alpha=0.5)
# ax.plot(New_Approx_Longs, np.flip(Aurora_Temps[0]), color='r', label='11:57 UTC')
# ax.fill_between(New_Approx_Longs, np.flip(Perr_Temps_AURORA1), np.flip(Nerr_Temps_AURORA1[0:23]), color='r', alpha=0.5)
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# #ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# #ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# #ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('Ionospheric Temperature from $H_{3}^{+}$ (K)', fontsize=25, labelpad=15)
# #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# #ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(350, 750)
# #ax2.set_ylim(0, 0.75)
# ax.legend(loc='upper right', fontsize=25)
# ax.grid(alpha = 0.5)
# #ax2.legend(loc='upper right', fontsize=25) 

# fig, ax = plt.subplots(figsize=(10,8))
# #ax2 = ax.twinx()
# # ax.set_title('~1 hr exposure around ' + '12:11' +' UTC on $11^{th}$ October 2016', fontsize=22, pad=25)
# ax.plot(New_Approx_Longs, np.flip(CD_AURORA)/1e16, color='b', label='8:59 UTC')
# ax.fill_between(New_Approx_Longs, np.flip(Perr_CD_AURORA)/1e16, np.flip(Nerr_CD_AURORA)/1e16, color='b', alpha=0.5)
# ax.plot(New_Approx_Longs, np.flip(Aurora_CDs[0])/1e16, color='r', label='11:57 UTC')
# ax.fill_between(New_Approx_Longs, np.flip(Perr_CD_AURORA1)/1e16, np.flip(Nerr_CD_AURORA1[0:23])/1e16, color='r', alpha=0.5)
# #ax.plot(np.arange(26), Planet_rotation_O4_90, color='green', ls = '--', label='Planetary Rotation (+250m/s)')
# #ax.plot(Approx_Longs, np.flip(Velocity_QAS[0:23]), color='b', ls= '--', label='Combined Q(1,0) and (3,0) IR Velocities', lw = 5)
# #ax.plot(np.arange(25), Velocity_O58 - Planet_rotation_O58, color='c', label='Q(3,0) IR Velocities')
# #ax.fill_between(np.arange(25), (Velocities_QAS[0,:] - Velocity_Err_QAS[0,:]), (Velocities_QAS[0,:] + Velocity_Err_QAS[0,:]), color='b', alpha=0.5)
# #ax.fill_between(Approx_Longs, np.flip(Velocity_Err_QAS_N[0:23]), np.flip(Velocity_Err_QAS_P[24:47]), color='b', alpha=0.5)
# #ax2.plot(np.arange(25), np.flip(ex_INTSQ38), color='m', label='Q(3,0) Intensity', lw = 5)
# #ax.xticks(np.arange(0, 181, 30), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 181, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$'], fontsize=15)
# # #plt.xticks(np.arange(0, 361, 30), ['0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$', '120$^\circ$', '150$^\circ$', '180$^\circ$', '210$^\circ$', '240$^\circ$', '270$^\circ$', '300$^\circ$', '330$^\circ$', '360$^\circ$'], fontsize=15)
# ax.set_xticks(np.linspace(-90, 90, 7), labels=['-90$^\circ$', '-60$^\circ$', '-30$^\circ$', '0$^\circ$', '30$^\circ$', '60$^\circ$', '90$^\circ$'], fontsize=20)
# ax.tick_params(axis='both', which='major', labelsize=20)
# #ax2.tick_params(axis='both', which='major', labelsize=20)
# # #plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
# ax.set_xlabel('Arbitrary Uranian Longitude ($^\circ$)', fontsize=25, labelpad=15)
# ax.set_ylabel('Ionospheric Column Density of $H_{3}^{+}$ (x 10$^{16}$ m$^{-2}$)', fontsize=25, labelpad=15)
# #ax2.set_ylabel(r'Intensity (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' ' + str("{:.2f}".format(round(np.nanmean(QAS_IntErr[28:51]), 2))) + ' ${\mu}Wm^{-2}sr^{-1}$', fontsize=25)
# #ax.vlines(np.linspace(-90, 90, 7), -50, 50, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# #ax.hlines(np.linspace(-3, 3, 7), -90, 90, colors = 'k', linestyles = 'dotted', alpha = 0.75)
# ax.set_xlim(-90, 90) #0, 3601
# ax.set_ylim(0.1, 7)
# #ax2.set_ylim(0, 0.75)
# ax.legend(loc='upper right', fontsize=25) #Correct for longitude!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# ax.grid(alpha = 0.7, linestyle = '--')
# #ax2.legend(loc='upper right', fontsize=25) 
