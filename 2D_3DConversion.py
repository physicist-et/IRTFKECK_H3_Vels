# -*- coding: utf-8 -*-
"""
Created on Fri Jun 18 16:17:50 2021

@author: snowy
"""
import numpy as np
import math
# from IRTFIntensityStep1 import All_IRTF_Data, Mapping_Q1A
import matplotlib.pyplot as plt
import csv
from matplotlib.patches import Ellipse

#%%
#11th October starts from 134.926872 so 0 will be this
#Testing script for how to get Uranus at Keck 2006 and IRTF 2016
uranus_seangleK = (-5.075900 + -5.065686)/2
uranus_seangleI = -27.342
#stretch yy to become a sphere
flattening =0.0229
uranus_seangleK_Rads = (uranus_seangleK*np.pi)/180
uranus_seangleI_Rads = (uranus_seangleI*np.pi)/180
losflattening=flattening*(1-np.sin(uranus_seangleK_Rads))
eq_po_ratioK=1-losflattening

y_Keck = np.flip(np.linspace(-16, 16, num=33))

losflattening=flattening*(1-np.sin((uranus_seangleI/180)*np.pi))
eq_po_ratioI=1-losflattening

y1_IRTF = 0.6479317709801423*2
y2_IRTF = -0.6479317709801423*2

y1_IRTFadj = y1_IRTF
y2_IRTFadj = y2_IRTF

#%% Now that Keck and IRTF are adjusted with y for flattening we can then start mapping latitude and longitude
#First Keck

x_IRTF = np.flip(np.linspace(-16, 16, num=33))
r1_IRTF = np.sqrt(x_IRTF[0]**2 + y1_IRTFadj**2)
r11_IRTF = np.sqrt(x_IRTF[-1]**2 + y1_IRTFadj**2)
r2_IRTF = np.sqrt(x_IRTF[0]**2 + y2_IRTFadj**2)
r21_IRTF = np.sqrt(x_IRTF[-1]**2 + y2_IRTFadj**2)
print(r1_IRTF)
print(r2_IRTF)

pheta1_IRTF = math.asin(r1_IRTF/16.380257074056157)
pheta2_IRTF = math.asin(r2_IRTF/16.380257074056157)
print(pheta1_IRTF)
print(pheta2_IRTF)

Contents1_IRTF = (x_IRTF[0]*np.sin(pheta1_IRTF))
Contents1_IRTF = Contents1_IRTF/((r1_IRTF*np.cos(pheta1_IRTF)*np.cos(uranus_seangleI_Rads))-(y1_IRTFadj*np.sin(uranus_seangleI_Rads)*np.sin(pheta1_IRTF)))
print(Contents1_IRTF)
Long1_IRTF = 0 - math.atan(Contents1_IRTF)
print((Long1_IRTF*180)/np.pi)

Contents2_IRTF = (np.cos(pheta1_IRTF)*np.sin(uranus_seangleI_Rads))+((y1_IRTFadj*np.sin(pheta1_IRTF)*np.cos(uranus_seangleI_Rads))/r1_IRTF)
print(Contents2_IRTF)
Lat1_IRTF = math.asin(Contents2_IRTF)
print(((Lat1_IRTF*180)/np.pi))

#%%So now lets get a grid of latitudes and longitudes Equator first
x_IRTF = np.flip(np.linspace(-16.5, 16.5, num=34))
b = 0
LongitudeK = []
uranus_seangleI = -27.345
#First lets sort Latitudes
for a in range(68):
    if b < 1 and a < 33:
        r_IRTF = np.sqrt(x_IRTF[a]**2 + y1_IRTF**2)
        if x_IRTF[a] > 16.380257074056157:
            r_IRTF = 16.380257074056157
            pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
            Contents1_IRTF = (16.380257074056157*np.sin(pheta_IRTF))
        else:
            pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
            Contents1_IRTF = (x_IRTF[a]*np.sin(pheta_IRTF))
        Contents1_IRTF = Contents1_IRTF/((r_IRTF*np.cos(pheta_IRTF)*np.cos(uranus_seangleI_Rads))-(y1_IRTF*np.sin(uranus_seangleI_Rads)*np.sin(pheta_IRTF)))
        Long_IRTF = 0 - math.atan(Contents1_IRTF)
        Long_IRTF = (Long_IRTF*180)/np.pi
        LongitudeK.append(Long_IRTF)
    elif b == 0 and a == 33:
        r_IRTF = 16.380257074056157
        pheta_IRTF = math.asin(1)
        Contents1_IRTF = (-16.380257074056157*np.sin(pheta_IRTF))
        Contents1_IRTF = Contents1_IRTF/((r_IRTF*np.cos(pheta_IRTF)*np.cos(uranus_seangleI_Rads))-(y1_IRTF*np.sin(uranus_seangleI_Rads)*np.sin(pheta_IRTF)))
        Long_IRTF = 0 - math.atan(Contents1_IRTF)
        Long_IRTF = (Long_IRTF*180)/np.pi
        LongitudeK.append(Long_IRTF)
        b = 1
    elif b == 1 and a >= 34:
        c = a - 34
        r_IRTF = np.sqrt(x_IRTF[c]**2 + y2_IRTF**2)
        if x_IRTF[c] > 16.380257074056157:
            r_IRTF = 16.380257074056157
            pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
            Contents1_IRTF = (16.380257074056157*np.sin(pheta_IRTF))
        elif x_IRTF[c] < -16.380257074056157:
            r_IRTF = 16.380257074056157
            pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
            Contents1_IRTF = (-16.380257074056157*np.sin(pheta_IRTF))
        else:
            pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
            Contents1_IRTF = (x_IRTF[c]*np.sin(pheta_IRTF))
        Contents1_IRTF = Contents1_IRTF/((r_IRTF*np.cos(pheta_IRTF)*np.cos(uranus_seangleI_Rads))-(y2_IRTF*np.sin(uranus_seangleI_Rads)*np.sin(pheta_IRTF)))
        Long_IRTF = 0 - math.atan(Contents1_IRTF)
        Long_IRTF = (Long_IRTF*180)/np.pi
        if c == 0 or c == 33:
            LongitudeK.append(-1*Long_IRTF)
        else:
            LongitudeK.append(Long_IRTF)

b = 0
LatitudeK = []

for a in range(68):
    if b < 1 and a < 33:
        r_IRTF = np.sqrt(x_IRTF[a]**2 + y1_IRTF**2)
        if r_IRTF > 16.380257074056157:
            r_IRTF = 16.380257074056157
        pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
        Contents2_IRTF = (np.cos(pheta_IRTF)*np.sin(uranus_seangleI_Rads))+((y1_IRTF*np.sin(pheta_IRTF)*np.cos(uranus_seangleI_Rads))/r_IRTF)
        Lat_IRTF = math.asin(Contents2_IRTF)
        Lat_IRTF = (Lat_IRTF*180)/np.pi
        LatitudeK.append(Lat_IRTF)
    elif b == 0 and a == 33:
        r_IRTF = np.sqrt(x_IRTF[a]**2 + y1_IRTF**2)
        if r_IRTF > 16.380257074056157:
            r_IRTF = 16.380257074056157
        pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
        Contents2_IRTF = (np.cos(pheta_IRTF)*np.sin(uranus_seangleI_Rads))+((y1_IRTF*np.sin(pheta_IRTF)*np.cos(uranus_seangleI_Rads))/r_IRTF)
        Lat_IRTF = math.asin(Contents2_IRTF)
        Lat_IRTF = (Lat_IRTF*180)/np.pi
        LatitudeK.append(Lat_IRTF)
        b = 1
    elif b == 1 and a >= 34:
        c = a - 34
        r_IRTF = np.sqrt(x_IRTF[c]**2 + y2_IRTF**2)
        if x_IRTF[c] > 16.380257074056157 or x_IRTF[c] < -16.380257074056157:
            r_IRTF = 16.380257074056157
        pheta_IRTF = math.asin(r_IRTF/16.380257074056157)
        Contents2_IRTF = (np.cos(pheta_IRTF)*np.sin(uranus_seangleI_Rads))+((y2_IRTF*np.sin(pheta_IRTF)*np.cos(uranus_seangleI_Rads))/r_IRTF)
        Lat_IRTF = math.asin(Contents2_IRTF)
        Lat_IRTF = (Lat_IRTF*180)/np.pi
        LatitudeK.append(Lat_IRTF)
        
plt.figure()
plt.plot(LongitudeK, LatitudeK, 'bo')
plt.xlabel('Uranian Longitude (ULS) (Degrees)', fontsize=15)
plt.ylabel('Uranian Latitude (Degrees)', fontsize=15)
plt.title('Expected slit mapping on a 2D surface', fontsize=20)

np.save('Longitudes5.npy', LongitudeK)
np.save('Latitudes5.npy', LatitudeK)

#%% Now lets work out 42 Latitudes and Longitudes along over 54 sets
#Testing script for how to get Uranus at Keck 2006 and IRTF 2016
LONGs = []
LATs = [] #Places_List needs sorting
Places_List = ['10:38', '10:47', '10:55', '11:03', '11:11', '11:15', '11:24', '11:32']

with open('HorizonKeck131014.txt', 'r') as file: #This reads in the file 
    for row in file:
        Data = row.rstrip("\n").split(' ')
        if Data[2] == Places_List[0]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[1]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[2]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[3]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[4]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[5]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[6]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))
        if Data[2] == Places_List[7]:
            LONGs.append(float(Data[6])-171.361857125)
            LATs.append(float(Data[8]))

Places_List = ['10:30']

with open('C:/Users/snowy/OneDrive/Documents/Python work/Keck 12OCT14/HorizonKeck121014.txt', 'r') as file: #This reads in the file 
    for row in file:
        Data = row.rstrip("\n").split(' ')
        if Data[2] == Places_List[0]:
            LONGs.append(float(Data[21])+223)
            LATs.append(float(Data[23]))

#%% Now to calculate 108 lines (Latitudes are all the same so same 2 lines repeated)
#Need to work out how to polygon each pixel! !!!!!!!!
from PIL import Image, ImageDraw
width = 360
height = 180

img = Image.new('L', (width, height), 0)

for i in range(33):
    polygon = [(LongitudeK[i]-np.nanmin(LongitudeK), LatitudeK[i]+90), (LongitudeK[1+i]-np.nanmin(LongitudeK), LatitudeK[1+i]+90), (LongitudeK[35+i]-np.nanmin(LongitudeK), LatitudeK[35+i]+90), (LongitudeK[34+i]-np.nanmin(LongitudeK), LatitudeK[34+i]+90)]
    ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
        
#ImageDraw.Draw(img).polygon(polygon, outline=1, fill=1)
#ImageDraw.Draw(img).polygon(polygon2, outline=1, fill=1)
Slit0 = np.flipud(np.array(img))

plt.figure()
plt.imshow(Slit0)

#%% So now to make 54 images of each slit 
SLITS = []
SLITShapes = [] # Use Mapping_Q1A for the pixel values remember the flicking about from North Equator and South

Mapping_Q1A = np.load('Ints_1hr_Keck.npy')
Mapping_Q1B = np.load('Ints_1hr_Keck12.npy')
Mean_Err_Q1B = np.nanmean(np.load('Ints_1hr_Keck12_Err.npy'))
Mean_Mapping_Q1A = np.nanmean(Mapping_Q1A, axis = 0)
Mean_Mapping_Q1B = np.nanmean(Mapping_Q1B, axis = 0)
Slits = []
SlitsB = []
SlitImg = [] #Assuming from Laurent's CML that 11:05 == 40 CML
SlitNo = []

for A in range(9):
    width = 360
    height = 180
    Pixels = []
    Pixels2 = []
    Slit1 = np.zeros((180, 360))
    Slit2 = np.zeros((180, 360))
    SLITShapesB= []
    Slit = np.zeros((180, 360))
    SLITSHAPE = np.zeros((180, 360))
    if A < 8:
        for i in range(30):
            ii = i
            img2 = Image.new('L', (width, height), 0)
            polygon = [(LongitudeK[ii+2]+float(LONGs[A])+180, LatitudeK[ii+2]+90), (LongitudeK[3+ii]+float(LONGs[A])+180, LatitudeK[3+ii]+90), (LongitudeK[37+ii]+float(LONGs[A])+180, LatitudeK[37+ii]+90), (LongitudeK[36+ii]+float(LONGs[A])+180, LatitudeK[36+ii]+90)]
            Fill_Data = 10**6*(Mapping_Q1A[A,i] - Mean_Mapping_Q1A[i])
            if i % 2 == 0:
                Pixel = ImageDraw.Draw(img2).polygon(polygon, outline=1, fill=1)
                PixelA = np.flipud(np.array(img2)*Fill_Data)
                PixelB = np.flipud(np.array(img2))
                Pixels.append(PixelA)
                Pixels2.append(PixelB)
            else:
                Pixel = ImageDraw.Draw(img2).polygon(polygon, outline=1, fill=1)
                PixelA = np.flipud(np.array(img2)*Fill_Data)
                PixelB = np.flipud(np.array(img2))
                Pixels.append(PixelA)
                Pixels2.append(PixelB)
        for i in range(30):
            Slit1 = np.add(Slit1, Pixels[i])
            Slit2 = np.add(Slit2, Pixels2[i])
        Slits.append(Slit1)
        SlitsB.append(Slit2)
        plt.figure()
        plt.imshow(Slit1/Slit2)
    else:
        for i in range(30):
            ii = i
            img2 = Image.new('L', (width, height), 0)
            polygon = [(LongitudeK[ii+2]+float(LONGs[A])-180, LatitudeK[ii+2]+90), (LongitudeK[3+ii]+float(LONGs[A])-180, LatitudeK[3+ii]+90), (LongitudeK[37+ii]+float(LONGs[A])-180, LatitudeK[37+ii]+90), (LongitudeK[36+ii]+float(LONGs[A])-180, LatitudeK[36+ii]+90)]
            Fill_Data = (Mapping_Q1B[30-i] - Mean_Mapping_Q1A[i]*10**6)
            if i % 2 == 0:
                Pixel = ImageDraw.Draw(img2).polygon(polygon, outline=1, fill=1)
                PixelA = np.flipud(np.array(img2)*Fill_Data)
                PixelB = np.flipud(np.array(img2))
                Pixels.append(PixelA)
                Pixels2.append(PixelB)
            else:
                Pixel = ImageDraw.Draw(img2).polygon(polygon, outline=1, fill=1)
                PixelA = np.flipud(np.array(img2)*Fill_Data)
                PixelB = np.flipud(np.array(img2))
                Pixels.append(PixelA)
                Pixels2.append(PixelB)
        for i in range(30):
            Slit1 = np.add(Slit1, Pixels[i])
            Slit2 = np.add(Slit2, Pixels2[i])

        Slits.append(Slit1)
        SlitsB.append(Slit2)
        plt.figure()
        plt.imshow(Slit1/Slit2)
        
SlitImgT = (Slits[0] + Slits[1] + Slits[2] + Slits[3] + Slits[4] + Slits[5] + Slits[6] + Slits[7])/8
SlitImg.append(SlitImgT)
SlitT = (SlitsB[0] + SlitsB[1] + SlitsB[2] + SlitsB[3] + SlitsB[4] + SlitsB[5] + SlitsB[6] + SlitsB[7])/8
SlitNo.append(SlitT)

plt.figure()
plt.imshow(SlitImgT/SlitT)

plt.figure()
plt.imshow(Slits[8]/SlitsB[8])

#%% Now we need to make two arrays (One that lets us know how many pixels intersect and the other to just keep adding the values together)
import scipy.ndimage as f

#Now to add all FIN_SLITS together and divide it by Slits Mask
#Lets seperate these out and then stack them together
IRTF_16_MAP = SlitImgT/SlitT
KECK_16_MAP_0 = IRTF_16_MAP[:,180:]
KECK_16_MAP_1 = IRTF_16_MAP[:,0:180]

# KECK_16_MAP = Slits[8]/SlitsB[8]
# KECK_16_MAP2_0 = KECK_16_MAP[:,180:]
# KECK_16_MAP2_1 = KECK_16_MAP[:,0:180]

KECK_16_MAP = np.hstack((KECK_16_MAP_0, KECK_16_MAP_1))
KECK_16_MAP = np.nan_to_num(KECK_16_MAP, nan=0)

# KECK_16_MAP2v1 = np.hstack((KECK_16_MAP2_1, KECK_16_MAP2_0))
# KECK_16_MAP2v1 = np.nan_to_num(KECK_16_MAP2v1, nan=0)

# FIN_KECK_16_MAP = (KECK_16_MAP + KECK_16_MAP2v1)/2
# FIN_KECK_16_MAP[FIN_KECK_16_MAP == 0] = np.nan
FIN_KECK_16_MAP = KECK_16_MAP
FIN_KECK_16_MAP[FIN_KECK_16_MAP == 0] = np.nan

plt.figure()
plt.imshow(FIN_KECK_16_MAP, vmax = 0.5360874077089963, vmin = -0.09651973255209997, cmap='gist_heat')
cbar = plt.colorbar()
cbar.set_label(r'Intensity Difference from Average Spectrum (${\mu}Wm^{-2}sr^{-1}$) ' + u'\u00B1' + ' 0.02 ${\mu}Wm^{-2}sr^{-1}$', fontsize=20) 
cbar.ax.tick_params(labelsize=20)
plt.title(r'Intensity Difference from Average $H_{3}^{+}$ Emissions at Uranus of $Q(1,0^{-})$ and $Q(3,0^{-})$ Emission Lines', pad = 100, fontsize = 30)
plt.xlabel('Longitude across Uranus (ULS) ($^\circ$)', fontsize=25, labelpad=11)
plt.ylabel('Latitude across Uranus (ULS) ($^\circ$)', fontsize=25, labelpad=11)
plt.xticks(np.arange(0, 361, 20), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$', '200$^\circ$', '220$^\circ$', '240$^\circ$', '260$^\circ$', '280$^\circ$', '300$^\circ$', '320$^\circ$', '340$^\circ$', '360$^\circ$'], fontsize=20)
#plt.xticks(np.arange(0, 181, 20), ['0$^\circ$', '20$^\circ$', '40$^\circ$', '60$^\circ$', '80$^\circ$', '100$^\circ$', '120$^\circ$', '140$^\circ$', '160$^\circ$', '180$^\circ$'], fontsize=11)
plt.yticks(np.arange(0, 181, 30), ['90$^\circ$N', '60$^\circ$N', '30$^\circ$N', '0$^\circ$', '30$^\circ$S', '60$^\circ$S', '90$^\circ$S'], fontsize=20)
#plt.vlines((0, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000, 3200, 3400, 3600), 0, 3600, colors = 'k', linestyles = 'dotted', alpha = 0.4)
#plt.vlines((0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280, 300, 320, 340, 360), 0, 360, colors = 'k', linestyles = 'dotted', alpha = 0.75)
#plt.hlines((0, 30, 60, 90, 120, 150, 180), 0, 360, colors = 'k', linestyles = 'dotted', alpha = 0.75)
plt.ylim(180, 0) #1801, -1
plt.xlim(0, 360) #0, 3601
