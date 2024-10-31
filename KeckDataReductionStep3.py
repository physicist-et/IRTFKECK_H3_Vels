# -*- coding: utf-8 -*-
"""
Created on Mon Jan 31 15:46:08 2022

@author: snowy
"""
import numpy as np
from KeckDataReductionStep1 import Spat_Ratio, A2Star1P, A2Star1N
from matplotlib.patches import Rectangle, Ellipse
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from PIL import Image, ImageDraw

#Need to set up Uranus a lot better (Expected brightness across Uranus plus shell from atmosphere (look to Henrik's thesis for starting points))
#%%
uranus_seangleI = (-35.773649 + -35.784010)/2
eqwidth=3.71875
equatorkm=25559*2
#stretch yy to become a sphere
flattening =0.0229
pix_to_arc = np.nanmean(Spat_Ratio)
uranus_seangleI_Rads = (uranus_seangleI*np.pi)/180
losflattening=flattening*(1-np.sin(uranus_seangleI_Rads))
eq_po_ratioI=1-losflattening
po_eq_ratioI=1/eq_po_ratioI
limbdist=0.5*eqwidth/pix_to_arc
A2 = np.nanmean(A2Star1P+A2Star1N)/2

#Now we make the background for the Uranus ball
Back_Img = np.zeros((100,100))
# fig,ax = plt.subplots(1)
# ax.set_aspect('equal')

# Show the image
# ax.imshow(Back_Img)

# Now, loop through coord arrays, and create a circle at each x,y pair
circ = Ellipse((100,100), po_eq_ratioI*limbdist, limbdist)
Convol_Img = gaussian_filter(circ, sigma=A2)
img = Image.new('L', (200, 200), 1)
Tot_Img = ImageDraw.Draw(img).ellipse([(86.6518096277065, 87.13619766424796), (113.3481903722935, 112.86380233575204)], outline=1, fill=0)
Convol_Img = gaussian_filter(img, sigma=0.64)

Slit_L = 24/pix_to_arc
Slit_W = 0.288/pix_to_arc
Slit = Rectangle((100-(Slit_W/2), 100-(Slit_L/2)), Slit_W, Slit_L, linewidth=0, facecolor='blue', alpha=0.5)

# Show the image
fig, ax = plt.subplots()
ax.imshow(img, cmap='gist_gray')
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')
ax.set_title('Approximate view of Uranus')
fig, ax = plt.subplots()
ax.imshow(Convol_Img, cmap='gist_gray')
ax.add_patch(Slit)
plt.show()
ax.set_xlabel('Pixels')
ax.set_ylabel('Pixels')
ax.set_title('Approximate view of Uranus and observation slit with a Gaussian blur from sigma of Star observations')

#Now do a gaussian on the Total_Img


#%% Now we want to put on the slit as well as blur Uranus to get a rough idea of how Urnaus will look
# Slit_L = 24/pix_to_arc
# Slit_W = 0.288/pix_to_arc

# polygon = [(100-Slit_W, 100-Slit_L), (100+Slit_W, 100+Slit_L)]
# ImageDraw.Draw(Convol_Img).polygon(polygon, outline=0, fill=2)

# plt.figure()
# plt.imshow(Convol_Img, cmap='gist_gray')
# plt.xlabel('Pixels')
# plt.ylabel('Pixels')
# plt.title('Approximate view of Uranus with a gaussian blur from sigma of Star observations')
