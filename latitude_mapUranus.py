# -*- coding: utf-8 -*-
"""
Created on Fri May 13 13:54:51 2022

@author: mnc8
"""


def Keck_spectra_viewer(array_to_view, title):
    """
    Input/(s):
        array_to_view --- the array of the spectra that we want to visualise
        title --- a title to use for the output plot, or whatever is outputted
    """
    #Modules to use
    import matplotlib.pyplot as plt          #For plots
    # from matplotlib_customs import grayscale_cmap        #For the correct colormap of the spectra
    import numpy as np
    from skimage import exposure            #For rescaling

    #Contrast stretching a la tvscl
    p_lower, p_higher = np.percentile(array_to_view, (0.5, 99.5))
    array_to_view = exposure.rescale_intensity(array_to_view, in_range=(p_lower, p_higher))

    #Visualise the array
    plt.figure(figsize=(7,5), dpi=200)
    # plt.imshow(array_to_view, cmap = grayscale_cmap('cubehelix'), origin = 'lower')
    plt.imshow(array_to_view, cmap = 'brg', origin = 'lower')
    plt.title(title, fontsize=14)
#    plt.colorbar(extend='both')
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    # plt.savefig('C:\\Users\\mnc8\\Google Drive\\PhD\\Thesis_Final_Draft\\Figures\\Methodology\\star_spectrum.pdf')
    plt.show()
    
    return


# def modelled_image(Day_Number):
Day_Number = 1
"""
Input/(s):
    Day_Number --- this input (1, 2, 3, 4, 5, 6, or 7) will determine the parameters that will be fed into the rest of the script, such as seeing etc.
Output/(s):
    fake1 --- the modelled image of reflected sunlight from the planet with line-of-sight correction factored in
    los1 --- the modelled line-of-sight correction factor for the observing night in question
    lat1 --- the latitude map of the modelled planet image for the night in question
    lon1 --- the longitude map of the modelled planet image for the night in question
"""
#Modules for import
import numpy as np
#from Keck_data_mapper import dist_ellipse

import scipy.ndimage as congrid
from scipy import signal

#Float the topmost variables depending on the inputted Day_Number value
if Day_Number == 1 or Day_Number == 2 or Day_Number == 3 or Day_Number == 4:
    slit_width = 0.288
elif Day_Number == 5 or Day_Number == 6 or Day_Number == 7:
    slit_width = 0.288

#Float all the other topmost variables that are common
#    planet_angular_diameter, planet_sub_Earth_angle, planet_CML, seeing, planet_equatorial_diameter, modelled_image_plate_scale, planet_flattening = 17.53142857142857, 31.844285714285714, 180., 0.5985714285714285, 116460., 0.025, 0.09796
planet_angular_diameter, planet_sub_Earth_angle, planet_CML, seeing, planet_equatorial_diameter, modelled_image_plate_scale, planet_flattening = 3.7056, -27.340, 180, 1.65, 25559., 0.1445, 0.02293
#New values following addition of Day Eight and Day Nine spectra. Ang. Dia. is average of Day Eight values, the widest limb-to-limb day of the dataset. SEL is average of all Night SELs.
# planet_angular_diameter, planet_sub_Earth_angle = 18.32718, 31.81

"""
c) dist_ellipse --- creates a grid of values needed for Tom's version of mapping scripts
"""
#DISTELLIPSE
def dist_ellipse(n, xc, yc, ratio, pa=0): # original implementation (like DIST_ELLIPSE IDL function). Raked from https://github.com/SKIRT/PTS/blob/master/magic/dist_ellipse.py  on  5th May 2020.

    """
    Input/(s):
    N = either  a scalar specifying the size of the N x N square output
              array, or a 2 element vector specifying the size of the
               M x N rectangular output array.
       XC,YC - Scalars giving the position of the ellipse center.   This does
               not necessarily have to be within the image
       RATIO - Scalar giving the ratio of the major to minor axis.   This
               should be greater than 1 for position angle to have its
               standard meaning.
    OPTIONAL INPUTS:
      POS_ANG - Position angle of the major axis in degrees, measured counter-clockwise
               from the Y axis.  For an image in standard orientation
               (North up, East left) this is the astronomical position angle.
               Default is 0 degrees.
    OUTPUT:
       IM - REAL*4 elliptical mask array, of size M x N.  THe value of each
               pixel is equal to the semi-major axis of the ellipse of center
                XC,YC, axial ratio RATIO, and position angle POS_ANG, which
               passes through the pixel.
    """
    #Modules for import
    import numpy as np
    
    ang = np.radians(pa + 90.)
    cosang = np.cos(ang)
    sinang = np.sin(ang)
    nx = n[1]
    ny = n[0]
    x = np.arange(-xc,nx-xc)
    y = np.arange(-yc,ny-yc)

    im = np.empty(n)
    xcosang = x*cosang
    xsinang = x*sinang

    for i in range(0, ny):

        xtemp = xcosang + y[i]*sinang
        ytemp = -xsinang + y[i]*cosang
        im[i,:] = np.sqrt((xtemp*ratio)**2 + ytemp**2)

    return im


#Start the copying and pasting in othe "stripped" version of the script proper
"""
mapcoords_width.pro
Very first thing done is a summation of the latitude and longitude array and converting them to degrees.
So, I need the latitude and longitude arrays from whichever script they're in...
They come from makeimage2
"""
#makeimage2
#Float topmost variables
sat_eqwidth, sat_seangle, sat_posangle, equatorkm = planet_angular_diameter, planet_sub_Earth_angle, 90., planet_equatorial_diameter
R1 = equatorkm/2.
R2, R3 = (R1+1000.)/R1, (R1+3000.)/R1
sat_cml, plate_scale, flattening = planet_CML, modelled_image_plate_scale, planet_flattening
losflattening = flattening * (1. - np.sin(np.deg2rad(sat_seangle)))
eq_po_ratio = 1. - losflattening
po_eq_ratio = 1./eq_po_ratio
limbdist = (0.5 * sat_eqwidth * R3)/plate_scale
#"MAPPING"
cox, coy = 100, 100
planetmask = np.zeros((200,200))
offsetxx, offsetyy = np.arange(200) - cox, np.arange(200) - coy
offsetx = np.zeros((2,200))
offsetx[0,:], offsetx[1,:] = offsetxx, offsetxx
offsetx = congrid.zoom(offsetx, [100,1], order=0, mode='nearest')
offsety = np.zeros((200,2))
offsety[:,0], offsety[:,1] = offsetyy, offsetyy
offsety = congrid.zoom(offsety, [1,100], order=0, mode='nearest')
im2 = dist_ellipse((200,200), cox, coy, po_eq_ratio, pa=(sat_posangle + 90.))
planetpos = np.where((im2 < (limbdist + 1.)))
planetmask[planetpos] = 1.
calctheseplanetpixels = np.where((planetmask > 0.))
ctpp_xxx, ctpp_yyy = np.zeros((200,200)), np.zeros((200,200))
ctpp_xxx[calctheseplanetpixels] = offsetx[calctheseplanetpixels]
ctpp_yyy[calctheseplanetpixels] = offsety[calctheseplanetpixels]
sat_pixel_radius = (0.5 * R3 * sat_eqwidth)/0.1445
latitude_array = np.zeros((4,200,200))
longitude_array = np.zeros((4,200,200))
#We now need the calculation of latitude and longitude, so nirspec_calc_lat_long_image
#%%
"""
nirspec_calc_lat_long_image
"""
#Float the variables at the top
x, y, jup_pixel_radius, jup_seangle, jup_cml, flattening = ctpp_xxx, ctpp_yyy, sat_pixel_radius, sat_seangle, sat_cml, flattening
R = jup_pixel_radius
#Float empty arrays for the four corner steps next. One day you have to devise a for loop for these...
x2_0 = x + 0.5
y2_0 = y + 0.5
x2_1 = x + 0.5
y2_1 = y - 0.5
x2_2 = x - 0.5
y2_2 = y - 0.5
x2_3 = x - 0.5
y2_3 = y + 0.5
#No coordinate transformations necessary since posangle is always 0 for now
#Straight onto the next step then...
xx_0 = x2_0
xx_1 = x2_1
xx_2 = x2_2
xx_3 = x2_3
#The yy's have to be stretched to become a sphere
losflattening = flattening * (1. - np.sin(np.deg2rad(jup_seangle)))
eq_po_ratio = 1. - losflattening
yy_0 = y2_0/eq_po_ratio
yy_1 = y2_1/eq_po_ratio
yy_2 = y2_2/eq_po_ratio
yy_3 = y2_3/eq_po_ratio
#Proper distance from centre
pp_0 = np.sqrt(xx_0**2 + yy_0**2)
pp_1 = np.sqrt(xx_1**2 + yy_1**2)
pp_2 = np.sqrt(xx_2**2 + yy_2**2)
pp_3 = np.sqrt(xx_3**2 + yy_3**2)
#Limit calculations next
goodvalues_0 = np.where(((pp_0/R) <= 0.998))
goodvalues_1 = np.where(((pp_1/R) <= 0.998))
goodvalues_2 = np.where(((pp_2/R) <= 0.998))
goodvalues_3 = np.where(((pp_3/R) <= 0.998))
if len(goodvalues_0[0]) > 0:
    #Float empty arrays
    lats_0, lons_0, cc_0, lons_x_0, lons_y_0 = np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0)))
    #Angular distance from centre
    cc_0[goodvalues_0] = np.arcsin(pp_0[goodvalues_0]/R)
    lons_x_0[goodvalues_0] = (xx_0[goodvalues_0] * np.sin(cc_0[goodvalues_0]))
    lons_y_0[goodvalues_0] = (pp_0[goodvalues_0] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_0[goodvalues_0])) - (yy_0[goodvalues_0] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_0[goodvalues_0]))
    lons_offset_0 = np.mod(np.arctan2(lons_y_0, lons_x_0), 2.*np.pi)
    lons_0 = np.mod(((np.deg2rad(jup_cml) + lons_offset_0) - (np.pi/2.)), 2.*np.pi)
    lats_0[goodvalues_0] = np.arcsin((np.cos(cc_0[goodvalues_0]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_0[goodvalues_0] * np.sin(cc_0[goodvalues_0]) * np.cos(np.deg2rad(jup_seangle)))/pp_0[goodvalues_0]))
if len(goodvalues_1[0]) > 0:
    #Float empty arrays
    lats_1, lons_1, cc_1, lons_x_1, lons_y_1 = np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1)))
    #Angular distance from centre
    cc_1[goodvalues_1] = np.arcsin(pp_1[goodvalues_1]/R)
    lons_x_1[goodvalues_1] = (xx_1[goodvalues_1] * np.sin(cc_1[goodvalues_1]))
    lons_y_1[goodvalues_1] = (pp_1[goodvalues_1] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_1[goodvalues_1])) - (yy_1[goodvalues_1] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_1[goodvalues_1]))
    lons_offset_1 = np.mod(np.arctan2(lons_y_1, lons_x_1), 2.*np.pi)
    lons_1 = np.mod(((np.deg2rad(jup_cml) + lons_offset_1) - (np.pi/2.)), 2.*np.pi)
    lats_1[goodvalues_1] = np.arcsin((np.cos(cc_1[goodvalues_1]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_1[goodvalues_1] * np.sin(cc_1[goodvalues_1]) * np.cos(np.deg2rad(jup_seangle)))/pp_1[goodvalues_1]))
if len(goodvalues_2[0]) > 0:
    #Float empty arrays
    lats_2, lons_2, cc_2, lons_x_2, lons_y_2 = np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2)))
    #Angular distance from centre
    cc_2[goodvalues_2] = np.arcsin(pp_2[goodvalues_2]/R)
    lons_x_2[goodvalues_2] = (xx_2[goodvalues_2] * np.sin(cc_2[goodvalues_2]))
    lons_y_2[goodvalues_2] = (pp_2[goodvalues_2] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_2[goodvalues_2])) - (yy_2[goodvalues_2] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_2[goodvalues_2]))
    lons_offset_2 = np.mod(np.arctan2(lons_y_2, lons_x_2), 2.*np.pi)
    lons_2 = np.mod(((np.deg2rad(jup_cml) + lons_offset_2) - (np.pi/2.)), 2.*np.pi)
    lats_2[goodvalues_2] = np.arcsin((np.cos(cc_2[goodvalues_2]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_2[goodvalues_2] * np.sin(cc_2[goodvalues_2]) * np.cos(np.deg2rad(jup_seangle)))/pp_2[goodvalues_2]))
if len(goodvalues_3[0]) > 0:
    #Float empty arrays
    lats_3, lons_3, cc_3, lons_x_3, lons_y_3 = np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3)))
    #Angular distance from centre
    cc_3[goodvalues_3] = np.arcsin(pp_3[goodvalues_3]/R)
    lons_x_3[goodvalues_3] = (xx_3[goodvalues_3] * np.sin(cc_3[goodvalues_3]))
    lons_y_3[goodvalues_3] = (pp_3[goodvalues_3] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_3[goodvalues_3])) - (yy_3[goodvalues_3] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_3[goodvalues_3]))
    lons_offset_3 = np.mod(np.arctan2(lons_y_3, lons_x_3), 2.*np.pi)
    lons_3 = np.mod(((np.deg2rad(jup_cml) + lons_offset_3) - (np.pi/2.)), 2.*np.pi)
    lats_3[goodvalues_3] = np.arcsin((np.cos(cc_3[goodvalues_3]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_3[goodvalues_3] * np.sin(cc_3[goodvalues_3]) * np.cos(np.deg2rad(jup_seangle)))/pp_3[goodvalues_3]))
#Now store each of these into empty arrays called latit and longit
latit = np.zeros((4,200,200))
longit = np.zeros((4,200,200))
latit[0,:,:], latit[1,:,:], latit[2,:,:], latit[3,:,:] = lats_0, lats_1, lats_2, lats_3
longit[0,:,:], longit[1,:,:], longit[2,:,:], longit[3,:,:] = lons_0, lons_1, lons_2, lons_3
#That's the final output from calc_lat_long_image so let's return to makeimage2

#%%
"""
makeimage2 (CONTINUED)
"""
#looks like we continue by placing the latit and longit outputs from above into the longitude and latitude arrays from earlier. This is done for each corner
cornerarray = np.zeros((200,200))
#Corner 0
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[0,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[0,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[0,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[0,:,:] = cornerarray
#Corner 1
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[1,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[1,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[1,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[1,:,:] = cornerarray
#Corner 2
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[2,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[2,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[2,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[2,:,:] = cornerarray
#Corner 3
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[3,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[3,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[3,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[3,:,:] = cornerarray
#The IDL reform statement is not needed here since there are no degenerate leading axes when extracting one of four corner data
#Tom then extracts the first corner from latitude_array and calls it fakeimage
fakeimage = latitude_array[0,:,:]
fakeimage = np.abs(fakeimage)
offlimb2 = np.where((fakeimage <= 0.))
fakeimage = np.cos((fakeimage*2))
fakeimage = fakeimage - np.amin(fakeimage)
fakeimage = fakeimage/np.amax(fakeimage) * 0.9
fakeimage = fakeimage + 0.1
fakeimage[offlimb2] = 0.
fakeimage = np.rot90(fakeimage, k=2)
#All the corners are then flipped in the latitde and longitude arrays
latitude_array[0,:,:] = np.rot90(latitude_array[0,:,:], k=2)
latitude_array[1,:,:] = np.rot90(latitude_array[1,:,:], k=2)
latitude_array[2,:,:] = np.rot90(latitude_array[2,:,:], k=2)
latitude_array[3,:,:] = np.rot90(latitude_array[3,:,:], k=2)
corner = np.arange(4)
for i in corner:
    longitude_array[i,:,:] = np.rot90(longitude_array[i,:,:], k=2)
#And that's it for makeimage2

#%%
"""
return to mapcoords_width.pro
"""
#Looks like Tom averages(?) the latitude and longitude corner arrays using IDL's total() and then converts them to degrees.
latarray = np.rad2deg(sum(latitude_array)/len(latitude_array))
lonarray = np.rad2deg(sum(longitude_array)/len(longitude_array))
#Need to go into makeimage3 now for the fakek thing which is the los_coor

"""
makeimage3
"""
#Float variables at the top
R2, R3, plate_scale = (R1+1000.)/R1, (R1+3000.)/R1, plate_scale
limbdist = (0.5 * sat_eqwidth)/plate_scale
losflattening = flattening * (1. - np.sin(np.deg2rad(sat_seangle)))
eq_po_ratio = 1. - losflattening
po_eq_ratio = 1./eq_po_ratio
naxis1, cox, coy = 200, 100, 100
#"MAPPING"
planetmask = np.zeros((200,200))
offsetxx, offsetyy = np.arange(200) - cox, np.arange(200) - coy
offsetx = np.zeros((2,200))
offsetx[0,:], offsetx[1,:] = offsetxx, offsetxx
offsetx = congrid.zoom(offsetx, [100,1], order=0, mode='nearest')
offsety = np.zeros((200,2))
offsety[:,0], offsety[:,1] = offsetyy, offsetyy
offsety = congrid.zoom(offsety, [1,100], order=0, mode='nearest')
im2 = dist_ellipse((200,200), cox, coy, po_eq_ratio, pa=(sat_posangle + 90.))
planetpos = np.where((im2 < (limbdist + 1.)))
planetmask[planetpos] = 1.
calctheseplanetpixels = np.where((planetmask > 0.))
ctpp_xxx, ctpp_yyy = np.zeros((200,200)), np.zeros((200,200))
ctpp_xxx[calctheseplanetpixels] = offsetx[calctheseplanetpixels]
ctpp_yyy[calctheseplanetpixels] = offsety[calctheseplanetpixels]
sat_pixel_radius = (0.5 * R3 * sat_eqwidth)/0.1445
latitude_array = np.zeros((4,200,200))
longitude_array = np.zeros((4,200,200))
#We now need the calculation of latitude and longitude, so nirspec_calc_lat_long_image

#%%
"""
nirspec_calc_lat_long_image
"""
#Float the variables at the top
x, y, jup_pixel_radius, jup_seangle, jup_cml, flattening = ctpp_xxx, ctpp_yyy, sat_pixel_radius, sat_seangle, sat_cml, flattening
R = jup_pixel_radius
#Float empty arrays for the four corner steps next. One day you have to devise a for loop for these...
x2_0 = x + 0.5
y2_0 = y + 0.5
x2_1 = x + 0.5
y2_1 = y - 0.5
x2_2 = x - 0.5
y2_2 = y - 0.5
x2_3 = x - 0.5
y2_3 = y + 0.5
#No coordinate transformations necessary since posangle is always 0 for now
#Straight onto the next step then...
xx_0 = x2_0
xx_1 = x2_1
xx_2 = x2_2
xx_3 = x2_3
#The yy's have to be stretched to become a sphere
losflattening = flattening * (1. - np.sin(np.deg2rad(jup_seangle)))
eq_po_ratio = 1. - losflattening
yy_0 = y2_0/eq_po_ratio
yy_1 = y2_1/eq_po_ratio
yy_2 = y2_2/eq_po_ratio
yy_3 = y2_3/eq_po_ratio
#Proper distance from centre
pp_0 = np.sqrt(xx_0**2 + yy_0**2)
pp_1 = np.sqrt(xx_1**2 + yy_1**2)
pp_2 = np.sqrt(xx_2**2 + yy_2**2)
pp_3 = np.sqrt(xx_3**2 + yy_3**2)
#Limit calculations next
goodvalues_0 = np.where(((pp_0/R) <= 0.998))
goodvalues_1 = np.where(((pp_1/R) <= 0.998))
goodvalues_2 = np.where(((pp_2/R) <= 0.998))
goodvalues_3 = np.where(((pp_3/R) <= 0.998))
if len(goodvalues_0[0]) > 0:
    #Float empty arrays
    lats_0, lons_0, cc_0, lons_x_0, lons_y_0 = np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0))), np.zeros((len(pp_0), len(pp_0)))
    #Angular distance from centre
    cc_0[goodvalues_0] = np.arcsin(pp_0[goodvalues_0]/R)
    lons_x_0[goodvalues_0] = (xx_0[goodvalues_0] * np.sin(cc_0[goodvalues_0]))
    lons_y_0[goodvalues_0] = (pp_0[goodvalues_0] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_0[goodvalues_0])) - (yy_0[goodvalues_0] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_0[goodvalues_0]))
    lons_offset_0 = np.mod(np.arctan2(lons_y_0, lons_x_0), 2.*np.pi)
    lons_0 = np.mod(((np.deg2rad(jup_cml) + lons_offset_0) - (np.pi/2.)), 2.*np.pi)
    lats_0[goodvalues_0] = np.arcsin((np.cos(cc_0[goodvalues_0]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_0[goodvalues_0] * np.sin(cc_0[goodvalues_0]) * np.cos(np.deg2rad(jup_seangle)))/pp_0[goodvalues_0]))
if len(goodvalues_1[0]) > 0:
    #Float empty arrays
    lats_1, lons_1, cc_1, lons_x_1, lons_y_1 = np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1))), np.zeros((len(pp_1), len(pp_1)))
    #Angular distance from centre
    cc_1[goodvalues_1] = np.arcsin(pp_1[goodvalues_1]/R)
    lons_x_1[goodvalues_1] = (xx_1[goodvalues_1] * np.sin(cc_1[goodvalues_1]))
    lons_y_1[goodvalues_1] = (pp_1[goodvalues_1] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_1[goodvalues_1])) - (yy_1[goodvalues_1] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_1[goodvalues_1]))
    lons_offset_1 = np.mod(np.arctan2(lons_y_1, lons_x_1), 2.*np.pi)
    lons_1 = np.mod(((np.deg2rad(jup_cml) + lons_offset_1) - (np.pi/2.)), 2.*np.pi)
    lats_1[goodvalues_1] = np.arcsin((np.cos(cc_1[goodvalues_1]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_1[goodvalues_1] * np.sin(cc_1[goodvalues_1]) * np.cos(np.deg2rad(jup_seangle)))/pp_1[goodvalues_1]))
if len(goodvalues_2[0]) > 0:
    #Float empty arrays
    lats_2, lons_2, cc_2, lons_x_2, lons_y_2 = np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2))), np.zeros((len(pp_2), len(pp_2)))
    #Angular distance from centre
    cc_2[goodvalues_2] = np.arcsin(pp_2[goodvalues_2]/R)
    lons_x_2[goodvalues_2] = (xx_2[goodvalues_2] * np.sin(cc_2[goodvalues_2]))
    lons_y_2[goodvalues_2] = (pp_2[goodvalues_2] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_2[goodvalues_2])) - (yy_2[goodvalues_2] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_2[goodvalues_2]))
    lons_offset_2 = np.mod(np.arctan2(lons_y_2, lons_x_2), 2.*np.pi)
    lons_2 = np.mod(((np.deg2rad(jup_cml) + lons_offset_2) - (np.pi/2.)), 2.*np.pi)
    lats_2[goodvalues_2] = np.arcsin((np.cos(cc_2[goodvalues_2]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_2[goodvalues_2] * np.sin(cc_2[goodvalues_2]) * np.cos(np.deg2rad(jup_seangle)))/pp_2[goodvalues_2]))
if len(goodvalues_3[0]) > 0:
    #Float empty arrays
    lats_3, lons_3, cc_3, lons_x_3, lons_y_3 = np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3))), np.zeros((len(pp_3), len(pp_3)))
    #Angular distance from centre
    cc_3[goodvalues_3] = np.arcsin(pp_3[goodvalues_3]/R)
    lons_x_3[goodvalues_3] = (xx_3[goodvalues_3] * np.sin(cc_3[goodvalues_3]))
    lons_y_3[goodvalues_3] = (pp_3[goodvalues_3] * np.cos(np.deg2rad(jup_seangle)) * np.cos(cc_3[goodvalues_3])) - (yy_3[goodvalues_3] * np.sin(np.deg2rad(jup_seangle)) * np.sin(cc_3[goodvalues_3]))
    lons_offset_3 = np.mod(np.arctan2(lons_y_3, lons_x_3), 2.*np.pi)
    lons_3 = np.mod(((np.deg2rad(jup_cml) + lons_offset_3) - (np.pi/2.)), 2.*np.pi)
    lats_3[goodvalues_3] = np.arcsin((np.cos(cc_3[goodvalues_3]) * np.sin(np.deg2rad(jup_seangle))) + ((yy_3[goodvalues_3] * np.sin(cc_3[goodvalues_3]) * np.cos(np.deg2rad(jup_seangle)))/pp_3[goodvalues_3]))
#Now store each of these into empty arrays called latit and longit
latit = np.zeros((4,200,200))
longit = np.zeros((4,200,200))
latit[0,:,:], latit[1,:,:], latit[2,:,:], latit[3,:,:] = lats_0, lats_1, lats_2, lats_3
longit[0,:,:], longit[1,:,:], longit[2,:,:], longit[3,:,:] = lons_0, lons_1, lons_2, lons_3
#That's the final output from calc_lat_long_image so let's return to makeimage2

#%%
"""
makeimage2
"""
#looks like we continue by placing the latit and longit outputs from above into the longitude and latitude arrays from earlier. This is done for each corner
cornerarray = np.zeros((200,200))
#Corner 0
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[0,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[0,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[0,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[0,:,:] = cornerarray
#Corner 1
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[1,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[1,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[1,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[1,:,:] = cornerarray
#Corner 2
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[2,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[2,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[2,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[2,:,:] = cornerarray
#Corner 3
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = latit[3,calctheseplanetpixels[0],calctheseplanetpixels[1]]
latitude_array[3,:,:] = cornerarray
cornerarray[:] = 0.
cornerarray[calctheseplanetpixels] = longit[3,calctheseplanetpixels[0],calctheseplanetpixels[1]]
longitude_array[3,:,:] = cornerarray
#fakeimage is then floated as an average of all the latitudes
fakeimage = sum(latitude_array)/len(latitude_array)
fakeimage = np.abs(fakeimage)
offlimb2 = np.where((fakeimage <= 0.))
fakeimage = np.cos((fakeimage*2))
fakeimage = fakeimage - np.amin(fakeimage)
fakeimage = fakeimage/np.amax(fakeimage) * 0.9
fakeimage = fakeimage + 0.1
fakeimage[offlimb2] = 0.
fakeimage = np.rot90(fakeimage, k=2)
los_coor = fakeimage

#%%
#"Effect los with smoothing"
a2 = (seeing/plate_scale)/(2*(np.sqrt(2*np.log(2))))
g_arraysize = a2 * 9
gaussian = np.zeros((int(g_arraysize), int(g_arraysize)))
a0 = 1.
a1 = g_arraysize/2.
g_array = np.arange(int(g_arraysize))           #To loop over
for x in g_array:
    for y in g_array:
        z = (x - a1)/a2
        z2 = (y - a1)/a2
        gaussian[x,y] = a0 * np.exp(-z**2 / 2) * np.exp(-z2**2 / 2)
#Tom then has a line which I don't think I can do in python (yet, anyway...)
los_coor = los_coor/los_coor[cox, coy]          #Maybe this is the line?
#It's los_coor=los_coor/los_coor(cox, coy)>0
#Next step is to convolve
#%%
los_coor = signal.convolve(los_coor, gaussian, mode='same')
#And that's it from makeimage3

#%%

"""
back to auaura_mapcoords_width
"""
#Tom calls the output from makeimage3 (i.e., los_coor) fakek in this script
fakek = los_coor
fsz = len(fakek)
#fakeimage = fakek[0:int(fsz/2),:]
#I now need the los variable which is from makeimage(1)

#%%
"""
makeimage1
"""
#Float the variables
eqwidth, sel, equatorkm, plate_scale = sat_eqwidth, sat_seangle, equatorkm, plate_scale
limbdist = (0.5 * eqwidth)/plate_scale
losflattening = flattening * (1. - np.sin(np.deg2rad(sel)))
eq_po_ratio = 1. - losflattening
po_eq_ratio = 1./eq_po_ratio
#"Add LOS effect"
im2 = dist_ellipse((naxis1, naxis1), cox, coy, po_eq_ratio, pa=90)
R2 = (R1 + 1000.)/R1
R3 = (R1 + 3000.)/R1
#Some conditions floated
rr1 = np.where(((im2/limbdist) < 1.))
rr2 = np.where(((im2/limbdist) >= 1.) & ((im2/limbdist) < R2))
rr3 = np.where(((im2/limbdist) >= R2) & ((im2/limbdist) < R3))
los_coor = (im2/limbdist)
los_coor[:] = 0
if len(rr1[0]) > 0:
    los_coor[rr1] = (np.sqrt(R3**2 - (im2/limbdist)[rr1]**2) - np.sqrt(R2**2 - (im2/limbdist)[rr1]**2))/(R3 - R2)
if len(rr2[0]) > 0:
    los_coor[rr2] = (np.sqrt(R3**2 - (im2/limbdist)[rr2]**2) - np.sqrt(R2**2 - (im2/limbdist)[rr2]**2))/(R3 - R2)       #;;;change from 2 to 1. because absorption reduces los effect - still not enough, but helps a little
if len(rr3[0]) > 0:
    los_coor[rr3] = 2. * (np.sqrt(R3**2 - (im2/limbdist)[rr3]**2))
#%%
#Next step is to multiply the los_coor by fakeimage
los_coor = los_coor * fakeimage
#%%
#"Effect los with smoothing"
a2 = (seeing/plate_scale)/(2*(np.sqrt(2*np.log(2))))
g_arraysize = a2 * 9
gaussian = np.zeros((int(g_arraysize), int(g_arraysize)))
a0 = 1.
a1 = g_arraysize/2.
g_array = np.arange(int(g_arraysize))           #To loop over
for x in g_array:
    for y in g_array:
        z = (x - a1)/a2
        z2 = (y - a1)/a2
        gaussian[x,y] = a0 * np.exp(-z**2 / 2) * np.exp(-z2**2 / 2)
        
#%%
#Tom then has a line which I don't think I can do in python (yet, anyway...)
los_coor = los_coor/los_coor[cox, coy]          #Maybe this is the line?
#It's los_coor=los_coor/los_coor(cox, coy)>0
#Next step is to convolve
los_coor = signal.convolve(los_coor, gaussian, mode='same')
los_coor = los_coor/los_coor[100,100]
#And that's the output from makeimage
#%%
"""
last return to mapcoords_width
"""
los = los_coor          #From makeimage1, above
fakeimage = fakek[:,:]
losimage = los[:,:]
latimage = latarray[:,:]
lonimage = lonarray[:,:]
fsz = len(fakeimage)
#Slit stuff next
slitwidth = slit_width   #0.288 (Day 5, 6, 7) or 0.432 (Day 1, 2, 3, 4)
iswidth = slitwidth/0.1445
fm2 = np.zeros((200,200))
fm2[:,:] = fakeimage
fakeimage = fm2
fsz1 = len(fakeimage[0])
fsz2 = len(fakeimage)
fm2 = np.zeros((200,200))
fm2[:,:] = losimage
losimage = fm2
fm2 = np.zeros((200,200))
fm2[:,:] = latimage
latimage = fm2
fm2 = np.zeros((200,200))
fm2[:,:] = lonimage
lonimage = fm2
#fakeslitarray stuff next
fakeslitarray = np.zeros((int(fsz2 - iswidth), fsz1))
losslitarray = np.zeros((int(fsz2 - iswidth), fsz1))
latslitarray = np.zeros((int(fsz2 - iswidth), fsz1))
lonslitarray = np.zeros((int(fsz2 - iswidth), fsz1))
#Next for loop
fsz2_iswidth_loop = np.arange(int(fsz2 - iswidth))
for y in fsz2_iswidth_loop:
    fakeslitarray[y,:] = sum(fakeimage[y:((y+int(iswidth))-1),:])   #Not entirely sure what this line is meant to do
    losslitarray[y,:] = sum(losimage[y:((y+int(iswidth))-1),:])/iswidth         #Same with this line too...
    latslitarray[y,:] = latimage[y+int(iswidth/2.),:]
    lonslitarray[y,:] = lonimage[y+int(iswidth/2.),:]
#That for loop was proper dodgy... Should really check with Tom what it's meant to do!
#Ok, so, the very last part of that script now and it looks like there'll be some congrid work going on here
#All bounds are (1500 - 150 * 4 = 900) and (1500 + 150 * 4 - 1 = ~2100) which together keep the central 1000 to 2000 and leave all the edges. Also, x-direction shrink factor is 800/(2100-900), i.e., desired/current
#%%
fsa2 = congrid.zoom(fakeslitarray[:,:], [1, 1], order=3, mode='nearest')       #0.18333333333333332 for 220, 0.25 for 300, 0.666666666 for 800
fsa3 = fsa2
for i in fsz2_iswidth_loop:
    fsa3[i,:] = fsa3[i,:]/np.amax(fsa3[i,:])
fake1 = fsa3
fsa2 = congrid.zoom(losslitarray[:,:], [1, 1], order=3, mode='nearest')       #0.18333333333333332
fsa3 = fsa2
for i in fsz2_iswidth_loop:
    fsa3[i,:] = fsa3[i,:]
los1 = fsa3
fsa2 = congrid.zoom(latslitarray[:,:], [1, 1], order=3, mode='nearest')       #0.18333333333333332
fsa3 = fsa2
for i in fsz2_iswidth_loop:
    fsa3[i,:] = fsa3[i,:]
lat1 = fsa3
fsa2 = congrid.zoom(lonslitarray[:,:], [1, 1], order=3, mode='nearest')       #0.18333333333333332
fsa3 = fsa2
for i in fsz2_iswidth_loop:
    fsa3[i,:] = fsa3[i,:]
lon1 = fsa3
#Visualise these
Keck_spectra_viewer(np.flip(fake1), 'fake1')
Keck_spectra_viewer(np.flip(los1), 'los1')
Keck_spectra_viewer(np.flip(lat1), 'lat1')
Keck_spectra_viewer(np.flip(lon1), 'lon1')

    # return fake1, los1, lat1, lon1
    
#%% Pull out the figures for Scam to get a rough idea of where Uranus is
from astropy.io import fits 
import matplotlib.pyplot as plt
import warnings

image_file1 = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/sci/NC.20141013.37870.fits.gz'
image_file2 = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/sci/NC.20141013.37883.fits.gz'
image_file3 = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/sci/NC.20141013.38321.fits.gz'
image_file4 = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/sci/NC.20141013.38337.fits.gz'

image_data1 = fits.getdata(image_file1, ext=0)
image_data2 = fits.getdata(image_file2, ext=0) #Focuses specifically on the array data contained within the fits file
image_data3 = fits.getdata(image_file3, ext=0)
image_data4 = fits.getdata(image_file4, ext=0)
warnings.filterwarnings('ignore', 'The following header keyword is invalid or follows an unrecognized non-standard convention',)
hdul = fits.open(image_file3, ignore_missing_end=True)
hdr = hdul[0].header

#print(hdr)

SCAM1 = (image_data1 - image_data2)/2
SCAM2 = (image_data3 - image_data4)/2

# plt.figure()
# plt.imshow(SCAM1, cmap='gist_gray', vmax = 250, vmin = -150)

import scipy.ndimage as f
SCAM3 = f.rotate(SCAM2, angle = 6.8100)
SCAM3 = f.shift(SCAM3, shift=(-58,-67)) #41 48
SCAM3 = f.zoom(SCAM3, zoom=0.183/0.1445)
# SCAM3 = f.rotate(SCAM3, angle = 8.1)
# plt.figure()
# plt.imshow(SCAM3[0:200,0:200], cmap='gist_gray', vmax = 250, vmin = -250)
# plt.imshow(np.flipud(lonslitarray), cmap='gist_gray', vmax = 10, vmin = -10, alpha = 0.5)
# plt.hlines((99,100), xmin=0, xmax=200)
# plt.xlim(0,200)

#Lets shift the positive and negative on top of each other
SCAM4 = f.rotate(SCAM2, angle = 8.1)
SCAM4 = f.shift(SCAM4, shift=(-7,-26))
SCAM4 = f.zoom(SCAM4, zoom=0.183/0.1445)

plt.figure()
plt.imshow(SCAM3[0:200,0:200], cmap='gist_gray', vmax = 250, vmin = -250, alpha = 0.5)
#plt.imshow(-1*SCAM4[0:200,0:200], vmax = 250, vmin = -250, alpha = 0.5)
plt.imshow(lonslitarray, cmap='gist_gray', vmax = 10, vmin = -10, alpha = 0.5)
plt.hlines((101,102), xmin=0, xmax=200)
plt.xlim(0,200)

#Hence hlines show where obs are from
LATS = latslitarray[101:103,83:116]
LONGS = np.hstack((lonslitarray[101,83:115], ((180 - lonslitarray[101,83]) + 180)))
LONGS1 = np.hstack((lonslitarray[102,83:115], ((180 - lonslitarray[102,83]) + 180)))
LONGS_Fin = np.vstack((LONGS, LONGS1))

# plt.figure()
# plt.imshow(SCAM3[0:200,0:200], cmap='gist_gray', vmax = 250, vmin = -250)
# plt.imshow(np.flipud(lonslitarray), cmap='gist_gray', vmax = 10, vmin = -10, alpha = 0.5)
# plt.hlines((99,100), xmin=0, xmax=200)
# plt.xlim(0,200)
# Uranus_array = lonslitarray + SCAM2[0:200,0:200]

# plt.figure()
# # plt.imshow(lonslitarray, cmap='gist_gray', vmax = )

# warnings.filterwarnings('ignore', 'The following header keyword is invalid or follows an unrecognized non-standard convention',)

# #Find the orientation of the slit
# for n in range(72): #We use this list to create a list which holds all the data from Order19
#     num = n + 45
#     if num < 100:
#         image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s00' + str(num) + '.fits'
#         hdul = fits.open(image_filei, ignore_missing_end=True)
#         hdr = hdul[0].header
#         SLIT_ANG.append(hdr['SLITANG'])
#         #print(hdr['UTC'])
#     elif num >= 100:
#         image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/KECK 13OCT14/oct13/oct13s0' + str(num) + '.fits'
#         hdul = fits.open(image_filei, ignore_missing_end=True)
#         hdr = hdul[0].header
#         SLIT_ANG.append(hdr['SLITANG'])
#         #print(hdr['UTC'])
#     else:
#         pass

#%%
#Saving lats and longs along with los to assist with los brightening
# np.save('LatitudeJup1116.npy', lat1)
# np.save('LongitudeJup1116.npy', lon1)
# np.save('LOSJup1116.npy', los1)

# #%% Now we work out the 11 slit positions on this model
# TIMES = np.load('OBS_Times.npy')
# Images = []
# TIMES_Kyle = []
# from astropy.io import fits

# for m in range(340):
#     if m < 9:
#         image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Kyle/ike.2016B041.161011.uranus.0000' + str(m+1) + '.a.fits'
#         image_datai = fits.getdata(image_filei, ext=0)
#         hdu = fits.open(image_filei)
#         hdr = hdu[0].header
#         TIMES_Kyle.append(hdr['TIME_OBS'])
#         IRTF_Data = image_datai
#         Images.append(IRTF_Data)
#     elif m < 99:
#         image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Kyle/ike.2016B041.161011.uranus.000' + str(m+1) + '.a.fits'
#         image_datai = fits.getdata(image_filei, ext=0)
#         hdu = fits.open(image_filei)
#         hdr = hdu[0].header
#         TIMES_Kyle.append(hdr['TIME_OBS'])
#         IRTF_Data = image_datai
#         Images.append(IRTF_Data)  
#     else:
#         image_filei = 'C:/Users/snowy/OneDrive/Documents/Python work/IRTF 11OCT16/Kyle/ike.2016B041.161011.uranus.00' + str(m+1) + '.a.fits'
#         image_datai = fits.getdata(image_filei, ext=0)
#         hdu = fits.open(image_filei)
#         hdr = hdu[0].header
#         TIMES_Kyle.append(hdr['TIME_OBS'])
#         IRTF_Data = image_datai
#         Images.append(IRTF_Data)

# #%%
# #We need the 343 to be sorted into 11 bins dependant on the stopping times and start times of the 11 sets above so gates would be 1, 4, 7, 10, etc.
# #Convert both to seconds 
# TIMES_Seconds_Kyle = []
# TIMES_Seconds = []
# for m in range(340):
#     Hours = int(TIMES_Kyle[m][0:2])*60*60
#     Minutes = int(TIMES_Kyle[m][3:5])*60
#     Seconds = float(TIMES_Kyle[m][6:])
#     TIMES_Seconds_Kyle.append(Hours+Minutes+Seconds)

# for m in range(45):
#     Hours = int(TIMES[m][0:2])*60*60
#     Minutes = int(TIMES[m][3:5])*60
#     Seconds = float(TIMES[m][6:])
#     TIMES_Seconds.append(Hours+Minutes+Seconds)
    

# #And now we can start gatekeeping!
# SetP = []
# Set1 = []
# Set2 = []
# Set3 = []
# Set4 = []
# Set5 = []
# Set6 = []
# Set7 = []
# Set8 = []
# Set9 = []
# Set10 = []
# Set11 = []
# SetA = []

# for m in range(340):
#     if TIMES_Seconds_Kyle[m] <= TIMES_Seconds[0]:
#         SetP.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[0] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[4]:
#         Set1.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[4] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[8]:
#         Set2.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[8] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[12]:
#         Set3.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[12] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[16]:
#         Set4.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[16] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[20]:
#         Set5.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[20] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[24]:
#         Set6.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[24] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[28]:
#         Set7.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[28] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[32]:
#         Set8.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[32] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[36]:
#         Set9.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[36] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[40]:
#         Set10.append(m)
#     elif TIMES_Seconds_Kyle[m] > TIMES_Seconds[40] and TIMES_Seconds_Kyle[m] <= TIMES_Seconds[44]:
#         Set11.append(m)
#     else:
#         SetA.append(m)
        
# #%% Lets start with Set1
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set1)):
#     Kyle_View_Set1 += Images[Set1[m]]

# import scipy.ndimage as f
# import matplotlib.pyplot as plt
    
# #Lets f.shift the view by approx amount it should be and then get the final image, put the amoutn Uranus should be and then work out how to take that over to the model
# #To confirm we have overlap lets do everything above the noise is 1 and everything below is 0

# #Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 7200] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,73.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 7200] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# #Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = (Kyle_View_Set1 + Kyle_View2_Set1)/2

# #%% Now we need measure how the image is to work out the middle
# import statistics
# def gauss_fit(x, a0, a1, a2, a3, a4, a5): # First write a guassian function credit to pen and pants IDL's Gaussfit in Python
#     z = (x - a1) / a2
#     y = a0 * np.exp(-z**2 / a2) + a3 + a4 * x + a5 * x**2
#     return y

# Sample_Set1_Kyle = Final_Set1_Kyle[255:315,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[6:37,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line) # 2.061685 pixels of slit
# #%%
# Blank = np.zeros((160,43))
# Blank2 = np.zeros((200,157))
# Set1_Figure_Final = np.concatenate((Set1_Figure, Blank), axis=0)
# Set1_Figure_Final = np.concatenate((Set1_Figure_Final, Blank2), axis=1)

# Set1_Figure_Final = f.shift(Set1_Figure_Final, shift=(80, 78.5), order=3, mode='wrap')
# Set_Figure_Finalrot = np.rot90(Set1_Figure_Final)


# # If we know the slit is 0.375 and the pixel scale is approx 0.12 then the slit is 

# fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
# ax.imshow(Set1_Figure_Final+Set_Figure_Finalrot, cmap='gist_gray')

# LONGS_VIEW = f.shift(lonslitarray[:,:], shift=(1,1), order=3, mode='wrap')
# LATS_VIEW = f.shift(latslitarray[:,:], shift=(1,1), order=3, mode='wrap')
# # LATS_VIEW[LATS_VIEW < 1e5 and LATS_VIEW > -1e5] = 0

# # MASK = LONGS_VIEW
# # MASK[MASK < 1] = 0
# # MASK[MASK > -1] = 1
# plt.figure()
# plt.imshow(LONGS_VIEW)
# plt.imshow(Set1_Figure_Final+Set_Figure_Finalrot, cmap='gist_gray', alpha=0.75)
# # plt.hlines((99,100,101), 0, 200, color='k')
# # plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# plt.figure()
# plt.imshow(LATS_VIEW)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# #So we know roughly where the lats and longitudes are for this first set
# Set1_LATS = LATS_VIEW[99:102,75:125] #Confirm where the planet starts and stops
# Set1_LONGS = LONGS_VIEW[99:102, 75:125]
# # print(LONGS_VIEW[99:101, 75:125])

# #%% Lets sort the final datapoint
# from scipy.optimize import curve_fit
# def poly_func(x, a , b, c):
#     y = a* x**2 + b*x + c
#     return y

# # gmodel = Model(poly_func)
# # result_L1 = gmodel.fit(Set1_LATS[0,12:38], x=Set1_LONGS[0,12:38], a=0, b=1, c=1, d=0)
# # pL1 = SimpleNamespace(**result_L1.best_values)
# # eL1 = np.sqrt(np.diag(result_L1.covar))

# # popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[0,12:38], Set1_LATS[0,12:38], p0=[4.48381659e-03, -1.61385393e+00,  1.09327948e+02])
# plt.figure()
# plt.plot(Set1_LONGS[0,12:39], Set1_LATS[0,12:39], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set1_LONGS[0,12]) + 180
# plt.plot(EXTRA, Set1_LATS[0,38], 'ro')

# Set1_LATS_TOP = Set1_LATS[0,12:39]
# Set1_LONGS_TOP = np.append(Set1_LONGS[0,12:38], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set1_LATS_TOP, Set1_LONGS_TOP, 'bo')

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[1,12:38], Set1_LATS[1,12:38], p0=[4.43377329e-03, -1.59934232e+00,  1.03391336e+02])
# plt.figure()
# plt.plot(Set1_LONGS[1,12:39], Set1_LATS[1,12:39], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set1_LONGS[1,12]) + 180
# plt.plot(EXTRA, Set1_LATS[1,38], 'ro')

# Set1_LATS_MID = Set1_LATS[1,12:39]
# Set1_LONGS_MID = np.append(Set1_LONGS[1,12:38], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set1_LATS_MID, Set1_LONGS_MID, 'bo')

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[2,13:37], Set1_LATS[2,13:37], p0=[4.31327778e-03, -1.55802499e+00,  9.49789733e+01])
# plt.figure()
# plt.plot(Set1_LONGS[2,13:38], Set1_LATS[2,13:38], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set1_LONGS[2,13]) + 180
# plt.plot(EXTRA, Set1_LATS[2,37], 'ro')

# Set1_LATS_BOT = Set1_LATS[2,13:38]
# Set1_LONGS_BOT = np.append(Set1_LONGS[2,13:37], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set1_LATS_BOT, Set1_LONGS_BOT, 'bo')

# #%% NORTH
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set5)):
#     Kyle_View_Set1 += Images[Set5[m]]
    
# #Lets f.shift the view by approx amount it should be and then get the final image, put the amoutn Uranus should be and then work out how to take that over to the model
# #To confirm we have overlap lets do everything above the noise is 1 and everything below is 0

# # Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 4e4] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,74.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 4e4] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = (Kyle_View_Set1 + Kyle_View2_Set1)/2

# plt.figure()
# plt.imshow(Final_Set1_Kyle, cmap='gist_gray', vmax = 1e6, vmin = 0)

# Sample_Set1_Kyle = Final_Set1_Kyle[255:315,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[13:37,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line)

# Blank = np.zeros((160,43))
# Blank2 = np.zeros((200,157))
# Set1_Figure_Final = np.concatenate((Set1_Figure, Blank), axis=0)
# Set1_Figure_Final = np.concatenate((Set1_Figure_Final, Blank2), axis=1)

# Set1_Figure_Final = f.shift(Set1_Figure_Final, shift=(80, 78.5), order=3, mode='wrap')
# Set_Figure_Finalrot = np.rot90(Set1_Figure_Final)


# # If we know the slit is 0.375 and the pixel scale is approx 0.12 then the slit is 

# fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
# ax.imshow(Set1_Figure_Final, cmap='gist_gray')

# LONGS_VIEW = f.shift(lonslitarray[:,899:1100], shift=(6,1), order=3, mode='wrap')
# LATS_VIEW = f.shift(latslitarray[:,899:1100], shift=(6,1), order=3, mode='wrap')
# # LATS_VIEW[LATS_VIEW < 1e5 and LATS_VIEW > -1e5] = 0

# # MASK = LONGS_VIEW
# # MASK[MASK < 1] = 0
# # MASK[MASK > -1] = 1
# plt.figure()
# plt.imshow(LONGS_VIEW)
# #plt.imshow(Set1_Figure_Final, cmap='gist_gray', alpha=0.75)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# plt.figure()
# plt.imshow(LATS_VIEW)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# Set5_LATS = LATS_VIEW[99:102, 75:125]
# Set5_LONGS = LONGS_VIEW[99:102, 75:125]

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set5_LONGS[0,14:37], Set5_LATS[0,14:37], p0=[4.48381659e-03, -1.61385393e+00,  1.09327948e+02])
# plt.figure()
# plt.plot(Set5_LONGS[0,14:39], Set5_LATS[0,14:39], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set5_LONGS[0,15]) + 180
# EXTRA2 = (180 - Set5_LONGS[0,14]) + 180
# plt.plot(EXTRA, Set5_LATS[0,37], 'ro')
# plt.plot(EXTRA2, Set5_LATS[0,38], 'ro')

# Set5_LATS_TOP = Set5_LATS[0,14:39]
# Set5_LONGS_TOP = np.append(Set5_LONGS[0,14:37], EXTRA)
# Set5_LONGS_TOP = np.append(Set5_LONGS_TOP, EXTRA2)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set5_LATS_TOP, Set5_LONGS_TOP, 'bo')

# #%%
# #popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[1,12:38], Set1_LATS[1,12:38], p0=[4.43377329e-03, -1.59934232e+00,  1.03391336e+02])
# plt.figure()
# plt.plot(Set5_LONGS[1,14:39], Set5_LATS[1,14:39], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set5_LONGS[1,14]) + 180
# plt.plot(EXTRA, Set5_LATS[1,38], 'ro')

# Set5_LATS_MID = Set5_LATS[1,14:39]
# Set5_LONGS_MID = np.append(Set5_LONGS[1,14:38], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set5_LATS_MID, Set5_LONGS_MID, 'bo')

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[2,13:37], Set1_LATS[2,13:37], p0=[4.31327778e-03, -1.55802499e+00,  9.49789733e+01])
# plt.figure()
# plt.plot(Set5_LONGS[2,14:38], Set5_LATS[2,14:38], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set5_LONGS[2,14]) + 180
# plt.plot(EXTRA, Set5_LATS[2,38], 'ro')

# Set5_LATS_BOT = Set5_LATS[2,14:39]
# Set5_LONGS_BOT = np.append(Set5_LONGS[2,14:38], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set5_LATS_BOT, Set5_LONGS_BOT, 'bo')


# #%% Set 6 EQUATOR
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set6)):
#     Kyle_View_Set1 += Images[Set6[m]]

# # Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 1.5e4] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,74.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 1.5e4] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = (Kyle_View_Set1 + Kyle_View2_Set1)/2

# plt.figure()
# plt.imshow(Final_Set1_Kyle, cmap='gist_gray')

# Sample_Set1_Kyle = Final_Set1_Kyle[255:315,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[12:31,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line)

# Blank = np.zeros((160,43))
# Blank2 = np.zeros((200,157))
# Set1_Figure_Final = np.concatenate((Set1_Figure, Blank), axis=0)
# Set1_Figure_Final = np.concatenate((Set1_Figure_Final, Blank2), axis=1)

# Set1_Figure_Final = f.shift(Set1_Figure_Final, shift=(80, 78.5), order=3, mode='wrap')
# Set_Figure_Finalrot = np.rot90(Set1_Figure_Final)

# fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
# ax.imshow(Set1_Figure_Final+Set_Figure_Finalrot, cmap='gist_gray')

# LONGS_VIEW = f.shift(lonslitarray[:,899:1100], shift=(1,0), order=3, mode='wrap')
# LATS_VIEW = f.shift(latslitarray[:,899:1100], shift=(1,0), order=3, mode='wrap')
# # LATS_VIEW[LATS_VIEW < 1e5 and LATS_VIEW > -1e5] = 0

# # MASK = LONGS_VIEW
# # MASK[MASK < 1] = 0
# # MASK[MASK > -1] = 1
# plt.figure()
# plt.imshow(LONGS_VIEW)
# plt.imshow(Set1_Figure_Final+Set_Figure_Finalrot, cmap='gist_gray', alpha=0.75)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# plt.figure()
# plt.imshow(LATS_VIEW)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# #So we know roughly where the lats and longitudes are for this first set
# Set6_LATS = LATS_VIEW[99:102, 75:125]
# Set6_LONGS = LONGS_VIEW[99:102, 75:125]

# #%% SOUTH
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set8)):
#     Kyle_View_Set1 += Images[Set8[m]]
# # FINISH Set 3 and 4 along with 8 before carrying out the other method.

# # Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 1.5e4] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,74.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 1.5e4] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = Kyle_View_Set1

# plt.figure()
# plt.imshow(Final_Set1_Kyle, cmap='gist_gray')

# Sample_Set1_Kyle = Final_Set1_Kyle[255:315,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[6:37,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line)

# Blank = np.zeros((160,43))
# Blank2 = np.zeros((200,157))
# Set1_Figure_Final = np.concatenate((Set1_Figure, Blank), axis=0)
# Set1_Figure_Final = np.concatenate((Set1_Figure_Final, Blank2), axis=1)

# Set1_Figure_Final = f.shift(Set1_Figure_Final, shift=(80, 78.5), order=3, mode='wrap')
# #Set_Figure_Finalrot = np.rot90(Set1_Figure_Final)

# fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
# ax.imshow(Set1_Figure_Final, cmap='gist_gray')

# LONGS_VIEW = f.shift(lonslitarray[:,899:1100], shift=(-2.5,1), order=3, mode='wrap')
# LATS_VIEW = f.shift(latslitarray[:,899:1100], shift=(-2.5,1), order=3, mode='wrap')
# # LATS_VIEW[LATS_VIEW < 1e5 and LATS_VIEW > -1e5] = 0

# # MASK = LONGS_VIEW
# # MASK[MASK < 1] = 0
# # MASK[MASK > -1] = 1
# plt.figure()
# plt.imshow(LONGS_VIEW)
# plt.imshow(Set1_Figure_Final, cmap='gist_gray', alpha=0.75)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# plt.figure()
# plt.imshow(LATS_VIEW)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# #So we know roughly where the lats and longitudes are for this first set
# Set8_LATS = LATS_VIEW[99:102, 75:125]
# Set8_LONGS = LONGS_VIEW[99:102, 75:125]

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set5_LONGS[0,14:37], Set5_LATS[0,14:37], p0=[4.48381659e-03, -1.61385393e+00,  1.09327948e+02])
# plt.figure()
# plt.plot(Set8_LONGS[0,14:39], Set8_LATS[0,14:39], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set8_LONGS[0,14]) + 180
# plt.plot(EXTRA, Set8_LATS[0,38], 'ro')

# Set8_LATS_TOP = Set8_LATS[0,14:39]
# Set8_LONGS_TOP = np.append(Set8_LONGS[0,14:38], EXTRA)

# plt.figure()
# plt.plot(Set8_LATS_TOP, Set8_LONGS_TOP, 'bo')

# #%%
# #popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[1,12:38], Set1_LATS[1,12:38], p0=[4.43377329e-03, -1.59934232e+00,  1.03391336e+02])
# plt.figure()
# plt.plot(Set8_LONGS[1,15:38], Set8_LATS[1,15:38], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set8_LONGS[1,15]) + 180
# plt.plot(EXTRA, Set8_LATS[1,37], 'ro')

# Set8_LATS_MID = Set8_LATS[1,15:38]
# Set8_LONGS_MID = np.append(Set8_LONGS[1,15:37], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set8_LATS_MID, Set8_LONGS_MID, 'bo')

# #%%
# # popt1, pcov1 = curve_fit(poly_func, Set1_LONGS[2,13:37], Set1_LATS[2,13:37], p0=[4.31327778e-03, -1.55802499e+00,  9.49789733e+01])
# plt.figure()
# plt.plot(Set8_LONGS[2,15:38], Set8_LATS[2,15:38], 'bo')
# #plt.plot(xdata, Att1, 'g--')
# #Fitdata = lin_func(A)
# # plt.plot(np.arange(360), poly_func(np.arange(360), *popt1), 'g--')
# EXTRA = (180 - Set8_LONGS[2,14]) + 180
# plt.plot(EXTRA, Set8_LATS[2,37], 'ro')

# Set8_LATS_BOT = Set8_LATS[2,15:38]
# Set8_LONGS_BOT = np.append(Set8_LONGS[2,15:37], EXTRA)
# # Set1_LONGS_TOP.extend(EXTRA)

# plt.figure()
# plt.plot(Set8_LATS_BOT, Set8_LONGS_BOT, 'bo')

# #%% Three positions done Find what the last position is
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set11)):
#     if Set11[m] < 257:
#         Kyle_View_Set1 += Images[Set11[m]]
#     else:
#         pass

# Kyle_View2_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set11)):
#     if Set11[m] >= 257 and Set11[m] < 263:
#         Kyle_View2_Set1 += Images[Set11[m]]
#     else:
#         pass

# Kyle_View3_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set11)):
#     if Set11[m] >= 263:
#         Kyle_View3_Set1 += Images[Set11[m]]
#     else:
#         pass
    
# # Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 1.5e4] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# #Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,74.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 1.5e4] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# #Final_Set1_Kyle = Kyle_View_Set1

# # plt.figure()
# # plt.imshow(Kyle_View_Set1, cmap='gist_gray')

# # plt.figure()
# # plt.imshow(Kyle_View2_Set1, cmap='gist_gray')

# # plt.figure()
# # plt.imshow(Kyle_View3_Set1, cmap='gist_gray')

# #%% Complete Set 3 and 4
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set3)):
#     Kyle_View_Set1 += Images[Set3[m]]

# Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 5000] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,73.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 5000] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = (Kyle_View_Set1 + Kyle_View2_Set1)/2

# plt.figure()
# plt.imshow(Final_Set1_Kyle, cmap='gist_gray')

# Sample_Set1_Kyle = Final_Set1_Kyle[255:305,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[6:37,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line)

# #%% SOUTH
# Blank = np.zeros((167,43))
# Blank2 = np.zeros((200,157))
# Set1_Figure_Final = np.concatenate((Set1_Figure, Blank), axis=0)
# Set1_Figure_Final = np.concatenate((Set1_Figure_Final, Blank2), axis=1)

# Set1_Figure_Final = f.shift(Set1_Figure_Final, shift=(80, 78.5), order=3, mode='wrap')
# #Set_Figure_Finalrot = np.rot90(Set1_Figure_Final)

# fig, ax = plt.subplots(subplot_kw={'aspect':'equal'})
# ax.imshow(Set1_Figure_Final, cmap='gist_gray')

# LONGS_VIEW = f.shift(lonslitarray[:,899:1100], shift=(-2.5,1), order=3, mode='wrap')
# LATS_VIEW = f.shift(latslitarray[:,899:1100], shift=(-2.5,1), order=3, mode='wrap')
# # LATS_VIEW[LATS_VIEW < 1e5 and LATS_VIEW > -1e5] = 0

# # MASK = LONGS_VIEW
# # MASK[MASK < 1] = 0
# # MASK[MASK > -1] = 1
# plt.figure()
# plt.imshow(LONGS_VIEW)
# plt.imshow(Set1_Figure_Final, cmap='gist_gray', alpha=0.75)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# plt.figure()
# plt.imshow(LATS_VIEW)
# plt.hlines((99,100,101), 0, 200, color='k')
# plt.vlines((99,100,101), 0, 200, color='k')
# plt.ylim(200,0)

# Set3_LATS = LATS_VIEW[99:102, 75:125]
# Set3_LONGS = LONGS_VIEW[99:102, 75:125]

# Set4_LATS = LATS_VIEW[99:102, 75:125]
# Set4_LONGS = LONGS_VIEW[99:102, 75:125]

# #%%
# Kyle_View_Set1 = np.zeros(np.shape(Images[0]))
# for m in range(len(Set4)):
#     Kyle_View_Set1 += Images[Set4[m]]

# Mask_Set1 = Kyle_View_Set1
# # Mask_Set1[Mask_Set1 < np.nanmean(Kyle_View_Set1)] = 0
# # Mask_Set1[Mask_Set1 > 5000] = 1
# # Mask_Set1[Mask_Set1 > 2] = 0

# Kyle_View2_Set1 = f.shift(Kyle_View_Set1, shift=(0,73.5), order=3, mode='wrap')
# # Mask_Set2 = Kyle_View2_Set1
# # Mask_Set2[Mask_Set2 < np.nanmean(Kyle_View2_Set1)] = 0
# # Mask_Set2[Mask_Set2 > 5000] = 1
# # Mask_Set2[Mask_Set2 > 2] = 0

# # Final_Mask = (Mask_Set1 + Mask_Set2)/2
# Final_Set1_Kyle = (Kyle_View_Set1 + Kyle_View2_Set1)/2

# plt.figure()
# plt.imshow(Final_Set1_Kyle, cmap='gist_gray')

# Sample_Set1_Kyle = Final_Set1_Kyle[255:305,240:305]
# Set1_Figure = f.zoom(Sample_Set1_Kyle, zoom = 0.6597392) #Use this figure to find out where the line intersects and then compare with the model!!! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# Line = []
# #Lets find the position of the line in this space
# for m in range(43):
#     a1_pointP = np.where(Set1_Figure == np.nanmin(Set1_Figure[6:37,m]))
#     Line.append(a1_pointP[0][0])

# Pos_Line = statistics.mode(Line)

# #%% Now we pull all the coordinates into an array so it can be pulled out easily
# NORTH_SLIT = []
# SOUTH_SLIT = []
# EQUATOR_SLIT = []

# NORTH_SLIT.append(Set5_LATS_TOP)
# NORTH_SLIT.append(Set5_LONGS_TOP)
# NORTH_SLIT.append(Set5_LATS_MID)
# NORTH_SLIT.append(Set5_LONGS_MID)
# NORTH_SLIT.append(Set5_LATS_BOT)
# NORTH_SLIT.append(Set5_LONGS_BOT)

# EQUATOR_SLIT.append(Set1_LATS_TOP)
# EQUATOR_SLIT.append(Set1_LONGS_TOP)
# EQUATOR_SLIT.append(Set1_LATS_MID)
# EQUATOR_SLIT.append(Set1_LONGS_MID)
# EQUATOR_SLIT.append(Set1_LATS_BOT)
# EQUATOR_SLIT.append(Set1_LONGS_BOT)

# SOUTH_SLIT.append(Set8_LATS_TOP)
# SOUTH_SLIT.append(Set8_LONGS_TOP)
# SOUTH_SLIT.append(Set8_LATS_MID)
# SOUTH_SLIT.append(Set8_LONGS_MID)
# SOUTH_SLIT.append(Set8_LATS_BOT)
# SOUTH_SLIT.append(Set8_LONGS_BOT)

# np.save('NORTH_SLIT_CONFIG.npy', NORTH_SLIT)
# np.save('EQUATOR_SLIT_CONFIG.npy', EQUATOR_SLIT)
#np.save('SOUTH_SLIT_CONFIG.npy', SOUTH_SLIT)