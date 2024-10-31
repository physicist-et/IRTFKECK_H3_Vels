# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 14:56:20 2022

@author: snowy
"""
import numpy as np

uranus_seangleI = (-35.773649 + -35.784010)/2
#stretch yy to become a sphere
flattening =0.0229
pix_to
uranus_seangleI_Rads = (uranus_seangleI*np.pi)/180
losflattening=flattening*(1-np.sin(uranus_seangleK_Rads))
eq_po_ratioK=1-losflattening

y_IRTF = [11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5, -6, -7, -8, -9, -10, -11]

y_Keckadj = y_Keck/eq_po_ratioK

losflattening=flattening*(1-np.sin((uranus_seangleI/180)*np.pi))
eq_po_ratioI=1-losflattening

y1_IRTF = 0.89263
y2_IRTF = -0.89263

y1_IRTFadj = y1_IRTF/eq_po_ratioI
y2_IRTFadj = y2_IRTF/eq_po_ratioI