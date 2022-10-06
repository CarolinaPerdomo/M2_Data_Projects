import numpy as np
import math
import healpy as hp
import matplotlib.pyplot as plt
from astropy import units as u
from astropy.io import fits as pyfits
from astropy import wcs


freq=['100', '143','217','353','545','857'] #list of strings
freq2=np.array([100.,143.,217.,353.,545.,857.]) #array of floats
path = '/home/abeelen/Planck/maps_2015/'
angle = 9.7/60 #9.7 arcmin, but we need it in radians


#Loop that, for each map (for each frequency), smooth it (9.7 arcmin, the precision of the 100GHz one, the least precise), write a new smoothed map, and degrade it (512 pixels), and write a new smoothed and degraded map.
"""
for i in range(6):
    m = hp.read_map(path+'HFI_SkyMap_'+freq[i]+'_2048_R2.00_full.fits')
    smooth = hp.smoothing(m, fwhm=np.radians(angle))
    hp.write_map('/home/adafne/Project/smooth_'+freq[i]+'.fits', smooth)
    degrad  = hp.ud_grade(smooth,nside_out=512, order_in='RING', pess=True)
    hp.write_map("/home/adafne/Project/degrad"+freq[i]+".fits", degrad)
"""

#To open and look at a map :

"""
m = hp.read_map(path+'HFI_SkyMap_353_2048_R2.00_full.fits')
hp.mollview(m, title="353",norm='hist')
plt.show()
"""

#Let's zoom on the patch around the cluster PSZ2 G006.76+30.45 (number 22 in the list).

"""
for k in range(6):
    m = hp.read_map("/home/adafne/Project/degrad"+freq[k]+".fits")
    y, x = np.indices((256,256)) #creates two indices for our patch
    w_patch = wcs.WCS(naxis=2) #creates a wcs object indicating the link btw our patch and the sky
    w_patch.wcs.crpix = [128,128]  #pixel in the center
    w_patch.wcs.cdelt = np.array([-0.01, 0.01])   #scale
    w_patch.wcs.crval = [6.76,30.45]   #localisation in the sky
    w_patch.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]  #type of coordinates
    theta,phi = w_patch.wcs_pix2world(y, x, 1)  #convert to angular coordinates
    theta = theta * math.pi /180.  #conversion degrees radians
    phi = phi * math.pi /180.      #conversion
    i = hp.ang2pix(512, math.pi/2.-phi, theta)   #find pixels of the maps corresponding to these coordinates
    patch = np.zeros((256,256)) #creates the patch
    patch[x,y] = m[i]  #fill with values of the map
    hdu = pyfits.PrimaryHDU(header=w_patch.to_header())
    hdu.data = patch
    hdu.writeto('/home/adafne/Project/patch_'+freq[k]+'.fits')
"""

#Visualise our patch
"""
m = pyfits.getdata("/home/adafne/Project/patch_100.fits")
plt.imshow(m)
"""


#test zoom on Andromeda

"""
m = hp.read_map("/home/adafne/Project/degrad857.fits")
y, x = np.indices((256,256)) 
w_patch = wcs.WCS(naxis=2) 
w_patch.wcs.crpix = [128,128]
w_patch.wcs.cdelt = np.array([-0.05, 0.05])
w_patch.wcs.crval = [121.17,-21.57]
w_patch.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]
theta,phi = w_patch.wcs_pix2world(y, x, 1)  
theta = theta * math.pi /180.  
phi = phi * math.pi /180.
i = hp.ang2pix(512, math.pi/2.-phi, theta)   
patch = np.zeros((256,256))
patch[x,y] = m[i]
hdu = pyfits.PrimaryHDU(header=w_patch.to_header())
hdu.data = patch
hdu.writeto('/home/adafne/Project/andromeda9.fits')
m = pyfits.getdata("/home/adafne/Project/andromeda9.fits")
plt.imshow(m)
"""

#We have to apply some conversion given in the table page 14 of http://planck.caltech.edu/pub/2013results/Planck_2013_results_09.pdf. The two last carts should be divided respectively by 58.04 and 2.27. The six should be multiplied by the coefficients given in the last column.

"""
conversion = np.array([-0.248, -0.359, 5.152, 0.161, 0.069/58.04, 0.038/2.27]) 

for i in range(6):
    data,header = pyfits.getdata('/home/adafne/Project/patch_'+freq[i]+'.fits', header = True)
    data = data*conversion[i]
    hdu = pyfits.PrimaryHDU()
    hdu.data = data
    hdu.writeto('/home/adafne/Project/patch_convert_'+freq[i]+'.fits')
"""


# We want now to make the linear combination. Let's write the covariant matrix.

#First let's get the temperature.
"""
Temp = []

for k in range(6):
    T = pyfits.getdata('/home/adafne/Project/patch_convert_'+freq[k]+'.fits')
    Temp.append(T)

#We need the average for each map.

Temp_average = []

for k in range(6):
    Ta = np.mean(Temp[k])
    Temp_average.append(Ta)

#Now let's compute the map-to-map covariant matrix.

Number_pixels = 12*(256.**2)

C = np.zeros((6,6),dtype=float)

for i in range(6):
    for j in range(6):
        C[i,j] = (1./Number_pixels)*(np.sum((Temp[i]-Temp_average[i])*(Temp[j]-Temp_average[j])))


#Let's compute the weights.

C_inv = np.linalg.inv(C) #inverse covariant matrix

C_inv_sum = 0. #sum of all elements of C_inv
for i in range(6):
    for j in range(6):
        C_inv_sum += C_inv[i,j]

C_inv_lign = np.zeros(6,dtype=float) #sum over a column only
for i in range(6):
    for j in range(6):
        C_inv_lign[i]+=C_inv[i,j]
        
w = np.zeros(6,dtype=float) #weights
for i in range(6):
    w[i] = C_inv_lign[i]/C_inv_sum


#Let's multiply each map by its weight.

T_weight = []  #List of all maps with their coefficient

for i in range(6):  #affect coefficient
    tw = Temp[i]*w[i]
    T_weight.append(tw)
    
T = T_weight[0]  #sum maps
for i in range(1,6):
    T += T_weight[i]


hdu = pyfits.PrimaryHDU()
hdu.data = T
hdu.writeto('/home/adafne/Project/patch_final.fits')
"""

#To see the cluster:
"""
T = pyfits.getdata("/home/adafne/Project/patch_143.fits")
plt.imshow(T)
plt.show
"""


#Let's do it for other clusters.

################################################################

#Cluster number 20, SNR = 23.17 GLON = 6.498554 GLAT = 50.5633107
#Cluster number 12, SNR = 17.36 GLON = 3.931727 GLAT = -59.4110211
#Cluster number 182, SNR = 12.92534 GLON = 46.887288 GLAT = 56.4836198
#Cluster number 29, SNR = 10.053  GLON = 8.4744769 GLAT = -56.342089
#Cluster number 104, SNR = 9.00659 GLON = 29.803595 GLAT = -17.4010584
#Cluster number 15, SNR = 9.07 GLON = 4.4519 GLAT = -19.556
#Cluster number 99 SNR = 8.74 GLON = 28.7759 GLAT = -33.56406
#Cluster number 163 SNR = 8.27 GLON = 44.77205 GLAT = -51.30328
#Cluster number 8 SNR = 8.07 GLON = 2.8277 GLAT = 39.2388
#Cluster number 100 SNR = 7.52 GLON = 28.8925 GLAT = 60.13629
#Cluster number 41 SNR = 7.01 GLON = 12.59322 GLAT = -20.1076
#Cluster number 196 SNR = 6.10 GLON = 49.09337 GLAT = 25.239339
#Cluster number 102 SNR = 5.03 GLON = 29.55865 GLAT = -60.1678
#Cluster number 292 SNR = 4.57 GLON = 69.35866 GLAT = -15.58959

#Let's zoom on the patch around the cluster PSZ2 G006.76+30.45 (number 22 in the list).

"""
for k in range(6):
    m = hp.read_map("/home/adafne/Project/degrad"+freq[k]+".fits")
    y, x = np.indices((256,256)) #creates two indices for our patch
    w_patch = wcs.WCS(naxis=2) #creates a wcs object indicating the link btw our patch and the sky
    w_patch.wcs.crpix = [128,128]  #pixel in the center
    w_patch.wcs.cdelt = np.array([-0.01, 0.01])   #scale
    w_patch.wcs.crval = [69.35866,-15.58959]   #localisation in the sky
    w_patch.wcs.ctype = ["GLON-TAN", "GLAT-TAN"]  #type of coordinates
    theta,phi = w_patch.wcs_pix2world(y, x, 1)  #convert to angular coordinates
    theta = theta * math.pi /180.  #conversion degrees radians
    phi = phi * math.pi /180.      #conversion
    i = hp.ang2pix(512, math.pi/2.-phi, theta)   #find pixels of the maps corresponding to these coordinates
    patch = np.zeros((256,256)) #creates the patch
    patch[x,y] = m[i]  #fill with values of the map
    hdu = pyfits.PrimaryHDU(header=w_patch.to_header())
    hdu.data = patch
    hdu.writeto('/home/adafne/Project/Cluster292/patch_'+freq[k]+'.fits')
"""

#Visualise our patch
"""
for i in range(6):
    m = pyfits.getdata("/home/adafne/Project/Cluster182/patch_"+freq[i]+".fits")
    plt.figure()
    plt.imshow(m)
    plt.show()
"""


#We have to apply some conversion given in the table page 14 of http://planck.caltech.edu/pub/2013results/Planck_2013_results_09.pdf. The two last carts should be divided respectively by 58.04 and 2.27. The six should be multiplied by the coefficients given in the last column.
"""
conversion = np.array([-0.248, -0.359, 5.152, 0.161, 0.069/58.04, 0.038/2.27]) 

for i in range(6):
    data,header = pyfits.getdata('/home/adafne/Project/Cluster292/patch_'+freq[i]+'.fits', header = True)
    data = data*conversion[i]
    hdu = pyfits.PrimaryHDU()
    hdu.data = data
    hdu.writeto('/home/adafne/Project/Cluster292/patch_convert_'+freq[i]+'.fits')
"""




# We want now to make the linear combination. Let's write the covariant matrix.

#First let's get the temperature.
"""
Temp = []

for k in range(6):
    T = pyfits.getdata('/home/adafne/Project/Cluster292/patch_convert_'+freq[k]+'.fits')
    Temp.append(T)

#We need the average for each map.

Temp_average = []

for k in range(6):
    Ta = np.mean(Temp[k])
    Temp_average.append(Ta)

#Now let's compute the map-to-map covariant matrix.

Number_pixels = 12*(256.**2)

C = np.zeros((6,6),dtype=float)

for i in range(6):
    for j in range(6):
        C[i,j] = (1./Number_pixels)*(np.sum((Temp[i]-Temp_average[i])*(Temp[j]-Temp_average[j])))


#Let's compute the weights.

C_inv = np.linalg.inv(C) #inverse covariant matrix

C_inv_sum = 0. #sum of all elements of C_inv
for i in range(6):
    for j in range(6):
        C_inv_sum += C_inv[i,j]

C_inv_lign = np.zeros(6,dtype=float) #sum over a column only
for i in range(6):
    for j in range(6):
        C_inv_lign[i]+=C_inv[i,j]
        
w = np.zeros(6,dtype=float) #weights
for i in range(6):
    w[i] = C_inv_lign[i]/C_inv_sum


#Let's multiply each map by its weight.

T_weight = []  #List of all maps with their coefficient

for i in range(6):  #affect coefficient
    tw = Temp[i]*w[i]
    T_weight.append(tw)
    
T = T_weight[0]  #sum maps
for i in range(1,6):
    T += T_weight[i]


hdu = pyfits.PrimaryHDU()
hdu.data = T
hdu.writeto('/home/adafne/Project/Cluster292/patch_final.fits')

"""

#To see the cluster:
"""
T = pyfits.getdata("/home/adafne/ILC_SZ_comp_2048.fits")
plt.imshow(T)
plt.show()
"""
#To compare with whole map :
"""
m = hp.read_map('/home/adafne/ILC_SZ_comp_2048.fits')
hp.mollzoom(m, title="353",norm='hist')
plt.show()
"""



