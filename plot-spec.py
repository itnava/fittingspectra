
# script for displaying spectra

import pyfits
import matplotlib.pyplot as plt


#function for redshifting wavebands
def redden(wave, vel):
    return wave * (1 + vel / (3 * 10 ** 5))
 
count = 0

#list of files to be plotted, saved in a file.
filename = [line.strip() for line in open("files_plot",'r')]

#List of emission line wavelengths and names for labelling
emission = [4861, 5007, 6300, 6548, 6563, 6584, 6724]
label = ['H$\beta', 'OIII', 'OI', 'NII', 'H\alpha', 'NII', 'SII']

#Switching on interactive plotting
plt.ion()

for path in filename:
    print path
    #script was crashing for some sources with non-standard names(probably have non-standard processing).
    if len(path) <= 45:
        with pyfits.open(path) as hdu:
            count += 1
            print count
            hdu[0].header['OBJECT']
            hdu[0].header['RFN']
            vel = hdu[0].header['VELOCITY']            
            ref_pixel = hdu[0].header['CRPIX1']
            coord_ref_pixel = hdu[0].header['CRVAL1']
            wave_pixel = hdu[0].header['CD1_1']
            wave0 = (coord_ref_pixel -((ref_pixel - 1) * wave_pixel))#/(1 + vel / (3 * 10 ** 5))
            counts = hdu[0].header['naxis1']
            waves_array = [wave0 + i * wave_pixel for i in range(counts)]
            sig0 = hdu[0].data[0][0]
            top_of_plot = max(sig0) - 0.05 * max(sig0)
            plt.figure(111)
            plt.xlabel('wavelength')
            plt.ylabel('counts')
            i = 0
            #annotation
            for val in emission:
                plt.axvline(x = redden(val, vel), color = 'r')
#                plt.annotate(label[i], xy = (redden(val, vel) + 2, top_of_plot + (-1)**i * 20))
                i += 1
#            plt.axvline(x = redden(5007, vel), color = 'r')
#            plt.axvline(x = redden(4861, vel), color = 'r')
#            plt.axvline(x = redden(6584, vel), color = 'r')
#            plt.axvline(x = redden(6548, vel), color = 'r')
#            plt.axvline(x = redden(6563, vel), color = 'r')
#            plt.axvline(x = redden(6300, vel), color = 'r')
#            plt.axvline(x = redden(6724, vel), color = 'r')
#            plt.annotate('OIII', xy = (redden(5007, vel) + 2, top_of_plot))
            plt.plot(waves_array, sig0)
            plt.draw()
            #introducing a button press for interactive plotting
            plt.waitforbuttonpress(timeout = -1)
            plt.clf()
            plt.close()
            
