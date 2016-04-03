
import numpy as np
import scipy.optimize
import lmfit
from lmfit import minimize, Parameters, report_fit, Parameter
from scipy.optimize import curve_fit
from scipy.integrate import simps
from matplotlib.pyplot import savefig
import matplotlib.pyplot as plt
import pyfits
from matplotlib import pyplot
from matplotlib.pyplot import subplot

# I have to first write a script to test on SDSS data, then on FAST.
#FAST data have a different structure and dimension, so the steps for extracting spectrum 
# will have to be redone.
#I cannot figure out how to read a FAST spectrum. Start with SDSS instead?
# include chi square fit. chisquare/dof 1+/- 0.05.
# just four gaussians in NII-Halpha. 
# preliminary fit to identify lines, mask them
# fit continuum and subtract.
# then do a fit for usable lines to get flux.
# fit to actual lines only. limit width
# constrain the width to be less than fwhm >< 20km/s, 12 A at the base
# final fit only on amplitude for narrow lines.
# let broad line halpha have narrow range in centroid and let width and amplitude vary.
# try a broad line sy1. to test the fit.



#defining a single gaussian and a multigaussian model
def func(x, a, b, c, d):
    return a*np.exp(-(x-b)**2/(2.*c**2)) + d

def resid_sing(params, x, data):
    a = params['a'].value
    b = params['b'].value
    c = params['c'].value
    d = params['d'].value
    
    model = a*np.exp(-(x-b)**2/(2.*c**2)) + d
    
    return data - model

def sing_gauss(x, a, b, c):
    return a*np.exp(-(x-b)**2/(2.*c**2))

def multi_gauss(x, a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, f) :
    return sing_gauss(x, a1, a2, a3) + sing_gauss(x, b1, b2, b3) + sing_gauss(x, c1, c2, c3) + sing_gauss(x, d1, d2, d3) + f

def resid_multi(params, x, data):
    a1 = params['a1'].value
    a2 = params['a2'].value
    a3 = params['a3'].value
    b1 = params['b1'].value
    b2 = params['b2'].value
    b3 = params['b3'].value
    c1 = params['c1'].value
    c2 = params['c2'].value
    c3 = params['c3'].value
    d1 = params['d1'].value
    d2 = params['d2'].value
    d3 = params['d3'].value
    f = params['f'].value
    
    model = sing_gauss(x, a1, a2, a3) + sing_gauss(x, b1, b2, b3) + sing_gauss(x, c1, c2, c3) + sing_gauss(x, d1, d2, d3) + f    

    return data - model

def errfunc_gauss(p, x, y):
    return y - func(p, x)

def errfunc_multi(p, x, y):
    return y - multi_gauss(p, x)
    
#reading the bandpass file. 
indexdef = open("jph_gauss_detect.bands",'r')

#reading the fits file. to be modified for a batch job.


filename = [line.strip() for line in open("files_test",'r')]


plt.ion()

for path in filename:
    print path
    if len(path) <= 60:
        with pyfits.open(path) as hdu:
            name = hdu[0].header['OBJECT']
            rfn = hdu[0].header['RFN']
            vel = hdu[0].header['VELOCITY']            
            ref_pixel = hdu[0].header['CRPIX1']
            coord_ref_pixel = hdu[0].header['CRVAL1']
            wave_pixel = hdu[0].header['CD1_1']
            wave0 = (coord_ref_pixel -((ref_pixel - 1) * wave_pixel))
            counts = hdu[0].header['naxis1']
            waves_array = [wave0 + i * wave_pixel for i in range(counts)] 
            #for subtracted data, data[0], for FAST data data[0][0]
            sig0 = hdu[0].data[0]
            z = vel  / (3 * 10**5)
            n248 = 6548 * (1+z)
            ha = 6562 * (1+z)
            n284 = 6584 * (1+z)
            #resetting detection counters
            detect = []
            #redshifting wavelengths of interest
            emission = [5007 * (1+z), 6300 * (1+z), 6724 * (1+z)]
            j = 0
            values = []
            indexdef.seek(0)
            for line in indexdef:
                bands_line = line.split()
                bands = list(bands_line)
                i = 1 
                for i in range(1,7):
                    bands[i] = float(bands_line[i]) * (1 + z) 
                    i += 1
                bands_1 =  float(bands[1]) - 100
                bands_6 =  float(bands[6]) + 100
                print bands_1, bands_6
                print min(waves_array), max(waves_array)
                waves = []
                flux = []
                y_fit = []
                y_residual = []
                sig_wave = []
                sig_count = []
                noise_wave = []
                noise_count = []
                #first guess for fit. func defined above, in the order (x, a,b,c,d)
                #a: amplitude, b: peak wavelength, c: band width, d: a constant
                i = 0
                for col in waves_array:
                #    #only interested in data between the bands specified in the .bands file
                    if bands_1 <= float(col) <= bands_6:
                        waves.append(col)
                        flux.append(sig0[waves_array.index(col)])
                    i += 1
                print len(waves), len(flux)
                fig = subplot(111)
                fig.plot(waves_array, sig0, 'r')
                fig.plot(waves, flux, 'b')
                fig.axis(xmin = min(waves_array), xmax = max(waves_array))
                x = np.array(waves)
                y = np.array(flux)
                do_flag = True
#                if bands[0] == "NII_Ha":
#                    p_guess=[10, n248, 2, 25, ha, 2, 25, ha, 8, 10, n284, 2, 25]
#                else:        
#                    p_guess=[10, emission[j], 6, 25]
#                    j += 1
#                do_flag = True   
                if bands[0] == "NII_Ha":
                    params = Parameters()
                    params.add('a1', value = 10., min = 0)
                    params.add('a2', value = ha, vary = False)
                    params.add('a3', value = 2., min = 0)
                    params.add('b1', value = 10., min = 0)
                    params.add('b2', value = ha, vary = False)
                    params.add('b3', value = 12., min = 0, max = 12)
                    params.add('c1', value = 10., min = 0)
                    params.add('c2', value = n248, vary = False)
                    params.add('c3', value = 2., min = 0, max = 12)
                    params.add('d1', value = 10., min = 0)
                    params.add('d2', value = n284, vary = False)
                    params.add('d3', value = 2., min = 0, max = 12)
                    params.add('f', value = 0., min = 0)
                    residual = minimize(resid_multi, params, args=(x, y))  
                    print report_fit(params)
                    a1 = params['a1'].value
                    a2 = params['a2'].value
                    a3 = params['a3'].value
                    b1 = params['b1'].value
                    b2 = params['b2'].value
                    b3 = params['b3'].value
                    c1 = params['c1'].value
                    c2 = params['c2'].value
                    c3 = params['c3'].value
                    d1 = params['d1'].value
                    d2 = params['d2'].value
                    d3 = params['d3'].value
                    f = params['f'].value
                    do_flag = True
                    values = [a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, f]
                    noise_wave = []
                    noise_count = []
                    sig_wave = []
                    sig_count = []
                    i = 0
                    y_fit = resid_multi(params, x, y) + y
                    """
                        for i in range(0, len(y)):
                            y_fit.append(multi_gauss(x[i], *values))
                            y_residual.append(y[i] - y_fit[i])
                            if bands[3] <= x[i] <= bands[4]:
                                sig_wave.append(x[i])
                                sig_count.append(y[i])
                            else :
                                noise_wave.append(x[i])
                                noise_count.append(y[i])
                            i += 1
                        print type(noise_count)
                        print type(sig_count)
                    """
                    fit= fig.plot(x, y_fit)
                    plt.setp(fit, color = 'black')
                    fig.plot(x, y_residual) 
                    fig.axis(xmin = min(waves_array), xmax = max(waves_array))
                    #detect.append(1)
                    print j
                    #                        do_flag = False
                    #   detect.append(0)
                    #  print j
                else:
                    print emission[j] , j
                    params = Parameters()
                    params.add('a', value = 10., min = 0.)
                    params.add('b', value = emission[j], vary = False)
                    params.add('c', value = 2., min = 0., max= 12.)
                    params.add('d', value = 0., min = 0.)
                    #try: 
                    residual = minimize(resid_sing, params, args=(x, y))    
                    print report_fit(params)
                    do_flag = True
                    a = params['a'].value
                    b = params['b'].value
                    c = params['c'].value
                    d = params['d'].value
                    values.append(a)
                    values.append(b)
                    values.append(c)
                    values.append(d)
                    noise_wave = []
                    noise_count = []
                    sig_wave = []
                    sig_count = []
                    i = 0
                    y_fit = resid_sing(params, x, y) + y
                    """
                        for i in range(0, len(y)):
                            y_fit.append(func(x[i], a, b, c, d))
                            y_residual.append(y[i] - y_fit[i])
                            if bands[3] <= x[i] <= bands[4]:
                                sig_wave.append(x[i])
                                sig_count.append(y)
                            else :
                                noise_wave.append(x[i])
                                noise_count.append(y)
                            i += 1
#                        print noise_count
#                        print sig_count
                    """
                    fit= fig.plot(x, y_fit)
                    plt.setp(fit, color = 'black')
                    fig.plot(x, y_residual) 
                    detect.append(1)
                    print j
                    #except:
                    #     do_flag = False
                    while(j<2): 
                         j+=1
                         print j
        plt.draw()
        plt.waitforbuttonpress(timeout = -1)
        plt.clf()
        plt.close()
        print detect
    else:
        print "Stop"
indexdef.close()

