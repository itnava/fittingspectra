
import numpy as np
import scipy.optimize
from scipy.optimize import curve_fit
from matplotlib.pyplot import savefig
import matplotlib
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
    
def multi_gauss(x, a1, a2, a3, b1, b2, b3, c1, c2, c3, d1, d2, d3, e1, e2, e3, f) :
    return func(x, a1, a2, a3, 0) + func(x, b1, b2, b3, 0) + func(x, c1, c2, c3, 0) + func(x, d1, d2, d3, 0) + func(x, e1, e2, e3, 0)+ f
 
#reading the bandpass file. 
indexdef = open("jph_mod_test.bands",'r')

#reading the fits file. to be modified for a batch job.
fitsfile = pyfits.open("spSpec-51984-0279-039.fit") #this one seems to be working

#fitsfile = pyfits.open("spSpec-51699-0349-514.fit") #works

#fitsfile = pyfits.open("spSpec-51909-0462-261.fit") #sort of

#fitsfile = pyfits.open("spSpec-51985-0295-071.fit")

#for SDSS, flux counts are in [0], continuum subtracted data in [1] and sigma in [2]
data = fitsfile[0].data
counts = [float(row) for row in data[0]]
sigma = [float(row) for row in data[2]]

#reading coeff0 and coeff1 to generate the wavelength table
coeff0 = fitsfile[0].header["COEFF0"]
coeff1 = fitsfile[0].header["COEFF1"]

# generating wavelength table. 
vac = [10.0**(coeff0 + coeff1 * i) for i in range(0,len(counts))]
# these are vacuum wavelengths. Have to convert them to air
air = [row/(1.0 + 2.735182 * 10**(-4) + 131.4182 / row**2 + (2.76249 * 10**8) / row**4) for row in vac]

#redshift
z = fitsfile[0].header["Z"]
#redshifting wavelengths of interest
hb = 4861 * (1+z)
o3 = 5007 * (1+z)
n284 = 6584 * (1+z)
ha = 6562 * (1+z)
n248 = 6548 * (1+z)

print hb, o3 , n248, ha, n284

# fitting Hbeta, skip first line
indexdef.readline()
line = indexdef.readline()
bands_hb = line.split()
bandshb_3 = float(bands_hb[3]) * (1 + z) - 100
bandshb_4 = float(bands_hb[4]) * (1 + z) + 100


# fitting OIII
line = indexdef.readline()
bands_o = line.split()
#redshifting and establishing bands for gaussian fit.
bandso_3 = float(bands_o[3]) * (1 + z) - 100
bandso_4 = float(bands_o[4]) * (1 + z) + 50

o34959 = 4959.0 *  (1 + z)
print o34959
print hb
print o3
print n248
print ha
print n284

#Creating lists for storing wavelength and counts for bands of interest.
waves_o = []
flux_o = []

#first guess for fit. func defined above, in the order (x, a,b,c,d)
#a: amplitude, b: peak wavelength, c: band width, d: a constant

po_guess=[10, o3, 6, 30]

i = 0

for col in air:
    #only interested in data between the bands based on the .bands file
    if bandso_3  <= float(col) <= bandso_4:
        waves_o.append(col)
        flux_o.append(counts[i])
    i += 1

#Need to plot data in the two areas of interest. 
fig = subplot(211)
fig.plot(waves_o, flux_o, linestyle = '--')

#defining x and y   
y = flux_o
x = waves_o
#Defining original data points for calculating final degrees of freedom and limits on how many fits.
fit_points = len(waves_o)
lim = 0.9 * fit_points 

fit_param = []


do_flag = True
#try-except loop to avoid fits with no convergence
try:
    popt_o, pcov_o = scipy.optimize.curve_fit(func, x, y, p0 = po_guess)
except:
    popt_o, pcov_o = po_guess, None
    do_flag = False
    print "no fit"
#print popt_o

if popt_o[0] <= 0:
    do_flag = False
else:
    fit_param.append(popt_o[1])
    fit_param.append(popt_o[2])

print fit_param

while do_flag == True:
#Calculating residuals
    for i in enumerate(y):
        yo_fit = func(x, *popt_o)
        yo_residual = y - yo_fit
    fig.plot(waves_o, yo_fit, linestyle="--")
    fig.plot(waves_o, yo_residual)
    if max(yo_residual) >= (3) * np.mean(abs(yo_residual)):
        for i in range(0, len(yo_residual)):
             if yo_residual[i] == max(yo_residual):
                 x_max = x[i]
                 i = len(yo_residual)
             else:
                 i += 1             
        po_guess = [max(yo_residual), x_max,10, popt_o[3]] 
        y = yo_residual + yo_fit[3]
        try:
            popt_o, pcov_o = scipy.optimize.curve_fit(func, x, y, p0 = po_guess)
            dof = fit_points - len(popt_o)
            print dof
            fit_points = dof
            if dof > lim:
#            do_flag = True
                fit_param.append(popt_o[1])
                fit_param.append(popt_o[2])
                print fit_param
            else:
                do_flag = False
                print "no fit"
        except:
            popt_o, pcov_o = po_guess, None
            do_flag = False
            print "no fit"
    else:
        do_flag = False
print fit_param
n = len(fit_param)
print n
mask_o = []
q = []
for i in range(0, n-1):
    for j in range(0, len(waves_o)):
         if fit_param[i] - fit_param[i + 1] <= waves_o[j] <= fit_param[i]+ fit_param[i+1]:
             mask_o.append(waves_o[j])
             q.append(10)
             j += 1
    i += 1

print mask_o
fig.plot(mask_o, q)

line = indexdef.readline()
bands_ha = line.split()
print bands_ha[3]
print bands_ha[4]

bandsha_3 = float(bands_ha[3]) * (1 + z) - 50
bandsha_4 = float(bands_ha[4]) * (1 + z) + 50

line = indexdef.readline()
bands_n248 = line.split()
print bands_n248[3]
print bands_n248[4]

bandsn248_3 = float(bands_n248[3]) * (1 + z) - 50
bandsn248_4 = float(bands_n248[4]) * (1 + z) + 50

line = indexdef.readline()
bands_n284 = line.split()
print bands_n284[3]
print bands_n284[4]

bandsn284_3 = float(bands_n284[3]) * (1 + z) - 50
bandsn284_4 = float(bands_n284[4]) * (1 + z) + 50

waves_nh = []
flux_nh = []

#first guess for fit. func defined above, in the order (x, a,b,c,d)
#a: amplitude, b: peak wavelength, c: band width, d: a constant

pnh_guess=[10, n248, 6, 10, n248, 6, 25, ha, 6, 25, ha, 8, 10, n284, 2, 35]

i = 0

for col in air:
    #only interested in data between the bands specified in the .bands file
    if bandsn248_3  <= float(col) <= bandsn284_4:
        waves_nh.append(col)
        flux_nh.append(counts[i])
    i += 1
fig = subplot(2,1,2)
fig.plot(waves_nh, flux_nh, linestyle = '--')
     
y = flux_nh
x = waves_nh

fit_points = len(waves_nh)
lim = 0.9 * fit_points

try:
    popt_nh, pcov_nh = scipy.optimize.curve_fit(multi_gauss, x, y, p0 = pnh_guess)
except:
    popt_nh, pcov_nh = pnh_guess, None
print popt_nh
#print pcov_nh
for i in enumerate(y):
    ynh_fit = multi_gauss(x, *popt_nh)
    ynh_residual = y - ynh_fit
fig.plot(waves_nh, ynh_fit, linestyle="--")
fig.plot(waves_nh, ynh_residual)
print max(ynh_residual)
print np.mean(abs(ynh_residual))
fit_points = len(waves_nh) - len(popt_nh)
print fit_points

do_flag = True
j = 0
for i in range(0, len(ynh_residual)):
        if ynh_residual[i] == max(ynh_residual):
             x_max = x[i]
             i = len(ynh_residual)
        else:
             i += 1             

pnh_guess = [max(ynh_residual), x_max,10, popt_nh[3]] 
y = ynh_residual + ynh_fit[3]

while do_flag == True:
    try:
        popt_nh, pcov_nh = scipy.optimize.curve_fit(func, x, y, p0 = pnh_guess)
    except:
        popt_nh, pcov_nh = pnh_guess, None
        do_flag = False
    print popt_nh
#    print pcov_nh
    for i in enumerate(y):
        ynh_fit = func(x, *popt_nh)
        ynh_residual = y - ynh_fit
    fig.plot(waves_nh, ynh_fit, linestyle="--")
    fig.plot(waves_nh, ynh_residual)
    dof = fit_points - len(popt_nh)
    print dof
    fit_points = dof
    if max(ynh_residual) >= (3) * np.mean(abs(ynh_residual)):
        for i in range(0, len(ynh_residual)):
             if ynh_residual[i] == max(ynh_residual):
                 x_max = x[i]
                 i = len(ynh_residual)
             else:
                 i += 1             
        pnh_guess = [max(ynh_residual), x_max,10, popt_nh[3]] 
        print pnh_guess
        y = ynh_residual + ynh_fit[3]
        if dof > lim:
            do_flag = True
            print "no fit"
        else:
            do_flag = False
    else:
        do_flag = False



#try:
#    popt_o, pcov_o = scipy.optimize.curve_fit(func, x, y, p0 = po_guess)
#except:
#    popt_o, pcov_o = p_guess, None

#print popt_o
#print pcov_o

#y = yo_residual

#for i in enumerate(y):
#    yo_fit = func(x, *popt_o)
#    yo_residual = y - yo_fit 

#if np.fabs(max(yo_residual)) >= 1.5* np.fabs(np.mean(yo_residual)):
#    print True

#print max(yo_residual)
#print np.mean(yo_residual)


"""

#first guess for fit. func defined above, in the order (x, a,b,c,d)
#a: amplitude, b: peak wavelength, c: band width, d: a constant

   
   
# recording initial function for plotting later. linspace returns evenly spaced samples
x_func = np.linspace(min(waves_nh), max(waves_nh), len(waves_nh))
initial_plot = multi_gauss(x_func, *pnh_guess)

# fitting function to data. Assuming 20% error in y for now. May need to come up with a more
#realistic way to estimate error, later. 
#Using try-except so that the script does not crash because it could not find a line.
#If line does not exist, the initial guess is returned with no covariance matrix
#maxfev is the number of iterations

    # Calculating residuals ( difference between data and fit). 
    #y_fit calculated by using the function described by the fitting parameters and x
sigma = [row for row in data[2]]
"""


savefig("test.png")


fitsfile.close()
indexdef.close()

