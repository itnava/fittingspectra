import scipy.optimize
import astropy.io.fits as pyfits
#from matplotlib.pyplot import savefig
import matplotlib.pyplot as plt
from matplotlib.pyplot import subplot

#from redshiftbandpass import redshift
from single_gauss import singlegauss
from multi_gauss import multigauss



"""
This program takes individual spectra, and tries to fit gussian emission lines at four different wavelengths.
We need to use a multigaussian profile to fit the NII + Halpha emission near 6560A and sig=ngle gaussian profiles for OIII 5007A, SII 6724A and OI 6300A.

The bandpass for the filters used for observation are defined in jph_gauss_detect.bands file. We expect that 99% of any line emission from the host galaxy will be observed within the bandpass defined here.
"""
indexdef = open("jph_gauss_detect.bands",'r')

"""
The paths for spectrum file of interest are added to input_list.

"""
filename = [line.strip() for line in open("sdss_file",'r')]

#Output file for line emission. We will use this to obtain ratios and classify the galaxies.
linewidth = open("linewidths.txt", 'w')


for path in filename:
    #Printing path to check for missing files.
    print(path)
    #Opening fits file
    with pyfits.open(path) as hdu:

        data = hdu[0].data
        counts = [float(row) for row in data[0]]
        sigma = [float(row) for row in data[2]]

        # reading coeff0 and coeff1 to generate the wavelength table
        coeff0 = hdu[0].header["COEFF0"]
        coeff1 = hdu[0].header["COEFF1"]

        # generating wavelength table for vaccuum wavelengths.
        vac = [10.0 ** (coeff0 + coeff1 * i) for i in range(0, len(counts))]
        # converting them to air wavelengths
        air = [row / (1.0 + 2.735182 * 10 ** (-4) + 131.4182 / row ** 2 + (2.76249 * 10 ** 8) / row ** 4) for row in vac]

        # redshift
        z = hdu[0].header["Z"]

        #for compatibility with fast data, keeping semilar variable names
        wavelength_array = air
        signal = counts

        #extracting details from header, we will use these later while plotting
        name = hdu[0].header['NAME']
#        rfn = hdu[0].header['RFN']


        #Converting rest wavelengths for emission lines into observed wavelengths.
        n2ha = [6548 * (1+z), 6562 * (1+z), 6584 * (1+z)]
        singleemissionlines = [4861 * (1 + z), 5007 * (1 + z), 6300 * (1+z), 6724 * (1+z)]
        j = 0

        #Reinitializing bandpass file.
        indexdef.seek(0)
        bands = []

        linewidth.write("%s " % (name))

        for line in indexdef:

            #We do the fitting for each bandpass separately.
            bands_line = line.split()

            # Initializing the arrays needed for modeling.
            waves = [] #Wavelength range of interest
            flux = []  #Flux recorded for wavelength range of interest
            y_fit = []  #Results of the fit
            y_residual = [] #Residuals obtained by subtracting signal from fit

            #Initializing the arrays for measuring signal and background
            signal_wave = []
            signal_count = []
            background_wave = []
            background_count = []
            x_fit = []

            #redshifting the bandwidth for the bandpass
            lowband = float(bands_line[1]) * (1+z) - 100
            highband = float(bands_line[6]) * (1 + z) + 100
            signalrange = [float(bands_line[3]) * (1+z), float(bands_line[4]) * (1+z)]

            #Selecting relevant data for modeling
            for value in wavelength_array:
                if lowband <= value <= highband:
                    waves.append(value)
                    flux.append(signal[wavelength_array.index(value)])

            #This flag indicates whether or not a fit was obtained.
            do_flag = False

            #Since the NII doublet and Halpha emission bandpasses overlap,they have to be fitted simultaneously using four gaussians. We know the central wavelength and the approximate widths which we will use for the initial guess since the curve_fit module requires a robust initial guess for convergence.

            if bands_line[0] == "NII_Ha":
                fitfunction = multigauss
                initial_guess = [10, n2ha[0], 2, 25, n2ha[1], 2, 25, n2ha[1], 8, 10, n2ha[2], 2, 25]
                try:
                    optimization, covariance = scipy.optimize.curve_fit(multigauss, waves, flux, p0 = initial_guess)
                    do_flag = True
                except:
                    optimization, covariance = initial_guess, "NONE"
                    do_flag = False

            # The remaining lines can be modeled as single gaussians, while they do have structure, the instrumentation did not have high enough resolution.

            else:
                initial_guess = [10, singleemissionlines[j], 6, 25 ]
                fitfunction = singlegauss
                print(singleemissionlines[j], j)
                j += 1

                try:
                    optimization, covariance = scipy.optimize.curve_fit(singlegauss, waves, flux, p0 = initial_guess)
                    print(optimization,initial_guess)
                    do_flag = True
                except:
                    optimization, covariance = initial_guess, "NONE"
                    do_flag = False

            #The try-except statements, allow us to fit the data where convergrence was possible, and assign a dummy value to the non-converging spectra. These can be manually examined later.

            # If we were able to fit the data, we want to measure the signal to noise ratio. Otherwise we set it to 0

            if optimization[0] <= 0 or covariance[0][0] == "NONE":
                do_flag = False
                linesignal = 0

            else:
                k = 0
                for k in range(0, len(flux)):
                    y = fitfunction(waves[k], *optimization)
                    y_fit.append(y ) #fitfunction, waves[k], *optimization)
                    y_residual.append(flux[k] - y_fit[k])
                    x_fit.append(waves[k])

                    if signalrange[0] <=  waves[k] <= signalrange[1]:
                        signal_wave.append(waves[k])
                        signal_count.append(y_fit[k])
                    else:
                        background_wave.append(waves[k])
                        background_count.append(y_fit[k])
                k += 1

                try:
                    lineemission = sum(signal_count) / (max(signal_wave) - min(signal_wave))
                    backgroundemission = sum(background_count) / ((max(background_count) - min(background_count)) - (max(signal_wave) - min(signal_wave)))

                    linesignal = lineemission / backgroundemission
                    print(linesignal)
                except:
                    linesignal = 0
                    print(0)

            #Writing the name of the line and the signal to noise ratio to output file

            linewidth.write("%s %f " % (bands_line[0], linesignal))
            #Plotting the data, fit and residual

            figure1 = subplot(311)
            figure1.plot(wavelength_array, signal, 'r')
            figure1.axis(xmin = min(wavelength_array), xmax = max(wavelength_array))
            figure2= subplot(312)
            figure2.plot(x_fit, y_residual, 'g')
            figure2.axis(xmin = min(wavelength_array), xmax = max(wavelength_array))
            figure3 = subplot(313)
            figure3.plot(x_fit, y_fit, 'b')
            figure3.axis(xmin = min(wavelength_array), xmax = max(wavelength_array))

        linewidth.write("\n")
        plt.draw()
        plt.waitforbuttonpress(timeout = -1)
        plt.clf()
        plt.close()

#Closing files.
linewidth.close()
indexdef.close()