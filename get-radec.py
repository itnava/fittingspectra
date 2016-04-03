
import pyfits

#FAST_files: 
fast_input = open("files",'r')
fast_id = open("Fast_id",'w')


for line in fast_input:
    path = line.strip()
    with pyfits.open(path) as hdu:
        name = hdu[0].header['OBJECT']
        ra = hdu[0].header['RA']
        ra_ha =ra.split(':')
        ra_deg = float(ra_ha[0]) *360 /24 + float(ra_ha[1]) * 360 / 24 /60 + float(ra_ha[2]) /3600 * 360/ 24
        dec = hdu[0].header['DEC']
        dec_ha = dec.split(':')
        dec_deg = float(dec_ha[0]) + float(dec_ha[1]) / 60 + float(dec_ha[2]) /3600
        velocity = hdu[0].header['VELOCITY']
        fast_id.write("%s %s %s %s %s\n" % (path, name, ra_deg, dec_deg, velocity))
        
fast_input.close()