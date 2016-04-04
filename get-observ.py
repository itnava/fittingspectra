#Converting coordinates to decimal from HA:Min:Sec and Deg: Min:Sec to help with plotting as well as for removing duplicates and observations of different areas on the same galaxy
import pyfits

def trunc(item,n):
    if type(item) is float:
        slen = len('%.*f' %(n, item)) #rounds coordinate to n decimals, measures length
        return str(item)[:slen] # truncates coordinate so that n decimals are included
    elif type(item) is str:
        return str(item)[:n] # truncates coord to length n

#FAST_files: 
fast_input = open("files",'r')
fast_unique = open("sources_unique.dat",'w')
fast_repeat = open("sources_repeat.dat",'w')


FLWO = 0
CTIO = 0

radec_dict = dict()
ra_dict = dict()
dec_dict = dict()

radec =[]
count = 0
z = []
repeat = 0
name_list = []
ra_list = []
dec_list = []

for line in fast_input:
    path = line.strip()
    with pyfits.open(path) as hdu:
        name = hdu[0].header['OBJECT']
        observatory = hdu[0].header['OBSERVAT']
        ra = hdu[0].header['RA']
        ra_ha =ra.split(':')
        ra_deg = float(ra_ha[0]) *360 /24 + float(ra_ha[1]) * 360 / 24 /60 + float(ra_ha[2]) /3600 * 360/ 24
        dec = hdu[0].header['DEC']
        dec_ha = dec.split(':')
        dec_deg = float(dec_ha[0]) + float(dec_ha[1]) / 60 + float(dec_ha[2]) /3600
        velocity = hdu[0].header['VELOCITY']
        if observatory == 'flwo1':
            FLWO += 1
        else:
            print name
            CTIO += 1
        radec_item = ("%s, %s" %(trunc(ra_deg,3), trunc(dec_deg,3)))
        name_trunc = trunc(name,17)
        ra_trunc = trunc(ra_deg,2)
        dec_trunc = trunc(dec_deg, 2)
#        if name_trunc not in name_list:
#            name_list.append(name_trunc)
#            fast_unique.write("%s %s %s %s %s\n" %(name, velocity, trunc(ra_deg,3), trunc(dec_deg,3), path))
#        else:
#            fast_repeat.write("%s %s %s %s %s\n" %(name, velocity, trunc(ra_deg, 3), trunc(dec_deg,3), path))
        if ra_trunc not in ra_list:
            ra_list.append(ra_trunc)
            dec_list.append(dec_trunc)
            z.append(velocity)
            fast_unique.write("%s %s %s %s %s\n" %(name, velocity, ra_trunc, dec_trunc, path))
        elif dec_trunc != dec_list[ra_list.index(ra_trunc)]:
            ra_list.append(ra_trunc)
            dec_list.append(dec_trunc)
            z.append(velocity)
            fast_unique.write("%s %s %s %s %s\n" %(name, velocity, ra_trunc, dec_trunc, path))
        else:
            fast_repeat.write("%s %s %s %s %s\n" %(name, velocity, ra_trunc, dec_trunc, path))
        

print count,repeat
print radec
print FLWO, CTIO
fast_input.close()
fast_unique.close()
fast_repeat.close()
