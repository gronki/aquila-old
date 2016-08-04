#!/bin/python
from aqutils import *
from os import read,write,path
from numpy import mean, median, std
import pyfits

for fn in read(0,10000).split('\n'):
    if fn == '':
        continue
    im = pyfits.getdata(fn)
    hdr = pyfits.getheader(fn)
    aql_print( 'read: ' + fn )
    hdr['HISTORY'] = ' EQUALIZED WITH GREYFLAT UTILITY'
    hdr['HISTORY'] = ' ...on ' + datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    (im2,hdr) = aql_cfaflat(im,hdr)
    hdr['CCDMEAN'] = mean(im2)
    hdr['CCDMED'] = median(im2)
    hdr['CCDSTD'] = std(im2)
    fn_split = path.splitext(fn)
    fn2 = fn_split[0] + "+eql" + fn_split[1]
    if path.exists(fn2):
        aql_print( u'File exists. Deleting...' )
        remove(fn2)
    pyfits.writeto(fn2,im2,header=hdr)
    aql_print( 'saved as: ' + fn2 )
    write(1,path.relpath(fn2)+'\n')
    aql_print( ' ==================================' )
