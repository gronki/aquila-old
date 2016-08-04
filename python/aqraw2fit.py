#!/bin/python

from aqutils import *
from os import read,write

for raw_fn in read(0,10000).split('\n'):
    if raw_fn == '':
        continue

    fit_fn = path.splitext(raw_fn)[0]+".fit"

    (imdata,h) = aql_raw2fits(raw_fn,fit_fn)

    h['CCDTYPE'] = ''

    # remove file in case of overwriting
    if path.exists(fit_fn):
        aql_print( u'File exists. Deleting...' )
        remove(fit_fn)

    write(1,path.relpath(fit_fn) + '\n')

    # save FITS
    pyfits.writeto(fit_fn,imdata,header=h)

    aql_print(u'result written to '+fit_fn)
    # write(1,fit_fn+'\n')

    aql_print(u' --------------------- ')
    aql_print(u'Conversion finished. ')
