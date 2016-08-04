#!/usr/bin/python

from numpy import array,mean,median,zeros,std
from os import getcwd, path, remove, read, write
import pyfits
from PIL import Image
import subprocess
import io
from datetime  import datetime
import json

"""
dcraw -c -o 0 -T -4 -t 0 -r 2 1 1.8 1 -q 3 light1.NEF > a.tif
dcraw -c -o 0 -T -4 -t 0 -d light1.NEF > b.tif
"""

# datetime.datetime.now().strftime("%Y/%m/%d %H:%M:%S")

def aql_print(s):
    # apply date-time information and print
    # %Y-%m-%d
    ss = ' [ ' + datetime.now().strftime("%H:%M:%S") + ' ] ' + s
    write(2, ss + '\n')

def aql_splitfn(fn):
    (fndir,fn0) = path.split(fn)
    (fname,fnext) = path.splitext(fn)
    a = fname.split('+')
    fnbase = a[0]
    fnsufx = ''
    if len(a) > 1:
        fnsufx = a[-1]
    return (fndir, fnbase, fnsufx, fnext)

def aql_joinfn(fndir, fnbase, fnsufx, fnext):
    fn0 = fnbase
    if fnsufx != '':
        fn0 = fn0 + '+' + fnsufx
    return path.join( fndir, fn0 + fnext )


def aql_loadraw(raw_fn):
    "Blabla"
    try:
        # read dcraw binary output
        f = io.BytesIO(subprocess.check_output(['dcraw','-c','-o','0','-T','-4','-t','0','-D',raw_fn]))
        # create PIL image
        im = Image.open(f)
        # extract number array and return
        return array(im)

    except subprocess.CalledProcessError:
        raise Exception(u'ERROR while reading file with DCRAW')

def aql_exifjson(raw_fn,tags=[]):
    try:
        # get EXIFTOOL Output as JSON
        j = subprocess.check_output(['exiftool','-json','-G0','-n',raw_fn])
        exif_full = (json.loads(j))[0]
        # if tags array was given, filter them
        exif_filt = exif_full
        if len(tags) > 0:
            exif_filt = {}
            for t in tags:
                if t in exif_full:
                    exif_filt[t] = exif_full[t]
        # return filtered or full array
        return exif_filt
    except subprocess.CalledProcessError:
        raise Exception(u'ERROR while reading file with EXIFTOOL')

def aql_addash(arr):
    r = []
    for s in arr:
        r.append( '-'+s )
    return r

def aql_exif2hdr(exif):
    hdr = pyfits.PrimaryHDU().header

    hdr['DATE'] = datetime.now().strftime("%Y-%m-%dT%H:%M:%S")


    # add HISTORY field
    hdr['HISTORY'] = 'Imported with AQUILA script.'
    hdr['HISTORY'] = ' ...on ' + datetime.now().strftime("%Y/%m/%d %H:%M:%S")
    hdr['HISTORY'] = 'Dominik Gronkiewicz    gronki@gmail.com'

    # add dump of entire EXIF
    # hdr['COMMENT'] = 'Some EXIF tags are dumped below.'
    # for k,v in exif.iteritems():
    #     hdr['COMMENT'] = (str(k) + '=' + str(v),'EXIF DUMP')

    if u'EXIF:DateTimeOriginal' in exif:
        hdr['DATE-OBS'] = 0
    if 'EXIF:Model' in exif:
        hdr['INSTRUME'] = exif['EXIF:Model']
        aql_print('Camera: %s' % hdr['INSTRUME'])
    if 'EXIF:ExposureTime' in exif:
        hdr['EXPTIME'] = exif['EXIF:ExposureTime']
        aql_print( u'Exposure: %0.1f sec' % hdr['EXPTIME'] )
    if 'EXIF:ISO' in exif:
        hdr['ISO'] = exif['EXIF:ISO']
        aql_print( u'ISO %d' % hdr['ISO'] )
    if 'MakerNotes:CameraTemperature' in exif:
        hdr['CCD_TEMP'] = exif['MakerNotes:CameraTemperature']
        aql_print( u'Detector temperature: %0.1f C' % exif['MakerNotes:CameraTemperature'] )
    # if u'EXIF:BitsPerSample' in exif:
    #     hdr['NATIVBIT'] = exif[u'EXIF:BitsPerSample']
    #     aql_print( u'%d bits per sample' % exif[u'EXIF:BitsPerSample'] )
    return hdr


def aql_cfaflat(im,h):
    im2 = zeros(im.shape)
    im2[0::2,0::2] = im[0::2,0::2] / median(im[0::2,0::2])
    im2[0::2,1::2] = im[0::2,1::2] / median(im[0::2,1::2])
    im2[1::2,0::2] = im[1::2,0::2] / median(im[1::2,0::2])
    im2[1::2,1::2] = im[1::2,1::2] / median(im[1::2,1::2])
    return (im2,h)

def aql_cfa2rggb(im,h):
    im2 = zeros([ 4, im.shape[0] / 2, im.shape[1] / 2 ])
    im2[0,:,:] =  im[1::2,0::2]
    im2[1,:,:] =  im[0::2,0::2]
    im2[2,:,:] =  im[1::2,1::2]
    im2[3,:,:] =  im[0::2,1::2]
    return (im2,h)

def aql_rggb2grey(im_rggb,h,cr=0.0,cg1=0.5,cg2=0.5,cb=0.0):
    im_grey = zeros(im_rggb[0].shape)
    im_grey =   cr  * im_rggb[0] \
            +   cg1 * im_rggb[1] \
            +   cg2 * im_rggb[2] \
            +   cb  * im_rggb[3]
    a = 'R%0.2f-G%0.2f-G%0.2f-B%0.2f' % ( cr, cg1, cg2, cb )
    h['FILTER'] = a.replace('.','').replace('-','')
    return (im_grey,h)

def aql_raw2fits(raw_fn, fit_fn):


    # tags to be extracted from FITS
    important_tags = [  u'File:FileName', \
                u'File:FileType', \
                u'File:MIMEType', \
                u'EXIF:Width', \
                u'EXIF:ImageHeight', \
                u'EXIF:BitsPerSample', \
                u'EXIF:Compression', \
                u'EXIF:Make', \
                u'EXIF:Model', \
                u'EXIF:ExposureProgram', \
                u'EXIF:ISO', \
                u'EXIF:ExposureTime', \
                u'EXIF:DateTimeOriginal', \
                u'EXIF:CreateDate', \
                u'MakerNotes:CameraTemperature']


    if raw_fn == '':
        return


    # read RAW as number array
    imdata = aql_loadraw(raw_fn)

    # read EXIF as JSON
    exif = aql_exifjson(raw_fn,tags=important_tags)

    # generate image header from EXIF
    h = aql_exif2hdr(exif)

    h['CCDMEAN'] = mean(imdata)
    h['CCDMED'] = median(imdata)
    h['CCDSTD'] = std(imdata)

    return (imdata,h)
