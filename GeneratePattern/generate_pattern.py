#!/usr/bin/env python
#
# This is a port of mda-newcaltag.C
# Note that template.ps must be in the cwd when executing this script

# 28/03/2012: [Ivo] Added white/black level and A4/A3/A0 paper sizes

from optparse import OptionParser
from math import sqrt
from time import localtime, strftime
from itertools import chain
import numpy as np
import scipy.io as sio
import string


class CRC:
    def __init__(self, idBits, crcBits):
        self.idBits = idBits
        self.crcBits = crcBits
        # from wikipedia
        genPolys = {1:0x1, 4:0x3, 5:0x15, 6:0x03, 7:0x09, 8:0x07,
                    10:0x233, 11:0x385, 12:0x80F, 15:0x4599,
                    16:0x1021, 24:0x5D6DCB, 32:0x04C11DB7}
        self.poly = genPolys[crcBits] | (2 ** crcBits)

    def encode(self, x):
        crc = x << self.crcBits
        div = self.poly << self.idBits
        while div >= self.poly:
            if (crc ^ div) < crc:
                crc ^= div
            div >>= 1
        return (x << self.crcBits) | crc

    def get_payload(self, x):
        return x >> self.crcBits

    def is_valid(self, x):
        return self.encode(self.get_payload(x)) == x


def countOnes_slow(x):
    ones = 0
    while True:
        if x & 1:
            ones += 1
        x = x >> 1
        if x == 0:
            break
    return ones
# precompute number of bits in all 8 bit uints
lut = [countOnes_slow(x) for x in range(256)]
def countOnes(x):
    ones = 0
    while x > 0:
        ones += lut[x & 0xFF]
        x >>= 8
    return ones


def hamming(x, y):
    return countOnes(x ^ y)


def numOnesOkay(x, bits):
    n = countOnes(x)
    return (n >= bits/4) and (n <= 3*bits/4)


def bin2dec(list):
    y = 0
    list.reverse()
    for i in range(len(list)):
        if list[i]:
            y += 2 ** i
    return y


#def rotateID(x, bits, k):
#    b = np.binary_repr(x, bits)
#    a = np.array(list(b), dtype=np.uint8)
#    n = int(sqrt(bits))
#    a = np.reshape(a, (n,n))
#    a = np.rot90(a, k)
#    a = a.flatten().tolist()
#    return bin2dec(a)


def rotateID90(x, bits):
    b = [x&(2**(i-1))>0 for i in range(bits,0,-1)]
    n = int(sqrt(bits))
    a = [b[i*n:(i+1)*n] for i in range(n)]
    a = zip(*a[::-1])
    b = list(chain(*[i for i in a]))
    b = [int(i) for i in b]
    s = [2**(len(b)-1-i)*x for (i,x) in enumerate(b)]
    return sum(s)


def markTooClose(id, codes, valid, minHamming):
    for idx, code in enumerate(codes):
        valid[idx] &= hamming(id, code) >= minHamming
    return valid


def computeCodes(idBits, crcBits, minHamming):
    crc = CRC(idBits, crcBits)
    numEntries = 2 ** idBits
    bits = idBits + crcBits
    codes = [crc.encode(i) for i in range(numEntries)]
    #print "initially, num codes =", len(codes)
    codes = filter(lambda x: numOnesOkay(x,idBits+crcBits), codes)
    #print "after bitcount filter =", len(codes)
    valid = [True for x in codes]
    for idx, x in enumerate(codes):
        if valid[idx]:
            for k in range(1,4):
                x = rotateID90(x, bits)
                valid = markTooClose(x, codes, valid, minHamming)
    comb = zip(codes, valid)
    codes = [x for (x,v) in comb if v]
    #print "after clash filter =", len(codes)
    ids = [crc.get_payload(x) for x in codes]
    return (codes, ids)


#Ivo: support different paper formats
#
#A4  -  8.3 x 11.7 in
#A3  - 11.7 x 16.5 in
#A0  - 33.1 x 46.8 in

def output(filename, markers, ids, nrows, ncols,
           blacklevel, whitelevel,
           scale, layout2,
           metric, minHamming, crcBits, idBits, width=8.5, height=11.0):
    if metric:
        scale /= 2.54
    landscape = nrows < ncols
    orientations = ["Portrait", "Landscape"]
    if landscape:
        comment = "% rotation to landscape mode"
        rotate = "%s\n%f 0 translate\n90 rotate\n" % (comment, width*72)
    else:
        rotate = ""
    codes = np.flipud(np.reshape(markers,(nrows,ncols)))
    ids = np.flipud(np.reshape(ids,(nrows,ncols)))
    markers = [str(i) for i in markers]
    markers = [markers[i:i+ncols] for i in range(0,len(markers),ncols)]
    markers = [string.join(i) for i in markers]
    markers = string.join(markers, "\n")
    template = open("template.ps", "r").read()
    template = string.Template(template)
    s = template.substitute(
            PAPERWIDTH_TIMES_72 = width*72,
            PAPERHEIGHT_TIMES_72 = height*72,
            ORIENTATION = orientations[landscape],
            CREATION_DATE = strftime("%Y-%m-%d %H:%M:%S", localtime()),
            LAYOUT = int(not layout2),
            MARKER_SIZE_IN_INCHES = scale,
            NUM_COLUMNS = ncols,
            NUM_ROWS = nrows,
            BLACK_LEVEL = blacklevel,
            WHITE_LEVEL = whitelevel,
            MARKER_IDS = markers,
            LANDSCAPE_ROTATION = rotate)
    dest = open(filename+".ps", "w")
    with open(filename+".ps", "w") as f:
        f.write(s)
    resCode = np.array([4.,4.])
    resMarker = np.array([8.,8.])
    resPattern = np.array([nrows,ncols], dtype=np.double)
    codes = np.array(codes, dtype=np.double)
    ids = np.array(ids, dtype=np.double)
    if metric:
        scale *= 2.54
    mDict = {'crcBits':crcBits, 'idBits':idBits, 'layout':int(layout2)+1,
             'minHamming':minHamming, 'ID':ids, 'CODE':codes,
             'resCode':resCode, 'resMarker':resMarker, 'resPattern':resPattern,
             'scale':scale}
    sio.savemat(filename+".mat", mDict, oned_as='row')




if __name__ == "__main__":
    p = OptionParser()
    p.add_option("--bits", dest="bits", type="int", default=16,
                 help="number of marker bits (must be square)")
    p.add_option("--crcbits", dest="crcBits", type="int", default=6,
                 help="number of checksum bits (must be < bits)")
    p.add_option("-r", "--rows", dest="rows", type="int", default=6,
                 help="number of rows")
    p.add_option("-c", "--cols", dest="cols", type="int", default=4,
                 help="number of columns")
    p.add_option("-b","--blacklevel", type="float",default=0.0,
                 help="black level of checkerboard (between 0 and 1)")
    p.add_option("-w","--whitelevel", type="float",default=1.0,
                 help="white level of checkerboard (between 0 and 1)")
    p.add_option("-s", "--scale", dest="scale", type="float", default=1.0,
                 help="marker size in inches")
    p.add_option("-m", "--minhamming", dest="minHamming", type="int", default=2,
                 help="minimum hamming distance between codes")
    p.add_option("-o", "--offset", dest="offset", type="int", default=0,
                 help="use codes starting from this index")
    p.add_option("--metric", dest="metric", action="store_true",
                 default=False, help="use centimeters instead of inches")
    p.add_option("--layout2", dest="layout2", action="store_true",
                 default=False, help="use rotated layout")
    p.add_option("-f", "--file", dest="filename", default="output",
                 help="output filenames will be FILENAME.ps and .mat")
    (opt, args) = p.parse_args()
    opt.idBits = opt.bits - opt.crcBits


    codes, ids = computeCodes(opt.idBits, opt.crcBits, opt.minHamming)
    nMarkers = opt.rows * opt.cols
    if nMarkers > len(codes):
        print "Error: insufficient marker IDs for requested resolution"
        exit()
    a, b = opt.offset, opt.offset + nMarkers
    if b > len(codes):
        print "Error: insufficient markers IDs from offset {}".format(a)
        exit()
    print( "Using codes {}...{} of {}...{}".format(a,b-1,0,len(codes)-1) )
    markers = codes[a:b]
    ids = ids[a:b]
    output(opt.filename, markers, ids, opt.rows, opt.cols,
           opt.blacklevel, opt.whitelevel,
           opt.scale,
           opt.layout2, opt.metric, opt.minHamming, opt.crcBits,
           opt.bits-opt.crcBits)

    print "Finished"

