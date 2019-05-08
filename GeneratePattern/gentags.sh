#!/bin/bash

# Written by Leon Shaner
# Without modification, the defaults will generate 64 pages of 6x6 codes,
# 

# How many pages of output of rows x cols codes do you want?
maxcodes=64

# How many rows by cols do you want?
rows=6
cols=6

# This number has to be just big enough to cover your rows x cols x maxcodes
# NOTE: 18 bits is enough for 64 pages of 6x6 codes
bits=18

# How many crc bits do you want?
crcbits=6

# What scale do you want?
scale=1

# What filename prefix do you want?
#
# By default you'll get caltag18-6x6-0000.ps(mat) thru caltag18-6x6-NNNN.ps(mat)
# Which is computed from bits, rows, and cols, above, but you can change it here
# as needed.
fprefix="caltag${bits}-${rows}x${cols}-"



#########################################
# Nothing to customize below this line...
#########################################

size=$((rows * cols))
endcode=$((size * maxcodes))

# for padding filename with zeros
digits=$(echo $endcode | wc -c)
digits=$((digits -1))


#########################################
# Main Loop
#########################################

# Note:  padname takes numbers like 0 and pads them to numbers like 0000.
#        This is helpful for later steps, manipulating the .ps files,
#        because the filenames will be already in lexicographical order,
#        even just using shell filename wildcard (*) expansion.

currcode=0
while [ $currcode -lt $endcode ] ; do
    padname=$(echo $currcode | awk "{printf (\"%0${digits}d\", \$0)}")
    echo genrate_pattern.py --bits=$bits --crcbits=$crcbits --rows=$rows --cols=$cols --scale=$scale --offset=$currcode --file=${fprefix}${padname}
    ./genrate_pattern.py --bits=$bits --crcbits=$crcbits --rows=$rows --cols=$cols --scale=$scale --offset=$currcode --file=${fprefix}${padname}
    currcode=$((currcode + size))
done
