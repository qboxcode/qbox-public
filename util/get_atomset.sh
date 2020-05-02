#!/bin/bash
# get_atomset: extract atomset from a Qbox sample
#
# use: get_atomset sample.xml
#
nlines=$(grep /atomset -m 1 -n $1 | cut -f1 -d: - )
head -$nlines $1
echo "</fpmd:sample>"
