#!/bin/bash
# modified from Bill Bristow to avoid anti-leeching measures
# reference: https://scivision.co/wget-download-an-http-directory-recursively/

date=$1
year=${date:0:4}
site=$2

case "$site" in
    PKR)
    type="DASC";;
    KAK)
    type="DASC";;
    TOO)
    type="CASC";;
esac

wget -nc -nd -nH -r -np --no-check-certificate -A.FITS -A.FIT \
 --random-wait --wait 1 -e robots=off \
 ftp://optics.gi.alaska.edu/$site/$type/RAW/$year/$date

